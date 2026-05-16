import { gaussianDensity, radialField } from "./coulomb";
import type { Fields } from "./fields";
import { idx } from "./fields";
import type { SimParams } from "./params";

export interface Charge {
  // Continuous position in grid units. zMid (== (Nz-1)/2) is provided
  // externally so the same Charge can be used regardless of grid changes.
  x: number;
  y: number;
  // Current velocity in grid units per simulation-time. Computed each step
  // from (target − position)/dt and clipped to V_MAX.
  vx: number;
  vy: number;
  // Drag target. When equal to (x, y) the charge is at rest.
  targetX: number;
  targetY: number;
  q: number;
  sigma: number;
  injectStep: number;
  injectTotal: number;
}

// Cap on drag speed in units of c (= 1 in normalized units). The wave
// equation cannot propagate information faster than c, so superluminal
// motion would violate causality and create numerical artifacts. A cap
// below c also bounds the per-step displacement (≤ V_MAX · dt grid cells),
// which keeps the motion deposit well-resolved on the Yee grid.
export const V_MAX = 0.5;

// Taper geometry — must match the injection's taper exactly. The motion
// deposit reuses these so that the same ρ_eff is being translated.
const TAPER_INNER_FACTOR = 3;
const TAPER_OUTER_FACTOR = 5;

// Effective volume charge density implicitly set up by the tapered
// injection: ρ_eff = ∇·(E_coulomb · w(r) · r̂) = ρ_gauss·w + E_coulomb·w'.
// In the taper region [3σ, 5σ] w' < 0, so ρ_eff has a negative "shell"
// term carrying ~−q. Total integral is 0 because the truncated field has
// zero flux at infinity, i.e. the injected configuration has zero net
// charge — the visible Gaussian core +q is balanced by the shell −q.
//
// Motion uses *this* density (not bare ρ_gauss) so that J = ρ_eff·v
// translates the entire injected configuration — core *and* shell —
// preserving continuity for the full ρ_eff. With bare ρ_gauss only the
// core moves and the −q shell is left at the old center, visible from
// outside as a Coulomb-like residual field pointing inward.
function effectiveInjectedDensity(r: number, sigma: number, q: number): number {
  const taperInner = TAPER_INNER_FACTOR * sigma;
  const taperOuter = TAPER_OUTER_FACTOR * sigma;
  if (r >= taperOuter) return 0;
  const rho = gaussianDensity(r, sigma, q);
  if (r <= taperInner) return rho;
  const taperWidth = taperOuter - taperInner;
  const t = (r - taperInner) / taperWidth;
  const w = 0.5 * (1 + Math.cos(Math.PI * t));
  const dwdr = -(Math.PI / (2 * taperWidth)) * Math.sin(Math.PI * t);
  return rho * w + radialField(r, sigma, q) * dwdr;
}

export function makeCharge(
  x: number,
  y: number,
  q: number,
  sigma = 2,
  injectTotal = 20,
): Charge {
  return {
    x,
    y,
    vx: 0,
    vy: 0,
    targetX: x,
    targetY: y,
    q,
    sigma,
    injectStep: 0,
    injectTotal,
  };
}

// During injection, we deposit fractional Coulomb-field increments such that
// after injectTotal steps the cumulative addition equals the analytic Coulomb
// field of the Gaussian charge. Per-step weight w_n = S(t_{n+1}) − S(t_n)
// where S is a quintic smoothstep (S(0)=0, S(1)=1, S' and S'' vanish at both
// endpoints). This eliminates the boxcar discontinuity of constant-weight
// injection, dropping the radiated high-frequency content by ~1/ω⁶ and
// suppressing the residual ripple that fills the grid after a click.
// Continuity ∇·J = −∂ρ/∂t is still satisfied exactly because ∇·E_Coulomb = ρ.
const smoothstep5 = (t: number): number => t * t * t * (t * (t * 6 - 15) + 10);

export function depositCurrent(
  f: Fields,
  p: SimParams,
  charges: Charge[],
  zMid: number,
): void {
  const { Nx, Ny, Nz, dt } = p;

  for (const c of charges) {
    // ============ Motion pass: J = ρ·v ============
    // For a translating Gaussian ρ(x, t) = ρ_0(x − x_c(t)),
    //   ∂ρ/∂t = −v · ∇ρ_0,  J = ρ_0 · v  ⇒  ∇·J = v · ∇ρ_0,
    // so continuity ∂ρ/∂t + ∇·J = 0 holds exactly when v is spatially
    // constant. Depositing −dt·J into E at each Yee position therefore
    // moves the implicit ρ (= ∇·E) rigidly without breaking Gauss's law.
    // Radiation appears only when v changes (acceleration / deceleration).
    const dxT = c.targetX - c.x;
    const dyT = c.targetY - c.y;
    const dist = Math.hypot(dxT, dyT);
    if (dist > 1e-10) {
      // Aim to reach target in one step, then clip to subluminal speed.
      let vx = dxT / dt;
      let vy = dyT / dt;
      const vMag = Math.hypot(vx, vy);
      if (vMag > V_MAX) {
        const s = V_MAX / vMag;
        vx *= s;
        vy *= s;
      }
      c.vx = vx;
      c.vy = vy;

      // Evaluate ρ at the midpoint position x_c + v·dt/2 (second-order
      // accurate in dt). The window is the same 5σ sphere as injection,
      // so density at the edge is exp(-25) ≈ 1.4e-11 — negligible.
      const cxMid = c.x + 0.5 * vx * dt;
      const cyMid = c.y + 0.5 * vy * dt;
      const halfW = 5 * c.sigma;
      const r2max = halfW * halfW;

      const iMin = Math.max(0, Math.floor(cxMid - halfW));
      const iMax = Math.min(Nx - 1, Math.ceil(cxMid + halfW));
      const jMin = Math.max(0, Math.floor(cyMid - halfW));
      const jMax = Math.min(Ny - 1, Math.ceil(cyMid + halfW));
      const kMin = Math.max(0, Math.floor(zMid - halfW));
      const kMax = Math.min(Nz - 1, Math.ceil(zMid + halfW));

      const negDtVx = -dt * vx;
      const negDtVy = -dt * vy;

      for (let k = kMin; k <= kMax; k++) {
        for (let j = jMin; j <= jMax; j++) {
          for (let i = iMin; i <= iMax; i++) {
            const id = idx(i, j, k, Nx, Ny);
            // Ex at (i+0.5, j, k): −dt · ρ_eff(r) · vx
            {
              const dx = i + 0.5 - cxMid;
              const dy = j - cyMid;
              const dz = k - zMid;
              const r2 = dx * dx + dy * dy + dz * dz;
              if (r2 <= r2max) {
                const rho = effectiveInjectedDensity(Math.sqrt(r2), c.sigma, c.q);
                f.Ex[id] += negDtVx * rho;
              }
            }
            // Ey at (i, j+0.5, k): −dt · ρ_eff(r) · vy
            {
              const dx = i - cxMid;
              const dy = j + 0.5 - cyMid;
              const dz = k - zMid;
              const r2 = dx * dx + dy * dy + dz * dz;
              if (r2 <= r2max) {
                const rho = effectiveInjectedDensity(Math.sqrt(r2), c.sigma, c.q);
                f.Ey[id] += negDtVy * rho;
              }
            }
            // Ez: vz = 0, no contribution.
          }
        }
      }

      // Advance position by v·dt.
      c.x += vx * dt;
      c.y += vy * dt;
    } else {
      c.vx = 0;
      c.vy = 0;
    }

    // ============ Injection pass (Coulomb field build-up) ============
    if (c.injectStep >= c.injectTotal) continue;
    const t0 = c.injectStep / c.injectTotal;
    const t1 = (c.injectStep + 1) / c.injectTotal;
    const factor = smoothstep5(t1) - smoothstep5(t0); // ∑ over all steps = 1
    // Spherical deposit region with a cosine taper near the outer surface.
    // A hard cutoff at r = halfW leaves an O(q/(4πhalfW²)) step in E_r at
    // the boundary; that step radiates broadband ripples via the discrete
    // Yee curl and leaves a non-radial O(dx²) static residue ("stuck"
    // inward arrows). Tapering w(r) from 1 at r ≤ taperInner to 0 at
    // r ≥ halfW spreads the implicit shell charge over a thickness ~σ,
    // band-limiting its spectrum and erasing the residue. Continuity
    // ∇·J = −∂ρ_eff/∂t still holds because ρ_eff is implicitly defined as
    // ∇·(E_Coulomb · w) and J is its time derivative.
    const halfW = 5 * c.sigma;
    const taperInner = 3 * c.sigma;
    const taperWidth = halfW - taperInner; // = 2 · c.sigma
    const r2max = halfW * halfW;
    const r2taper = taperInner * taperInner;

    const taperWeight = (r2: number, r: number): number => {
      if (r2 <= r2taper) return 1;
      const t = (r - taperInner) / taperWidth;
      return 0.5 * (1 + Math.cos(Math.PI * t));
    };

    const iMin = Math.max(0, Math.floor(c.x - halfW));
    const iMax = Math.min(Nx - 1, Math.ceil(c.x + halfW));
    const jMin = Math.max(0, Math.floor(c.y - halfW));
    const jMax = Math.min(Ny - 1, Math.ceil(c.y + halfW));
    const kMin = Math.max(0, Math.floor(zMid - halfW));
    const kMax = Math.min(Nz - 1, Math.ceil(zMid + halfW));

    for (let k = kMin; k <= kMax; k++) {
      for (let j = jMin; j <= jMax; j++) {
        for (let i = iMin; i <= iMax; i++) {
          const id = idx(i, j, k, Nx, Ny);
          // Ex at (i+0.5, j, k)
          const dxX = i + 0.5 - c.x;
          const dyX = j - c.y;
          const dzX = k - zMid;
          const r2X = dxX * dxX + dyX * dyX + dzX * dzX;
          if (r2X > 0 && r2X <= r2max) {
            const rX = Math.sqrt(r2X);
            const Er = radialField(rX, c.sigma, c.q);
            f.Ex[id] += factor * Er * (dxX / rX) * taperWeight(r2X, rX);
          }
          // Ey at (i, j+0.5, k)
          const dxY = i - c.x;
          const dyY = j + 0.5 - c.y;
          const dzY = k - zMid;
          const r2Y = dxY * dxY + dyY * dyY + dzY * dzY;
          if (r2Y > 0 && r2Y <= r2max) {
            const rY = Math.sqrt(r2Y);
            const Er = radialField(rY, c.sigma, c.q);
            f.Ey[id] += factor * Er * (dyY / rY) * taperWeight(r2Y, rY);
          }
          // Ez at (i, j, k+0.5)
          const dxZ = i - c.x;
          const dyZ = j - c.y;
          const dzZ = k + 0.5 - zMid;
          const r2Z = dxZ * dxZ + dyZ * dyZ + dzZ * dzZ;
          if (r2Z > 0 && r2Z <= r2max) {
            const rZ = Math.sqrt(r2Z);
            const Er = radialField(rZ, c.sigma, c.q);
            f.Ez[id] += factor * Er * (dzZ / rZ) * taperWeight(r2Z, rZ);
          }
        }
      }
    }
    c.injectStep++;
  }
}
