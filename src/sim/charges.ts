import type { Fields } from "./fields";
import { idx } from "./fields";
import type { SimParams } from "./params";

// A click creates a dipole. A wire current J_x flows from A=(xA,yA) to
// B=(xB,yB), and continuity ∂ρ/∂t = −∇·J deposits +q at B and −q at A.
// After injection, each end can be dragged independently; the motion current
// J = ρ·v at the moving end radiates an EM wave while the static end stays.
export interface DipoleCharge {
  xA: number; yA: number; vxA: number; vyA: number;  // −q end
  xB: number; yB: number; vxB: number; vyB: number;  // +q end
  q: number;
  sigma: number;
  injectStep: number;
  injectTotal: number;
}

// Default half-separation in units of sigma. 2σ keeps the two charge discs
// visually distinct while the dipole still fits comfortably in the domain.
const DEFAULT_HALF_SEP_FACTOR = 2;

export function makeDipole(
  cx: number,
  cy: number,
  q: number,
  sigma = 2,
  injectTotal = 10,
): DipoleCharge {
  const d = DEFAULT_HALF_SEP_FACTOR * sigma;
  return {
    xA: cx - d, yA: cy, vxA: 0, vyA: 0,
    xB: cx + d, yB: cy, vxB: 0, vyB: 0,
    q, sigma, injectStep: 0, injectTotal,
  };
}

// Per-step injection weight w_n = S(t_{n+1}) − S(t_n) where S is a quintic
// smoothstep — S′ and S″ vanish at both endpoints, so the radiated spectrum
// from the current pulse decays as ω⁻⁸ instead of ω⁻² (boxcar).
const smoothstep5 = (t: number): number => t * t * t * (t * (t * 6 - 15) + 10);

// Transverse Gaussian support. exp(−16) ≈ 1.1e−7 — negligible truncation.
const TRANSVERSE_HALF_W_FACTOR = 4;

export function depositCurrent(
  f: Fields,
  p: SimParams,
  charges: DipoleCharge[],
  zMid: number,
): void {
  const { Nx, Ny, Nz } = p;

  for (const c of charges) {
    if (c.injectStep >= c.injectTotal) continue;

    const t0 = c.injectStep / c.injectTotal;
    const t1 = (c.injectStep + 1) / c.injectTotal;
    const factor = smoothstep5(t1) - smoothstep5(t0);

    const ax = c.xA;
    const bx = c.xB;
    const wireCy = c.yA;  // wire y-center (yA == yB on new placement)
    const sigma2 = c.sigma * c.sigma;
    const halfW = TRANSVERSE_HALF_W_FACTOR * c.sigma;

    // Ex sits at (i+½, j, k); wire is active when ax ≤ i+½ ≤ bx.
    const iMin = Math.max(0, Math.ceil(ax - 0.5));
    const iMax = Math.min(Nx - 1, Math.floor(bx - 0.5));
    if (iMin > iMax) continue;

    const jMin = Math.max(0, Math.floor(wireCy - halfW));
    const jMax = Math.min(Ny - 1, Math.ceil(wireCy + halfW));
    const kMin = Math.max(0, Math.floor(zMid - halfW));
    const kMax = Math.min(Nz - 1, Math.ceil(zMid + halfW));

    // Transverse normalization: discrete ∇·(ΔE) at the B-end cell sums to
    // exactly q over (j, k) — and −q at the A-end — after all injectTotal
    // steps, so charge conservation is exact to machine precision.
    let aEff = 0;
    for (let k = kMin; k <= kMax; k++) {
      const dz = k - zMid;
      for (let j = jMin; j <= jMax; j++) {
        const dy = j - wireCy;
        aEff += Math.exp(-(dy * dy + dz * dz) / sigma2);
      }
    }
    if (aEff <= 0) continue;

    // Ampère: ∂E/∂t = ∇×B − J. Injecting J_x > 0 (A → B) means subtracting
    // from Ex. Only at the two wire endpoints does ΔEx jump, producing discrete
    // ρ = ±q·G_yz/aEff.
    const coeff = factor * c.q / aEff;
    for (let k = kMin; k <= kMax; k++) {
      const dz = k - zMid;
      for (let j = jMin; j <= jMax; j++) {
        const dy = j - wireCy;
        const g = Math.exp(-(dy * dy + dz * dz) / sigma2);
        const delta = coeff * g;
        for (let i = iMin; i <= iMax; i++) {
          f.Ex[idx(i, j, k, Nx, Ny)] -= delta;
        }
      }
    }

    c.injectStep++;
  }
}

// Advance endpoint positions by v·dt. Must be called AFTER depositMotion so
// that J is deposited at the position where the charge is (before the move),
// satisfying continuity ∂ρ/∂t = −∇·J for each step.
export function advancePositions(charges: DipoleCharge[], dt: number): void {
  for (const c of charges) {
    c.xA += c.vxA * dt;
    c.yA += c.vyA * dt;
    c.xB += c.vxB * dt;
    c.yB += c.vyB * dt;
  }
}

// Deposit J = q·ρ(x)·v for each endpoint that is moving. Unlike wire
// injection (Ex only), motion can be in any direction → both Ex and Ey.
export function depositMotion(
  f: Fields,
  p: SimParams,
  charges: DipoleCharge[],
  zMid: number,
): void {
  const { Nx, Ny, Nz, dt } = p;
  for (const c of charges) {
    if (c.vxB !== 0 || c.vyB !== 0)
      _depositEndpoint(f, Nx, Ny, Nz, dt, c.xB, c.yB, zMid, c.sigma, +c.q, c.vxB, c.vyB);
    if (c.vxA !== 0 || c.vyA !== 0)
      _depositEndpoint(f, Nx, Ny, Nz, dt, c.xA, c.yA, zMid, c.sigma, -c.q, c.vxA, c.vyA);
  }
}

function _depositEndpoint(
  f: Fields,
  Nx: number,
  Ny: number,
  Nz: number,
  dt: number,
  xE: number,
  yE: number,
  zMid: number,
  sigma: number,
  q: number,
  vx: number,
  vy: number,
): void {
  const sigma2 = sigma * sigma;
  const halfW = 4 * sigma;
  const iMin = Math.max(0, Math.floor(xE - halfW));
  const iMax = Math.min(Nx - 1, Math.ceil(xE + halfW));
  const jMin = Math.max(0, Math.floor(yE - halfW));
  const jMax = Math.min(Ny - 1, Math.ceil(yE + halfW));
  const kMin = Math.max(0, Math.floor(zMid - halfW));
  const kMax = Math.min(Nz - 1, Math.ceil(zMid + halfW));

  // 3D normalization using integer-cell positions. The half-cell shift
  // introduces < 1% error for sigma ≥ 2 — acceptable for the demo.
  let vEff = 0;
  for (let k = kMin; k <= kMax; k++) {
    const dz = k - zMid;
    for (let j = jMin; j <= jMax; j++) {
      const dy = j - yE;
      for (let i = iMin; i <= iMax; i++) {
        const dx = i - xE;
        vEff += Math.exp(-(dx * dx + dy * dy + dz * dz) / sigma2);
      }
    }
  }
  if (vEff <= 0) return;

  const cX = dt * q * vx / vEff;
  const cY = dt * q * vy / vEff;
  for (let k = kMin; k <= kMax; k++) {
    const dz = k - zMid;
    for (let j = jMin; j <= jMax; j++) {
      const dy = j - yE;
      for (let i = iMin; i <= iMax; i++) {
        const dx = i - xE;
        const g = Math.exp(-(dx * dx + dy * dy + dz * dz) / sigma2);
        f.Ex[idx(i, j, k, Nx, Ny)] -= cX * g;
        f.Ey[idx(i, j, k, Nx, Ny)] -= cY * g;
      }
    }
  }
}
