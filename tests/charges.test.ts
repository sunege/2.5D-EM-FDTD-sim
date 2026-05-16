import { describe, expect, it } from "vitest";
import { depositCurrent, makeCharge, V_MAX } from "../src/sim/charges";
import { gaussianDensity, radialField } from "../src/sim/coulomb";
import { createFields, idx } from "../src/sim/fields";
import type { SimParams } from "../src/sim/params";

const p: SimParams = {
  Nx: 48, Ny: 48, Nz: 21,
  dx: 1, dy: 1, dz: 1, c: 1, dt: 0.4,
};

describe("smoothstep injection", () => {
  it("per-step weights sum to 1", () => {
    // Inject into a fresh field repeatedly and verify the accumulated Ex at
    // a probe point matches the analytic Coulomb field after injectTotal steps.
    const f = createFields(p);
    const cx = 24, cy = 24;
    const zMid = (p.Nz - 1) / 2;
    const charges = [makeCharge(cx, cy, 1, 2, 20)];
    for (let n = 0; n < 20; n++) depositCurrent(f, p, charges, zMid);

    // Probe Ex at (cx+0.5+3, cy, zMid): position (cx+3.5, cy, zMid).
    const i = cx + 3, j = cy, k = zMid;
    const got = f.Ex[idx(i, j, k, p.Nx, p.Ny)];
    const dx = i + 0.5 - cx;
    const r = dx; // dy=dz=0
    const expected = radialField(r, 2, 1) * (dx / r);
    expect(got).toBeCloseTo(expected, 6);
  });

  it("first-step weight is small (smooth start, not boxcar)", () => {
    const f = createFields(p);
    const zMid = (p.Nz - 1) / 2;
    const charges = [makeCharge(24, 24, 1, 2, 20)];
    depositCurrent(f, p, charges, zMid);
    const i = 27, j = 24, k = zMid;
    const after1 = f.Ex[idx(i, j, k, p.Nx, p.Ny)];

    // Reset and run all 20 steps; the first-step contribution should be a
    // small fraction of the total (much less than the 1/20 = 5% of the old
    // boxcar profile, because quintic-smoothstep starts flat).
    const f2 = createFields(p);
    const charges2 = [makeCharge(24, 24, 1, 2, 20)];
    for (let n = 0; n < 20; n++) depositCurrent(f2, p, charges2, zMid);
    const final = f2.Ex[idx(i, j, k, p.Nx, p.Ny)];

    const ratio = after1 / final;
    expect(ratio).toBeGreaterThan(0);
    expect(ratio).toBeLessThan(0.005); // boxcar gave 1/20 = 0.05; we want ≪ that
  });

  it("stops depositing after injectTotal steps", () => {
    const f = createFields(p);
    const zMid = (p.Nz - 1) / 2;
    const charges = [makeCharge(24, 24, 1, 2, 20)];
    for (let n = 0; n < 20; n++) depositCurrent(f, p, charges, zMid);
    const snapshot = new Float32Array(f.Ex);
    for (let n = 0; n < 10; n++) depositCurrent(f, p, charges, zMid);
    for (let i = 0; i < snapshot.length; i++) {
      expect(f.Ex[i]).toBe(snapshot[i]);
    }
  });
});

describe("motion deposit (J = ρ·v)", () => {
  it("does nothing when target equals current position", () => {
    const f = createFields(p);
    const zMid = (p.Nz - 1) / 2;
    const c = makeCharge(24, 24, 1, 2, 5);
    const charges = [c];
    for (let n = 0; n < 5; n++) depositCurrent(f, p, charges, zMid); // finish injection
    const snapshot = new Float32Array(f.Ex);
    // target is still (24, 24) after makeCharge, so motion pass is a no-op.
    for (let n = 0; n < 10; n++) depositCurrent(f, p, charges, zMid);
    for (let i = 0; i < snapshot.length; i++) {
      expect(f.Ex[i]).toBe(snapshot[i]);
    }
    expect(c.vx).toBe(0);
    expect(c.vy).toBe(0);
  });

  it("advances position toward the target at sub-V_MAX speed", () => {
    const f = createFields(p);
    const zMid = (p.Nz - 1) / 2;
    const c = makeCharge(24, 24, 1, 2, 5);
    const charges = [c];
    for (let n = 0; n < 5; n++) depositCurrent(f, p, charges, zMid);
    // Slow steady drag: per-step desired displacement well below V_MAX·dt = 0.2.
    const step = 0.05;
    for (let n = 0; n < 10; n++) {
      c.targetX = c.x + step;
      depositCurrent(f, p, charges, zMid);
    }
    // After 10 steps of step=0.05, charge should track tightly.
    expect(c.x).toBeCloseTo(24 + 10 * step, 6);
    expect(c.vx).toBeCloseTo(step / p.dt, 6);
  });

  it("clips per-step displacement to V_MAX · dt", () => {
    const f = createFields(p);
    const zMid = (p.Nz - 1) / 2;
    const c = makeCharge(24, 24, 1, 2, 5);
    const charges = [c];
    for (let n = 0; n < 5; n++) depositCurrent(f, p, charges, zMid);
    c.targetX = 100; // far beyond V_MAX·dt
    const before = c.x;
    depositCurrent(f, p, charges, zMid);
    const stepDist = c.x - before;
    expect(stepDist).toBeCloseTo(V_MAX * p.dt, 6);
    expect(Math.hypot(c.vx, c.vy)).toBeCloseTo(V_MAX, 6);
  });

  it("motion carries the taper shell along, leaving no residual at the old position", () => {
    // After dragging far enough that the new charge's support no longer
    // overlaps the old taper shell, the discrete divergence at points in
    // the *old* taper shell must collapse to near zero. With the bare
    // ρ_gauss motion deposit this would leave the −q shell behind as a
    // residual of order |E_coulomb(4σ)·w'(4σ)| ≈ 5e-4.
    const zMid = (p.Nz - 1) / 2;
    const kc = Math.round(zMid);
    const f = createFields(p);
    const c = makeCharge(24, 24, 1, 2, 10);
    for (let n = 0; n < 10; n++) depositCurrent(f, p, [c], zMid);
    // Slow drag (v = 0.25, no clipping). 100 steps × 0.1 = 10 cells = 5σ.
    for (let n = 0; n < 100; n++) {
      c.targetX = c.x + 0.1;
      depositCurrent(f, p, [c], zMid);
    }
    expect(c.x).toBeCloseTo(34, 6);
    const divAt = (i: number, j: number, k: number): number =>
      f.Ex[idx(i, j, k, p.Nx, p.Ny)] - f.Ex[idx(i - 1, j, k, p.Nx, p.Ny)]
      + f.Ey[idx(i, j, k, p.Nx, p.Ny)] - f.Ey[idx(i, j - 1, k, p.Nx, p.Ny)]
      + f.Ez[idx(i, j, k, p.Nx, p.Ny)] - f.Ez[idx(i, j, k - 1, p.Nx, p.Ny)];
    // Probe at (16, 24, 10): r=8 from old (in old taper shell at 4σ),
    // r=18 from new (well beyond the new charge's 5σ support). Without
    // the shell-aware motion deposit, divergence here is ≈ -5e-4.
    const divResidual = divAt(16, 24, kc);
    expect(Math.abs(divResidual)).toBeLessThan(1e-4);
  });

  it("motion deposit translates ρ_eff = ∇·E (continuity check)", () => {
    // The motion deposit J = ρ·v guarantees ∇·J = −∂ρ/∂t, so the discrete
    // divergence of E (which equals the implicit charge density ρ_eff)
    // should track the charge's position. After dragging by 3 cells, the
    // peak of ρ_eff must move with it. We don't expect E itself to look
    // like a translated Coulomb field — the motion deposit also writes
    // curl, which would radiate away under full FDTD evolution but
    // persists in this isolated test.
    const zMid = (p.Nz - 1) / 2;
    const kc = Math.round(zMid);
    const f = createFields(p);
    const c = makeCharge(24, 24, 1, 2, 10);
    for (let n = 0; n < 10; n++) depositCurrent(f, p, [c], zMid);
    for (let n = 0; n < 30; n++) {
      c.targetX = c.x + 0.1;
      depositCurrent(f, p, [c], zMid);
    }
    expect(c.x).toBeCloseTo(27, 6);
    const divAt = (i: number, j: number, k: number): number =>
      f.Ex[idx(i, j, k, p.Nx, p.Ny)] - f.Ex[idx(i - 1, j, k, p.Nx, p.Ny)]
      + f.Ey[idx(i, j, k, p.Nx, p.Ny)] - f.Ey[idx(i, j - 1, k, p.Nx, p.Ny)]
      + f.Ez[idx(i, j, k, p.Nx, p.Ny)] - f.Ez[idx(i, j, k - 1, p.Nx, p.Ny)];
    const rhoMax = gaussianDensity(0, 2, 1);
    // New center now hosts the ρ peak — finite-difference discretization
    // gives ~5–15% error vs. analytic ρ_max, so use a band of 15%.
    const divNew = divAt(27, 24, kc);
    expect(divNew).toBeGreaterThan(0.85 * rhoMax);
    expect(divNew).toBeLessThan(1.15 * rhoMax);
    // Old center: ρ has moved away. Expected ρ_gauss(3) / ρ_max ≈ 0.105.
    const divOld = divAt(24, 24, kc);
    expect(divOld).toBeGreaterThan(0); // still inside the Gaussian's tail
    expect(divOld).toBeLessThan(0.2 * rhoMax);
  });
});
