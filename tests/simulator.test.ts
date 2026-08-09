import { describe, expect, it } from "vitest";
import { idx } from "../src/sim/fields";
import { Simulator } from "../src/sim/simulator";

describe("Simulator", () => {
  it("starts with zero fields and no charges", () => {
    const s = new Simulator();
    expect(s.charges.length).toBe(0);
    expect(s.energy()).toBe(0);
    expect(s.stepCount).toBe(0);
  });

  it("places a dipole that produces outward Ex on the +q side", () => {
    const s = new Simulator({
      Nx: 64,
      Ny: 64,
      Nz: 21,
      dx: 1,
      dy: 1,
      dz: 1,
      c: 1,
      dt: 0.4,
    });
    s.addCharge(32, 32, 1);
    // Default sigma=2 → d=4, so B (+q end) sits at x=36 and A (−q) at x=28.
    // Run through injection + propagation so the near-field is established.
    for (let n = 0; n < 30; n++) s.step();
    const kMid = 10;
    // Probe just outside B (Yee position 36.5, ≈0.5 cells past the +q end).
    // The +q dominates over the far −q, so Ex must point outward (+x).
    const ex = s.fields.Ex[idx(36, 32, kMid, s.params.Nx, s.params.Ny)];
    expect(ex).toBeGreaterThan(0);
    expect(s.energy()).toBeGreaterThan(0);
    expect(Number.isFinite(s.energy())).toBe(true);
  });

  it("does not blow up after many steps with a single charge", () => {
    const s = new Simulator({
      Nx: 48,
      Ny: 48,
      Nz: 21,
      dx: 1,
      dy: 1,
      dz: 1,
      c: 1,
      dt: 0.4,
    });
    s.addCharge(24, 24, 1);
    for (let n = 0; n < 500; n++) s.step();
    expect(Number.isFinite(s.energy())).toBe(true);
    // Energy should be modest; the radiation pulse has left and only the
    // residual Coulomb field remains.
    expect(s.energy()).toBeLessThan(1000);
  });
});
