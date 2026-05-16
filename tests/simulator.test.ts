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

  it("places a charge at the z mid-plane and generates a radial field", () => {
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
    // Run through the injection window plus a few steps for the field to settle.
    for (let n = 0; n < 30; n++) s.step();
    const cx = 32, cy = 32, kMid = 10;
    // Probe Ex a few cells to the +x side of the charge — should be positive (outward).
    const ex = s.fields.Ex[idx(cx + 4, cy, kMid, s.params.Nx, s.params.Ny)];
    expect(ex).toBeGreaterThan(0);
    // Energy after injection should be finite and positive.
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
