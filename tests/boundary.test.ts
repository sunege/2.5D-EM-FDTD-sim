import { describe, expect, it } from "vitest";
import { applyMurAbc, createMurState, snapshotBoundary } from "../src/sim/boundary";
import { createFields, idx } from "../src/sim/fields";
import { totalEnergy, updateB, updateE } from "../src/sim/fdtd";
import type { SimParams } from "../src/sim/params";

function makeRadialPulse(N: number, dt = 0.4) {
  const p: SimParams = { Nx: N, Ny: N, Nz: 11, dx: 1, dy: 1, dz: 1, c: 1, dt };
  const f = createFields(p);
  const cx = N / 2, cy = N / 2, cz = (p.Nz - 1) / 2;
  const w = 4;
  for (let k = 0; k < p.Nz; k++) {
    for (let j = 0; j < p.Ny; j++) {
      for (let i = 0; i < p.Nx; i++) {
        const r2 = (i - cx) ** 2 + (j - cy) ** 2 + (k - cz) ** 2;
        f.Ez[idx(i, j, k, p.Nx, p.Ny)] = Math.exp(-r2 / (w * w));
      }
    }
  }
  return { p, f };
}

describe("Mur 1st-order ABC", () => {
  it("absorbs an outgoing radial pulse (energy drops below 10% of peak)", () => {
    const { p, f } = makeRadialPulse(40);
    const m = createMurState(p);
    const e0 = totalEnergy(f);
    let eMax = e0;
    // Run long enough for the entire pulse to leave the grid.
    // Wave travel time across radius (20 cells) at c*dt=0.4 ≈ 50 steps.
    const steps = 200;
    for (let n = 0; n < steps; n++) {
      snapshotBoundary(f, p, m);
      updateB(f, p);
      updateE(f, p);
      applyMurAbc(f, p, m);
      const e = totalEnergy(f);
      if (e > eMax) eMax = e;
    }
    const eFinal = totalEnergy(f);
    expect(eFinal / eMax).toBeLessThan(0.1);
  });

  it("preserves fields when r=0 (c*dt = dx exact magic angle)", () => {
    // When c*dt=dx, r=0 and Mur 1st order is exact for normal incidence.
    const p: SimParams = { Nx: 24, Ny: 24, Nz: 11, dx: 1, dy: 1, dz: 1, c: 1, dt: 1 };
    // Note: dt=1 violates 3D Courant; this test only checks the formula.
    const m = createMurState(p);
    expect(m.rx).toBe(0);
    expect(m.ry).toBe(0);
    expect(m.rz).toBe(0);
  });
});
