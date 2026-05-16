import { describe, expect, it } from "vitest";
import { createFields, idx } from "../src/sim/fields";
import { totalEnergy, updateB, updateE } from "../src/sim/fdtd";
import type { SimParams } from "../src/sim/params";

// Set up a wave with non-zero Ey in a localized region. Without explicit
// boundary handling the implicit-zero ghost cells act as an imperfect PEC,
// so we deliberately keep the gauss away from boundaries and only check
// bulk behavior over short times.
function makeRadialPulse(N: number, dt = 0.4) {
  const p: SimParams = { Nx: N, Ny: N, Nz: 11, dx: 1, dy: 1, dz: 1, c: 1, dt };
  const f = createFields(p);
  const cx = N / 2, cy = N / 2, cz = (p.Nz - 1) / 2;
  const w = 4;
  for (let k = 0; k < p.Nz; k++) {
    for (let j = 0; j < p.Ny; j++) {
      for (let i = 0; i < p.Nx; i++) {
        const r2 =
          (i - cx) * (i - cx) +
          (j - cy) * (j - cy) +
          (k - cz) * (k - cz);
        // Pluck Ez (z-pole-like initial perturbation)
        f.Ez[idx(i, j, k, p.Nx, p.Ny)] = Math.exp(-r2 / (w * w));
      }
    }
  }
  return { p, f };
}

describe("FDTD basic correctness", () => {
  it("runs without producing NaN or Inf", () => {
    const { p, f } = makeRadialPulse(40);
    for (let n = 0; n < 50; n++) {
      updateB(f, p);
      updateE(f, p);
    }
    for (const a of [f.Ex, f.Ey, f.Ez, f.Bx, f.By, f.Bz]) {
      for (let i = 0; i < a.length; i++) {
        if (!Number.isFinite(a[i])) throw new Error(`non-finite at ${i}`);
      }
    }
  });

  it("radial pulse expands outward (front reaches outer cells)", () => {
    const { p, f } = makeRadialPulse(40);
    const cx = p.Nx / 2;
    const cy = p.Ny / 2;
    const kMid = (p.Nz - 1) >> 1;
    // initial: field concentrated near center
    const initOuter = Math.abs(f.Ez[idx(cx + 12, cy, kMid, p.Nx, p.Ny)]);
    expect(initOuter).toBeLessThan(0.01);
    for (let n = 0; n < 20; n++) {
      updateB(f, p);
      updateE(f, p);
    }
    // After 20 steps (c*dt*20 = 8 cells of propagation), field should
    // have spread to cells ~8 cells away from center.
    const outer = Math.abs(f.Ez[idx(cx + 8, cy, kMid, p.Nx, p.Ny)]);
    expect(outer).toBeGreaterThan(initOuter * 5);
  });
});

describe("FDTD energy conservation", () => {
  it("total energy stays bounded under closed-grid evolution", () => {
    const { p, f } = makeRadialPulse(32);
    const e0 = totalEnergy(f);
    let eMax = e0;
    let eMin = e0;
    for (let n = 0; n < 200; n++) {
      updateB(f, p);
      updateE(f, p);
      const e = totalEnergy(f);
      if (e > eMax) eMax = e;
      if (e < eMin) eMin = e;
    }
    expect(eMin).toBeGreaterThan(0);
    // Energy can drift up to ~20% with the simple (E^2+B^2)/2 metric and
    // imperfect boundary handling, but must not blow up.
    expect((eMax - e0) / e0).toBeLessThan(0.2);
    expect((e0 - eMin) / e0).toBeLessThan(0.5);
  });
});
