import { describe, expect, it } from "vitest";
import { depositCurrent, makeDipole } from "../src/sim/charges";
import { createFields, idx, type Fields } from "../src/sim/fields";
import type { SimParams } from "../src/sim/params";

const p: SimParams = {
  Nx: 48, Ny: 48, Nz: 21,
  dx: 1, dy: 1, dz: 1, c: 1, dt: 0.4,
};

// Discrete divergence ∇·E at cell corner (i, j, k). For our wire injection
// only Ex is modified, so Ey and Ez terms contribute zero — but we still
// include them so the helper is the genuine 3-component divergence.
const divAt = (f: Fields, i: number, j: number, k: number): number =>
  f.Ex[idx(i, j, k, p.Nx, p.Ny)] - f.Ex[idx(i - 1, j, k, p.Nx, p.Ny)]
  + f.Ey[idx(i, j, k, p.Nx, p.Ny)] - f.Ey[idx(i, j - 1, k, p.Nx, p.Ny)]
  + f.Ez[idx(i, j, k, p.Nx, p.Ny)] - f.Ez[idx(i, j, k - 1, p.Nx, p.Ny)];

describe("dipole wire current injection", () => {
  it("accumulates +q at B and −q at A after full injection", () => {
    // cx=24, sigma=2 → d = 2σ = 4 → A at x=20, B at x=28.
    // The discrete divergence shells live at i=20 (A) and i=28 (B).
    const f = createFields(p);
    const zMid = (p.Nz - 1) / 2;
    const c = makeDipole(24, 24, 1, 2, 10);
    for (let n = 0; n < 10; n++) depositCurrent(f, p, [c], zMid);

    let qB = 0, qA = 0;
    for (let k = 1; k < p.Nz; k++) {
      for (let j = 1; j < p.Ny; j++) {
        qB += divAt(f, 28, j, k);
        qA += divAt(f, 20, j, k);
      }
    }
    expect(qB).toBeCloseTo(1, 5);
    expect(qA).toBeCloseTo(-1, 5);
  });

  it("net charge over the entire domain is zero", () => {
    // Sum of ∇·E over all cells telescopes to a boundary flux. The wire
    // injection only modifies Ex deep inside the domain, so the boundary
    // E stays zero and the total charge must be exactly 0.
    const f = createFields(p);
    const zMid = (p.Nz - 1) / 2;
    const c = makeDipole(24, 24, 1, 2, 10);
    for (let n = 0; n < 10; n++) depositCurrent(f, p, [c], zMid);

    let total = 0;
    for (let k = 1; k < p.Nz; k++) {
      for (let j = 1; j < p.Ny; j++) {
        for (let i = 1; i < p.Nx; i++) {
          total += divAt(f, i, j, k);
        }
      }
    }
    expect(total).toBeCloseTo(0, 5);
  });

  it("Ex is unchanged outside the wire region in x", () => {
    // Wire from ax=20 to bx=28: Ex at (i+½, j, k) is modified only for
    // i ∈ [20, 27]. Probe i=19 (Yee 19.5, outside left) and i=28 (Yee 28.5,
    // outside right) — both must still be exactly zero.
    const f = createFields(p);
    const zMid = (p.Nz - 1) / 2;
    const c = makeDipole(24, 24, 1, 2, 10);
    for (let n = 0; n < 10; n++) depositCurrent(f, p, [c], zMid);

    const k = Math.round(zMid);
    expect(f.Ex[idx(19, 24, k, p.Nx, p.Ny)]).toBe(0);
    expect(f.Ex[idx(28, 24, k, p.Nx, p.Ny)]).toBe(0);
  });

  it("first-step weight is small (smooth start, not boxcar)", () => {
    // Smoothstep starts flat, so the first step deposits much less than
    // 1/N_inj of the total. Boxcar would give 1/20 = 5%; quintic smoothstep
    // gives ≈ S(1/20) − S(0) ≈ 1.16e−3 = 0.12%.
    const f1 = createFields(p);
    const zMid = (p.Nz - 1) / 2;
    const c1 = makeDipole(24, 24, 1, 2, 20);
    depositCurrent(f1, p, [c1], zMid);
    const k = Math.round(zMid);
    const after1 = f1.Ex[idx(23, 24, k, p.Nx, p.Ny)];

    const f2 = createFields(p);
    const c2 = makeDipole(24, 24, 1, 2, 20);
    for (let n = 0; n < 20; n++) depositCurrent(f2, p, [c2], zMid);
    const final = f2.Ex[idx(23, 24, k, p.Nx, p.Ny)];

    const ratio = after1 / final;
    expect(ratio).toBeGreaterThan(0);
    expect(ratio).toBeLessThan(0.005);
  });

  it("stops depositing after injectTotal steps", () => {
    const f = createFields(p);
    const zMid = (p.Nz - 1) / 2;
    const c = makeDipole(24, 24, 1, 2, 10);
    for (let n = 0; n < 10; n++) depositCurrent(f, p, [c], zMid);
    const snapshot = new Float32Array(f.Ex);
    for (let n = 0; n < 10; n++) depositCurrent(f, p, [c], zMid);
    for (let i = 0; i < snapshot.length; i++) {
      expect(f.Ex[i]).toBe(snapshot[i]);
    }
  });

  it("flips sign correctly with negative q", () => {
    // c.q < 0 → wire current effectively reverses: B gets −|q|, A gets +|q|.
    const f = createFields(p);
    const zMid = (p.Nz - 1) / 2;
    const c = makeDipole(24, 24, -1, 2, 10);
    for (let n = 0; n < 10; n++) depositCurrent(f, p, [c], zMid);

    let qB = 0, qA = 0;
    for (let k = 1; k < p.Nz; k++) {
      for (let j = 1; j < p.Ny; j++) {
        qB += divAt(f, 28, j, k);
        qA += divAt(f, 20, j, k);
      }
    }
    expect(qB).toBeCloseTo(-1, 5);
    expect(qA).toBeCloseTo(1, 5);
  });
});
