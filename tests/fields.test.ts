import { describe, expect, it } from "vitest";
import { createFields, idx, zeroFields } from "../src/sim/fields";
import { assertStable, courantLimit, defaultParams } from "../src/sim/params";

describe("idx", () => {
  it("maps (0,0,0) to 0", () => {
    expect(idx(0, 0, 0, 8, 8)).toBe(0);
  });

  it("increments by 1 along x", () => {
    expect(idx(1, 0, 0, 8, 8) - idx(0, 0, 0, 8, 8)).toBe(1);
  });

  it("increments by Nx along y", () => {
    expect(idx(0, 1, 0, 8, 8) - idx(0, 0, 0, 8, 8)).toBe(8);
  });

  it("increments by Nx*Ny along z", () => {
    expect(idx(0, 0, 1, 8, 8) - idx(0, 0, 0, 8, 8)).toBe(64);
  });

  it("reaches Nx*Ny*Nz - 1 at the top corner", () => {
    expect(idx(7, 7, 4, 8, 8)).toBe(8 * 8 * 5 - 1);
  });
});

describe("createFields", () => {
  it("creates 6 Float32Array buffers of Nx*Ny*Nz", () => {
    const f = createFields(defaultParams);
    const expected = defaultParams.Nx * defaultParams.Ny * defaultParams.Nz;
    expect(f.Ex).toBeInstanceOf(Float32Array);
    expect(f.Ex.length).toBe(expected);
    expect(f.Ey.length).toBe(expected);
    expect(f.Ez.length).toBe(expected);
    expect(f.Bx.length).toBe(expected);
    expect(f.By.length).toBe(expected);
    expect(f.Bz.length).toBe(expected);
  });

  it("zeroFields resets buffers to zero", () => {
    const f = createFields(defaultParams);
    f.Ex[0] = 1.5;
    f.Bz[10] = -2.0;
    zeroFields(f);
    expect(f.Ex[0]).toBe(0);
    expect(f.Bz[10]).toBe(0);
  });
});

describe("courant", () => {
  it("default params satisfy stability", () => {
    expect(defaultParams.c * defaultParams.dt).toBeLessThan(courantLimit(defaultParams));
    expect(() => assertStable(defaultParams)).not.toThrow();
  });

  it("rejects unstable dt", () => {
    expect(() => assertStable({ ...defaultParams, dt: 1.0 })).toThrow(/Courant/);
  });
});
