import type { SimParams } from "./params";

export interface Fields {
  Ex: Float32Array;
  Ey: Float32Array;
  Ez: Float32Array;
  Bx: Float32Array;
  By: Float32Array;
  Bz: Float32Array;
}

export function idx(i: number, j: number, k: number, Nx: number, Ny: number): number {
  return i + Nx * (j + Ny * k);
}

export function createFields(p: SimParams): Fields {
  const n = p.Nx * p.Ny * p.Nz;
  return {
    Ex: new Float32Array(n),
    Ey: new Float32Array(n),
    Ez: new Float32Array(n),
    Bx: new Float32Array(n),
    By: new Float32Array(n),
    Bz: new Float32Array(n),
  };
}

export function zeroFields(f: Fields): void {
  f.Ex.fill(0);
  f.Ey.fill(0);
  f.Ez.fill(0);
  f.Bx.fill(0);
  f.By.fill(0);
  f.Bz.fill(0);
}
