export interface SimParams {
  Nx: number;
  Ny: number;
  Nz: number;
  dx: number;
  dy: number;
  dz: number;
  c: number;
  dt: number;
}

export const defaultParams: SimParams = {
  Nx: 128,
  Ny: 128,
  Nz: 21,
  dx: 1,
  dy: 1,
  dz: 1,
  c: 1,
  dt: 0.4,
};

export function courantLimit(p: SimParams): number {
  return 1 / Math.sqrt(1 / (p.dx * p.dx) + 1 / (p.dy * p.dy) + 1 / (p.dz * p.dz));
}

export function assertStable(p: SimParams): void {
  const limit = courantLimit(p);
  if (p.c * p.dt > limit) {
    throw new Error(
      `Courant condition violated: c*dt=${p.c * p.dt} > ${limit.toFixed(4)}`,
    );
  }
}
