import type { Fields } from "./fields";
import { idx } from "./fields";
import type { SimParams } from "./params";

// Mur 1st-order absorbing boundary condition.
//
// At each of the 6 faces we apply, to each tangential E component:
//   E^{n+1}[b] = E^n[adj] + r * (E^{n+1}[adj] - E^n[b])
//   r = (c·dt − Δ) / (c·dt + Δ)   (Δ is the spacing normal to that face)
//
// Only the *interior* of each face is updated; edges/corners (where two or
// three faces meet) are left to the natural curl evolution to avoid
// double-application instabilities.

export interface MurState {
  rx: number;
  ry: number;
  rz: number;
  xMin: { Ey0: Float32Array; Ey1: Float32Array; Ez0: Float32Array; Ez1: Float32Array };
  xMax: { Ey0: Float32Array; Ey1: Float32Array; Ez0: Float32Array; Ez1: Float32Array };
  yMin: { Ex0: Float32Array; Ex1: Float32Array; Ez0: Float32Array; Ez1: Float32Array };
  yMax: { Ex0: Float32Array; Ex1: Float32Array; Ez0: Float32Array; Ez1: Float32Array };
  zMin: { Ex0: Float32Array; Ex1: Float32Array; Ey0: Float32Array; Ey1: Float32Array };
  zMax: { Ex0: Float32Array; Ex1: Float32Array; Ey0: Float32Array; Ey1: Float32Array };
}

export function createMurState(p: SimParams): MurState {
  const slabYZ = () => new Float32Array(p.Ny * p.Nz);
  const slabXZ = () => new Float32Array(p.Nx * p.Nz);
  const slabXY = () => new Float32Array(p.Nx * p.Ny);
  return {
    rx: (p.c * p.dt - p.dx) / (p.c * p.dt + p.dx),
    ry: (p.c * p.dt - p.dy) / (p.c * p.dt + p.dy),
    rz: (p.c * p.dt - p.dz) / (p.c * p.dt + p.dz),
    xMin: { Ey0: slabYZ(), Ey1: slabYZ(), Ez0: slabYZ(), Ez1: slabYZ() },
    xMax: { Ey0: slabYZ(), Ey1: slabYZ(), Ez0: slabYZ(), Ez1: slabYZ() },
    yMin: { Ex0: slabXZ(), Ex1: slabXZ(), Ez0: slabXZ(), Ez1: slabXZ() },
    yMax: { Ex0: slabXZ(), Ex1: slabXZ(), Ez0: slabXZ(), Ez1: slabXZ() },
    zMin: { Ex0: slabXY(), Ex1: slabXY(), Ey0: slabXY(), Ey1: slabXY() },
    zMax: { Ex0: slabXY(), Ex1: slabXY(), Ey0: slabXY(), Ey1: slabXY() },
  };
}

export function snapshotBoundary(f: Fields, p: SimParams, m: MurState): void {
  const { Nx, Ny, Nz } = p;

  for (let k = 0; k < Nz; k++) {
    for (let j = 0; j < Ny; j++) {
      const s = j + Ny * k;
      m.xMin.Ey0[s] = f.Ey[idx(0, j, k, Nx, Ny)];
      m.xMin.Ey1[s] = f.Ey[idx(1, j, k, Nx, Ny)];
      m.xMin.Ez0[s] = f.Ez[idx(0, j, k, Nx, Ny)];
      m.xMin.Ez1[s] = f.Ez[idx(1, j, k, Nx, Ny)];
      m.xMax.Ey0[s] = f.Ey[idx(Nx - 1, j, k, Nx, Ny)];
      m.xMax.Ey1[s] = f.Ey[idx(Nx - 2, j, k, Nx, Ny)];
      m.xMax.Ez0[s] = f.Ez[idx(Nx - 1, j, k, Nx, Ny)];
      m.xMax.Ez1[s] = f.Ez[idx(Nx - 2, j, k, Nx, Ny)];
    }
  }

  for (let k = 0; k < Nz; k++) {
    for (let i = 0; i < Nx; i++) {
      const s = i + Nx * k;
      m.yMin.Ex0[s] = f.Ex[idx(i, 0, k, Nx, Ny)];
      m.yMin.Ex1[s] = f.Ex[idx(i, 1, k, Nx, Ny)];
      m.yMin.Ez0[s] = f.Ez[idx(i, 0, k, Nx, Ny)];
      m.yMin.Ez1[s] = f.Ez[idx(i, 1, k, Nx, Ny)];
      m.yMax.Ex0[s] = f.Ex[idx(i, Ny - 1, k, Nx, Ny)];
      m.yMax.Ex1[s] = f.Ex[idx(i, Ny - 2, k, Nx, Ny)];
      m.yMax.Ez0[s] = f.Ez[idx(i, Ny - 1, k, Nx, Ny)];
      m.yMax.Ez1[s] = f.Ez[idx(i, Ny - 2, k, Nx, Ny)];
    }
  }

  for (let j = 0; j < Ny; j++) {
    for (let i = 0; i < Nx; i++) {
      const s = i + Nx * j;
      m.zMin.Ex0[s] = f.Ex[idx(i, j, 0, Nx, Ny)];
      m.zMin.Ex1[s] = f.Ex[idx(i, j, 1, Nx, Ny)];
      m.zMin.Ey0[s] = f.Ey[idx(i, j, 0, Nx, Ny)];
      m.zMin.Ey1[s] = f.Ey[idx(i, j, 1, Nx, Ny)];
      m.zMax.Ex0[s] = f.Ex[idx(i, j, Nz - 1, Nx, Ny)];
      m.zMax.Ex1[s] = f.Ex[idx(i, j, Nz - 2, Nx, Ny)];
      m.zMax.Ey0[s] = f.Ey[idx(i, j, Nz - 1, Nx, Ny)];
      m.zMax.Ey1[s] = f.Ey[idx(i, j, Nz - 2, Nx, Ny)];
    }
  }
}

// Apply Mur on each face's interior. The well-known Mur 1st-order corner
// instability is suppressed by a thin damping layer near each face.
const DAMP_DEPTH = 4;
const DAMP_FLOOR = 0.92;

export function applyMurAbc(f: Fields, p: SimParams, m: MurState): void {
  const { Nx, Ny, Nz } = p;

  for (let k = 1; k < Nz - 1; k++) {
    for (let j = 1; j < Ny - 1; j++) {
      const s = j + Ny * k;
      const b0 = idx(0, j, k, Nx, Ny);
      const b1 = idx(1, j, k, Nx, Ny);
      f.Ey[b0] = m.xMin.Ey1[s] + m.rx * (f.Ey[b1] - m.xMin.Ey0[s]);
      f.Ez[b0] = m.xMin.Ez1[s] + m.rx * (f.Ez[b1] - m.xMin.Ez0[s]);
      const c0 = idx(Nx - 1, j, k, Nx, Ny);
      const c1 = idx(Nx - 2, j, k, Nx, Ny);
      f.Ey[c0] = m.xMax.Ey1[s] + m.rx * (f.Ey[c1] - m.xMax.Ey0[s]);
      f.Ez[c0] = m.xMax.Ez1[s] + m.rx * (f.Ez[c1] - m.xMax.Ez0[s]);
    }
  }

  for (let k = 1; k < Nz - 1; k++) {
    for (let i = 1; i < Nx - 1; i++) {
      const s = i + Nx * k;
      const b0 = idx(i, 0, k, Nx, Ny);
      const b1 = idx(i, 1, k, Nx, Ny);
      f.Ex[b0] = m.yMin.Ex1[s] + m.ry * (f.Ex[b1] - m.yMin.Ex0[s]);
      f.Ez[b0] = m.yMin.Ez1[s] + m.ry * (f.Ez[b1] - m.yMin.Ez0[s]);
      const c0 = idx(i, Ny - 1, k, Nx, Ny);
      const c1 = idx(i, Ny - 2, k, Nx, Ny);
      f.Ex[c0] = m.yMax.Ex1[s] + m.ry * (f.Ex[c1] - m.yMax.Ex0[s]);
      f.Ez[c0] = m.yMax.Ez1[s] + m.ry * (f.Ez[c1] - m.yMax.Ez0[s]);
    }
  }

  for (let j = 1; j < Ny - 1; j++) {
    for (let i = 1; i < Nx - 1; i++) {
      const s = i + Nx * j;
      const b0 = idx(i, j, 0, Nx, Ny);
      const b1 = idx(i, j, 1, Nx, Ny);
      f.Ex[b0] = m.zMin.Ex1[s] + m.rz * (f.Ex[b1] - m.zMin.Ex0[s]);
      f.Ey[b0] = m.zMin.Ey1[s] + m.rz * (f.Ey[b1] - m.zMin.Ey0[s]);
      const c0 = idx(i, j, Nz - 1, Nx, Ny);
      const c1 = idx(i, j, Nz - 2, Nx, Ny);
      f.Ex[c0] = m.zMax.Ex1[s] + m.rz * (f.Ex[c1] - m.zMax.Ex0[s]);
      f.Ey[c0] = m.zMax.Ey1[s] + m.rz * (f.Ey[c1] - m.zMax.Ey0[s]);
    }
  }

  applyBoundaryDamping(f, p);
}

function dampFactor(d: number): number {
  if (d >= DAMP_DEPTH) return 1;
  const t = d / DAMP_DEPTH; // 0 at wall, 1 at inner edge of layer
  return DAMP_FLOOR + (1 - DAMP_FLOOR) * t * t;
}

function applyBoundaryDamping(f: Fields, p: SimParams): void {
  const { Nx, Ny, Nz } = p;
  for (let k = 0; k < Nz; k++) {
    const dz = Math.min(k, Nz - 1 - k);
    for (let j = 0; j < Ny; j++) {
      const dy = Math.min(j, Ny - 1 - j);
      for (let i = 0; i < Nx; i++) {
        const dx = Math.min(i, Nx - 1 - i);
        const d = Math.min(dx, dy, dz);
        if (d >= DAMP_DEPTH) continue;
        const a = dampFactor(d);
        const id = idx(i, j, k, Nx, Ny);
        f.Ex[id] *= a;
        f.Ey[id] *= a;
        f.Ez[id] *= a;
        f.Bx[id] *= a;
        f.By[id] *= a;
        f.Bz[id] *= a;
      }
    }
  }
}
