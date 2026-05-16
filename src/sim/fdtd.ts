import type { Fields } from "./fields";
import { idx } from "./fields";
import type { SimParams } from "./params";

// Interior-only B updates: skip the outermost B cells in the direction
// where the curl stencil would reach outside the array. Those skipped B
// cells stay at their previous value, which is harmless because they sit
// in the absorbing region overwritten by the ABC of the adjacent face on
// the next E step.
export function updateB(f: Fields, p: SimParams): void {
  const { Nx, Ny, Nz, dx, dy, dz, dt } = p;
  const { Ex, Ey, Ez, Bx, By, Bz } = f;
  const invDx = 1 / dx;
  const invDy = 1 / dy;
  const invDz = 1 / dz;

  // Bx needs Ez(j+1), Ey(k+1): skip j=Ny-1 and k=Nz-1.
  for (let k = 0; k < Nz - 1; k++) {
    for (let j = 0; j < Ny - 1; j++) {
      for (let i = 0; i < Nx; i++) {
        const c000 = idx(i, j, k, Nx, Ny);
        const cj1 = idx(i, j + 1, k, Nx, Ny);
        const ck1 = idx(i, j, k + 1, Nx, Ny);
        const curlEx = (Ez[cj1] - Ez[c000]) * invDy - (Ey[ck1] - Ey[c000]) * invDz;
        Bx[c000] -= dt * curlEx;
      }
    }
  }

  // By needs Ex(k+1), Ez(i+1): skip i=Nx-1 and k=Nz-1.
  for (let k = 0; k < Nz - 1; k++) {
    for (let j = 0; j < Ny; j++) {
      for (let i = 0; i < Nx - 1; i++) {
        const c000 = idx(i, j, k, Nx, Ny);
        const ci1 = idx(i + 1, j, k, Nx, Ny);
        const ck1 = idx(i, j, k + 1, Nx, Ny);
        const curlEy = (Ex[ck1] - Ex[c000]) * invDz - (Ez[ci1] - Ez[c000]) * invDx;
        By[c000] -= dt * curlEy;
      }
    }
  }

  // Bz needs Ey(i+1), Ex(j+1): skip i=Nx-1 and j=Ny-1.
  for (let k = 0; k < Nz; k++) {
    for (let j = 0; j < Ny - 1; j++) {
      for (let i = 0; i < Nx - 1; i++) {
        const c000 = idx(i, j, k, Nx, Ny);
        const ci1 = idx(i + 1, j, k, Nx, Ny);
        const cj1 = idx(i, j + 1, k, Nx, Ny);
        const curlEz = (Ey[ci1] - Ey[c000]) * invDx - (Ex[cj1] - Ex[c000]) * invDy;
        Bz[c000] -= dt * curlEz;
      }
    }
  }
}

// Update E only on cells where every Yee-stencil neighbor is in-bounds;
// tangential-E values on the 6 boundary slabs are NOT updated here. Those
// will be set by the Mur ABC step (or any equivalent boundary handler).
export function updateE(f: Fields, p: SimParams): void {
  const { Nx, Ny, Nz, dx, dy, dz, dt } = p;
  const { Ex, Ey, Ez, Bx, By, Bz } = f;
  const invDx = 1 / dx;
  const invDy = 1 / dy;
  const invDz = 1 / dz;

  // Ex at (i+0.5, j, k): needs Bz(j-1), By(k-1). Tangential on y/z faces.
  for (let k = 1; k < Nz - 1; k++) {
    for (let j = 1; j < Ny - 1; j++) {
      for (let i = 0; i < Nx; i++) {
        const c000 = idx(i, j, k, Nx, Ny);
        const cjm = idx(i, j - 1, k, Nx, Ny);
        const ckm = idx(i, j, k - 1, Nx, Ny);
        const curlBx = (Bz[c000] - Bz[cjm]) * invDy - (By[c000] - By[ckm]) * invDz;
        Ex[c000] += dt * curlBx;
      }
    }
  }

  // Ey at (i, j+0.5, k): tangential on x/z faces.
  for (let k = 1; k < Nz - 1; k++) {
    for (let j = 0; j < Ny; j++) {
      for (let i = 1; i < Nx - 1; i++) {
        const c000 = idx(i, j, k, Nx, Ny);
        const cim = idx(i - 1, j, k, Nx, Ny);
        const ckm = idx(i, j, k - 1, Nx, Ny);
        const curlBy = (Bx[c000] - Bx[ckm]) * invDz - (Bz[c000] - Bz[cim]) * invDx;
        Ey[c000] += dt * curlBy;
      }
    }
  }

  // Ez at (i, j, k+0.5): tangential on x/y faces.
  for (let k = 0; k < Nz; k++) {
    for (let j = 1; j < Ny - 1; j++) {
      for (let i = 1; i < Nx - 1; i++) {
        const c000 = idx(i, j, k, Nx, Ny);
        const cim = idx(i - 1, j, k, Nx, Ny);
        const cjm = idx(i, j - 1, k, Nx, Ny);
        const curlBz = (By[c000] - By[cim]) * invDx - (Bx[c000] - Bx[cjm]) * invDy;
        Ez[c000] += dt * curlBz;
      }
    }
  }
}

export function totalEnergy(f: Fields): number {
  let sum = 0;
  const { Ex, Ey, Ez, Bx, By, Bz } = f;
  for (let i = 0; i < Ex.length; i++) {
    sum +=
      Ex[i] * Ex[i] +
      Ey[i] * Ey[i] +
      Ez[i] * Ez[i] +
      Bx[i] * Bx[i] +
      By[i] * By[i] +
      Bz[i] * Bz[i];
  }
  return 0.5 * sum;
}
