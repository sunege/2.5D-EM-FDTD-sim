# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

3D Maxwell FDTD simulator with a "thin-3D / 2.5D" configuration (compute is full 3D on a Yee grid with `Nz ≈ 21`; visualization renders the `z = Nz/2` slice). Built as a high-school physics teaching aid — the goals (in `要件.md`) are to make **finite-speed field propagation**, **static fields emerging from radiation**, and **radiation from accelerating charges** visible and intuitive. Educational clarity is the design driver; performance, full-3D rendering, and physical units are explicitly out of scope (see `要件.md` §15).

Stack: Vite + TypeScript + Canvas2D + Vitest, no other runtime deps. Normalized units (`dx = dy = dz = 1`, `c = 1`, `μ₀ = ε₀ = 1`).

## Commands

```bash
npm run dev          # Vite dev server (localhost:5173)
npm run build        # production build to dist/
npm test             # vitest run (single pass)
npm run test:watch   # vitest watch mode
npx tsc --noEmit     # type-check only
npx vitest run tests/<file>.test.ts   # run a single test file
```

## Architecture

### Step pipeline — order matters

`Simulator.step()` is a tight sequence:

```
snapshotBoundary  →  updateB  →  updateE  →  depositCurrent  →  applyMurAbc
```

- `snapshotBoundary` saves `E^n` at the boundary slabs (boundary + adjacent inner slab, for each tangential component on each of 6 faces). Mur uses these one step later, so it must run *before* anything mutates E.
- `depositCurrent` adds `−dt·J` to E *between* the curl update and the ABC. Inserting it anywhere else changes Gauss-law consistency.
- `applyMurAbc` is the last writer of E_tangential at the 6 boundary slabs; intermediate values produced by `updateE` for those cells are intentionally discarded.

### Yee grid layout

All 6 components are stored as `Nx*Ny*Nz` `Float32Array` (SoA). Indexing is `idx = i + Nx*(j + Ny*k)`. Continuous positions are encoded by *which component* is being read:

```
Ex at (i+½, j,   k)        Bx at (i,   j+½, k+½)
Ey at (i,   j+½, k)        By at (i+½, j,   k+½)
Ez at (i,   j,   k+½)      Bz at (i+½, j+½, k)
```

Every curl in `fdtd.ts`, every J deposition offset in `charges.ts`, and every Mur reference in `boundary.ts` follows this exact staggering — when adding code that touches E/B by index, verify the half-cell offset for the component before adjusting any other code.

### Boundary handling (load-bearing — read before changing)

This part is where the algorithm is fragile and well-known to be unstable; the current shape is the result of trial and error documented in this section.

**`updateE` and `updateB` are interior-only.** They deliberately do *not* update:

- E_tangential cells on the 6 outer faces (those slabs are owned by Mur).
- B cells where the curl stencil would step outside the array (`j = Ny−1` for Bx, etc.). Those B cells stay frozen each step; they sit inside the damping layer described below and harmlessly fade.

This convention removes any need for ghost-cell-zero hacks in the curl loops; previous versions that used implicit `ghost = 0` were unstable because PEC ghosts at the lower face were inconsistent with Mur at the upper face.

**Mur 1st-order ABC** is applied only on the *interior* of each face (`j: 1..Ny−2, k: 1..Nz−2` for x-faces, etc.). Edge and corner cells are skipped to avoid the classic 1st-order Mur corner instability that would otherwise cause exponential blow-up over a few hundred steps.

**Boundary damping layer** (`applyBoundaryDamping` in `boundary.ts`, `DAMP_DEPTH = 4`, `DAMP_FLOOR = 0.92`) runs after Mur each step and multiplies all 6 fields by a quadratic ramp inside the 4-cell layer adjacent to each face. This suppresses both the residual oblique-incidence reflection of Mur and the slow corner instability. Without it the demo accumulates ~50× energy growth over 2000 steps; with it the energy floor is ~0.05% of initial pulse, stable indefinitely.

### Charge injection model

`depositCurrent` (`src/sim/charges.ts`) is the physical core of "click to place a charge":

```
J(x, t) = −w(t)/dt · E_Coulomb_radial(|x − x_c|) · r̂
```

where `E_Coulomb_radial` is the analytic field of a 3D Gaussian charge distribution (`src/sim/coulomb.ts`, uses an `erf` approximation and a Taylor series for `u = r/σ < 0.1` to avoid catastrophic cancellation). `w(t)` is the **quintic smoothstep weight per step** (`S(t) = t³(6t² − 15t + 10)` evaluated at step boundaries). The smoothstep, not a constant `1/N_inj`, is essential: constant weighting radiates a boxcar in time and fills the grid with long-lived high-frequency ripples (radiated spectrum decays as `1/ω²`; smoothstep gives `1/ω⁸`).

Continuity `∇·J = −∂ρ/∂t` holds exactly because `∇·E_Coulomb = ρ`. After `injectTotal` steps the cumulative E increment equals the analytic Coulomb field, so the static field "builds up via radiation" — which is the educational point.

Charges currently have no velocity; the planned drag/Lorentz extension (`要件.md` §10.3, §14.1) will add `vx, vy` and a `J = ρv` term using the same `radialField` profile.

### Rendering

`Renderer.draw()` paints in order: Bz heatmap (via `ImageData` on a `Nx × Ny` offscreen canvas, then pixelated upscale) → E vectors (decimated by `vectorStride`) → charge markers. The bipolar colormap (blue/white/red) lives in `colormap.ts`.

**Color scale is absolute, not adaptive.** `bzScale` and `eScale` in `RendererOptions` (default `0.02`) are passed through `tanh` so that `|Bz| = bzScale` saturates to fully red/blue. The previous EMA-of-peak normalization caused a new strong charge to "white-out" the residual ripples from an earlier charge — which looked wrong because the prior fields were still numerically present. Keep this absolute; if a charge with much larger `q` is added later, tune `bzScale` rather than reintroducing adaptive scaling.

Y is flipped at render time (`(Ny − j)`) so that `+y` is up on screen; this affects both the heatmap and the click→grid mapping in `ui/input.ts`.

## Things that look wrong but are intentional

- Outer-edge B cells (`Bx` at `j = Ny−1`, etc.) are never updated by `updateB`. This is correct — they fall inside the boundary damping layer and decay each step.
- The boundary damping layer scales *every* field component, including `Bz` which is what the heatmap shows. This means a wave entering the damping layer visibly fades before "reaching" the edge; that's the intended absorber visualization.
- `radialField(r, σ, q)` switches between closed-form and Taylor series at `r/σ < 0.1`. Don't lower the threshold — `erf(u) − (2u/√π)·exp(−u²)` cancels to single-digit precision for small `u` in Float64.

## Tests

Tests are pure Node (`environment: 'node'` in `vite.config.ts`) — no DOM, no Canvas. They cover the simulation core only; the renderer and UI layer are verified by `tsc --noEmit` plus manual browser testing via `npm run dev`.

Stability tests (`fdtd.test.ts`, `boundary.test.ts`, `simulator.test.ts`) deliberately run for hundreds to thousands of steps and check that energy stays bounded. If you change boundary code or update equations, run these — silent instabilities take ~200 steps to surface.
