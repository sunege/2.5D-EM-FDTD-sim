import { Renderer } from "../render/renderer";
import { defaultParams } from "../sim/params";
import type { DipoleCharge } from "../sim/charges";
import { Simulator } from "../sim/simulator";
import { eventToGrid } from "./input";

const PIXELS_PER_CELL = 5;
const CHARGE_SIGMA = 4;
const CHARGE_INJECT_STEPS = 10;
const STEPS_PER_FRAME = 1;
const BZ_SCALE = 0.0004;
const E_SCALE = 0.0002;

// Hit radius in grid cells for dragging a charge endpoint.
const HIT_RADIUS_CELLS = 1.5;

// Maximum drag speed in units of c (c=1 in sim units). Capping at sub-luminal
// speed prevents the PIC ghost-charge artifact: if v >> c, the charge
// teleports multiple σ in one step and J is deposited at the old position
// while the charge density appears at the new position, leaving ∇·E ≠ 0.
const MAX_DRAG_SPEED = 0.5;

function requireEl<T extends HTMLElement>(id: string, ctor: new () => T): T {
  const el = document.getElementById(id);
  if (!(el instanceof ctor)) throw new Error(`#${id} not found or wrong type`);
  return el;
}

function hitChargeEndpoint(
  sim: Simulator,
  gx: number,
  gy: number,
): { charge: DipoleCharge; endpoint: "A" | "B" } | null {
  for (const c of sim.charges) {
    const hr = Math.max(HIT_RADIUS_CELLS, c.sigma * 0.8);
    if (Math.hypot(gx - c.xB, gy - c.yB) <= hr) return { charge: c, endpoint: "B" };
    if (Math.hypot(gx - c.xA, gy - c.yA) <= hr) return { charge: c, endpoint: "A" };
  }
  return null;
}

export function startApp(): void {
  const canvas = requireEl("sim-canvas", HTMLCanvasElement);
  const qSlider = requireEl("q-slider", HTMLInputElement);
  const qValue = requireEl("q-value", HTMLSpanElement);
  const pauseBtn = requireEl("pause-btn", HTMLButtonElement);
  const resetBtn = requireEl("reset-btn", HTMLButtonElement);

  const sim = new Simulator(defaultParams);
  const renderer = new Renderer(canvas, sim, {
    pixelsPerCell: PIXELS_PER_CELL,
    vectorStride: 4,
    bzScale: BZ_SCALE,
    eScale: E_SCALE,
  });

  let paused = false;

  // Drag state: tracks which endpoint is being dragged and the pointer's
  // current grid position. Velocity is recomputed each frame in tick() as
  // (target − physics_position) / dt so the charge reaches the pointer in
  // exactly STEPS_PER_FRAME FDTD steps without position jumps.
  let dragging: {
    charge: DipoleCharge;
    endpoint: "A" | "B";
    pointerId: number;
    targetGx: number;
    targetGy: number;
  } | null = null;

  const formatQ = (q: number) => (q >= 0 ? `+${q.toFixed(1)}` : q.toFixed(1));
  const refreshQDisplay = () => {
    qValue.textContent = formatQ(Number(qSlider.value));
  };
  refreshQDisplay();
  qSlider.addEventListener("input", refreshQDisplay);

  pauseBtn.addEventListener("click", () => {
    paused = !paused;
    pauseBtn.textContent = paused ? "再開" : "一時停止";
    pauseBtn.classList.toggle("active", paused);
  });

  resetBtn.addEventListener("click", () => {
    sim.reset();
    dragging = null;
  });

  canvas.addEventListener("pointerdown", (evt) => {
    const { gx, gy } = eventToGrid(canvas, evt, sim.params.Nx, sim.params.Ny, PIXELS_PER_CELL);

    // Drag existing endpoint if close enough; otherwise place a new charge.
    const hit = hitChargeEndpoint(sim, gx, gy);
    if (hit) {
      canvas.setPointerCapture(evt.pointerId);
      dragging = { charge: hit.charge, endpoint: hit.endpoint, pointerId: evt.pointerId, targetGx: gx, targetGy: gy };
      return;
    }

    const base = Number(qSlider.value);
    const flip = evt.shiftKey || evt.button === 2;
    const q = flip ? -base : base;
    if (q === 0) return;
    sim.addCharge(gx, gy, q, CHARGE_SIGMA, CHARGE_INJECT_STEPS);
  });

  canvas.addEventListener("pointermove", (evt) => {
    if (!dragging || evt.pointerId !== dragging.pointerId) return;
    const { gx, gy } = eventToGrid(canvas, evt, sim.params.Nx, sim.params.Ny, PIXELS_PER_CELL);
    // Only update the target; velocity is computed in tick() each frame so
    // the physics position smoothly tracks the pointer via advancePositions.
    dragging.targetGx = gx;
    dragging.targetGy = gy;
  });

  const endDrag = (pointerId: number) => {
    if (!dragging || dragging.pointerId !== pointerId) return;
    // Zero velocity so advancePositions does not carry the charge further.
    const c = dragging.charge;
    if (dragging.endpoint === "B") { c.vxB = 0; c.vyB = 0; }
    else { c.vxA = 0; c.vyA = 0; }
    dragging = null;
  };
  canvas.addEventListener("pointerup", (evt) => endDrag(evt.pointerId));
  canvas.addEventListener("pointercancel", (evt) => endDrag(evt.pointerId));

  // Prevent the context menu so right-click can place negative charges.
  canvas.addEventListener("contextmenu", (e) => e.preventDefault());

  document.addEventListener("keydown", (e) => {
    if (e.code === "Space" && !e.repeat) {
      e.preventDefault();
      pauseBtn.click();
    }
  });

  function tick(): void {
    if (!paused) {
      // Recompute velocity each frame so the physics position reaches the
      // pointer target in exactly STEPS_PER_FRAME steps, without jumping.
      if (dragging) {
        const c = dragging.charge;
        const dtFrame = sim.params.dt * STEPS_PER_FRAME;
        if (dragging.endpoint === "B") {
          let vx = (dragging.targetGx - c.xB) / dtFrame;
          let vy = (dragging.targetGy - c.yB) / dtFrame;
          const spd = Math.hypot(vx, vy);
          if (spd > MAX_DRAG_SPEED) { vx *= MAX_DRAG_SPEED / spd; vy *= MAX_DRAG_SPEED / spd; }
          c.vxB = vx; c.vyB = vy;
        } else {
          let vx = (dragging.targetGx - c.xA) / dtFrame;
          let vy = (dragging.targetGy - c.yA) / dtFrame;
          const spd = Math.hypot(vx, vy);
          if (spd > MAX_DRAG_SPEED) { vx *= MAX_DRAG_SPEED / spd; vy *= MAX_DRAG_SPEED / spd; }
          c.vxA = vx; c.vyA = vy;
        }
      }
      for (let i = 0; i < STEPS_PER_FRAME; i++) sim.step();
    }
    renderer.draw(sim);
    requestAnimationFrame(tick);
  }
  requestAnimationFrame(tick);
}
