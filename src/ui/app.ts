import { Renderer } from "../render/renderer";
import type { Charge } from "../sim/charges";
import { defaultParams } from "../sim/params";
import { Simulator } from "../sim/simulator";
import { eventToGrid } from "./input";

const PIXELS_PER_CELL = 5;
const CHARGE_SIGMA = 4;
const CHARGE_INJECT_STEPS = 10;
const STEPS_PER_FRAME = 1;
const BZ_SCALE = 0.00002;
const E_SCALE = 0.0004;

function requireEl<T extends HTMLElement>(id: string, ctor: new () => T): T {
  const el = document.getElementById(id);
  if (!(el instanceof ctor)) throw new Error(`#${id} not found or wrong type`);
  return el;
}

// Hit-test: returns the topmost charge at (gx, gy), or null. Radius matches
// the visual marker drawn by the renderer (max(2, 0.8·σ) grid cells).
function hitChargeAt(charges: Charge[], gx: number, gy: number): Charge | null {
  for (let k = charges.length - 1; k >= 0; k--) {
    const c = charges[k];
    const radius = Math.max(2, c.sigma * 0.8);
    if (Math.hypot(c.x - gx, c.y - gy) <= radius) return c;
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
  let dragging: Charge | null = null;
  let dragPointerId: number | null = null;

  const formatQ = (q: number) => (q >= 0 ? `+${q.toFixed(1)}` : q.toFixed(1));
  const refreshQDisplay = () => {
    qValue.textContent = formatQ(Number(qSlider.value));
  };
  refreshQDisplay();
  qSlider.addEventListener("input", refreshQDisplay);

  const endDrag = () => {
    if (!dragging) return;
    // Snap target back to current position so the motion deposit stops on
    // the next step — this produces a small deceleration radiation pulse,
    // which is the physically correct behavior of "letting go" of a charge.
    dragging.targetX = dragging.x;
    dragging.targetY = dragging.y;
    dragging = null;
    dragPointerId = null;
    canvas.classList.remove("dragging");
  };

  pauseBtn.addEventListener("click", () => {
    paused = !paused;
    pauseBtn.textContent = paused ? "再開" : "一時停止";
    pauseBtn.classList.toggle("active", paused);
    // While paused, no FDTD steps run; let go of any drag so the charge
    // doesn't drift the moment we unpause.
    if (paused) endDrag();
  });

  resetBtn.addEventListener("click", () => {
    endDrag();
    sim.reset();
  });

  canvas.addEventListener("pointerdown", (evt) => {
    const { gx, gy } = eventToGrid(canvas, evt, sim.params.Nx, sim.params.Ny, PIXELS_PER_CELL);

    // Drag an existing charge if the press lands on its marker (and we're
    // not paused — no deposit happens while paused, so dragging is futile).
    if (!paused) {
      const hit = hitChargeAt(sim.charges, gx, gy);
      if (hit) {
        dragging = hit;
        dragPointerId = evt.pointerId;
        canvas.setPointerCapture(evt.pointerId);
        canvas.classList.add("dragging");
        return;
      }
    }

    // Otherwise place a new charge using the slider value.
    const base = Number(qSlider.value);
    const flip = evt.shiftKey || evt.button === 2;
    const q = flip ? -base : base;
    if (q === 0) return;
    sim.addCharge(gx, gy, q, CHARGE_SIGMA, CHARGE_INJECT_STEPS);
  });

  canvas.addEventListener("pointermove", (evt) => {
    const { gx, gy } = eventToGrid(canvas, evt, sim.params.Nx, sim.params.Ny, PIXELS_PER_CELL);
    if (dragging) {
      dragging.targetX = gx;
      dragging.targetY = gy;
      return;
    }
    // Hover cursor feedback so users discover that charges are draggable.
    const hovering = hitChargeAt(sim.charges, gx, gy) !== null;
    canvas.classList.toggle("hover-charge", hovering && !paused);
  });

  const releaseHandler = (evt: PointerEvent) => {
    if (dragPointerId !== null && evt.pointerId === dragPointerId) endDrag();
  };
  canvas.addEventListener("pointerup", releaseHandler);
  canvas.addEventListener("pointercancel", releaseHandler);

  // Prevent the context menu so right-click can place negative charges.
  canvas.addEventListener("contextmenu", (e) => e.preventDefault());

  function tick(): void {
    if (!paused) {
      for (let i = 0; i < STEPS_PER_FRAME; i++) sim.step();
    }
    renderer.draw(sim);
    requestAnimationFrame(tick);
  }
  requestAnimationFrame(tick);
}
