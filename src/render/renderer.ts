import { idx } from "../sim/fields";
import type { Simulator } from "../sim/simulator";
import { bipolar } from "./colormap";

export interface RendererOptions {
  pixelsPerCell: number;
  vectorStride: number; // draw an E arrow every N cells
  // Absolute heatmap/vector normalization. |Bz|/bzScale (resp. |E|/eScale)
  // is passed through tanh before mapping to color/length. Choose so that a
  // typical click-charge's peak radiation is around 1.0 here.
  bzScale: number;
  eScale: number;
}

export class Renderer {
  readonly canvas: HTMLCanvasElement;
  readonly opts: RendererOptions;
  private ctx: CanvasRenderingContext2D;
  private imageData: ImageData;
  private slice: HTMLCanvasElement;
  private sliceCtx: CanvasRenderingContext2D;

  constructor(canvas: HTMLCanvasElement, sim: Simulator, opts: Partial<RendererOptions> = {}) {
    this.canvas = canvas;
    this.opts = {
      pixelsPerCell: opts.pixelsPerCell ?? 5,
      vectorStride: opts.vectorStride ?? 8,
      bzScale: opts.bzScale ?? 0.01,
      eScale: opts.eScale ?? 0.02,
    };
    const { Nx, Ny } = sim.params;
    canvas.width = Nx * this.opts.pixelsPerCell;
    canvas.height = Ny * this.opts.pixelsPerCell;
    const ctx = canvas.getContext("2d");
    if (!ctx) throw new Error("Canvas 2D context not available");
    this.ctx = ctx;
    this.ctx.imageSmoothingEnabled = false;

    // Render heatmap at native cell resolution first, then upscale.
    this.slice = document.createElement("canvas");
    this.slice.width = Nx;
    this.slice.height = Ny;
    const sctx = this.slice.getContext("2d");
    if (!sctx) throw new Error("Slice canvas 2D context not available");
    this.sliceCtx = sctx;
    this.imageData = this.sliceCtx.createImageData(Nx, Ny);
  }

  draw(sim: Simulator): void {
    this.drawHeatmap(sim);
    this.ctx.imageSmoothingEnabled = false;
    this.ctx.drawImage(this.slice, 0, 0, this.canvas.width, this.canvas.height);
    this.drawVectors(sim);
    this.drawCharges(sim);
  }

  private drawHeatmap(sim: Simulator): void {
    const { Nx, Ny } = sim.params;
    const kMid = Math.floor(sim.zMid);
    const Bz = sim.fields.Bz;
    const data = this.imageData.data;

    const inv = 1 / this.opts.bzScale;
    for (let j = 0; j < Ny; j++) {
      for (let i = 0; i < Nx; i++) {
        const v = Bz[idx(i, j, kMid, Nx, Ny)] * inv;
        const [r, g, b] = bipolar(Math.tanh(v));
        const p = 4 * (i + Nx * (Ny - 1 - j)); // flip y so +y is up
        data[p] = r;
        data[p + 1] = g;
        data[p + 2] = b;
        data[p + 3] = 255;
      }
    }
    this.sliceCtx.putImageData(this.imageData, 0, 0);
  }

  private drawVectors(sim: Simulator): void {
    const { Nx, Ny } = sim.params;
    const kMid = Math.floor(sim.zMid);
    const { Ex, Ey } = sim.fields;
    const px = this.opts.pixelsPerCell;
    const stride = this.opts.vectorStride;

    const inv = 1 / this.opts.eScale;
    const arrowLen = stride * px * 0.675;

    // Collect arrow geometry first so we can draw outline pass then fill pass.
    type Arrow = { cx: number; cy: number; dx: number; dy: number; head: number; ang: number };
    const arrows: Arrow[] = [];
    for (let j = stride / 2; j < Ny; j += stride) {
      for (let i = stride / 2; i < Nx; i += stride) {
        const cx = (i + 0.5) * px;
        const cy = (Ny - j - 0.5) * px;
        const ii = Math.floor(i);
        const jj = Math.floor(j);
        const exV = Ex[idx(ii, jj, kMid, Nx, Ny)];
        const eyV = Ey[idx(ii, jj, kMid, Nx, Ny)];
        const m = Math.hypot(exV, eyV) * inv;
        if (m < 0.05) continue;
        const sat = Math.tanh(m);
        const mag = Math.hypot(exV, eyV);
        const dx = (exV / mag) * arrowLen * sat;
        const dy = -(eyV / mag) * arrowLen * sat;
        const ang = Math.atan2(dy, dx);
        const head = arrowLen * 0.28 * sat;
        arrows.push({ cx, cy, dx, dy, head, ang });
      }
    }

    const drawArrows = () => {
      this.ctx.beginPath();
      for (const { cx, cy, dx, dy, head, ang } of arrows) {
        // Each segment starts with moveTo so they are never joined — prevents
        // the round lineJoin at the tip from making the arrow look curved.
        this.ctx.moveTo(cx, cy);
        this.ctx.lineTo(cx + dx, cy + dy);
        this.ctx.moveTo(cx + dx, cy + dy);
        this.ctx.lineTo(
          cx + dx - head * Math.cos(ang - 0.4),
          cy + dy - head * Math.sin(ang - 0.4),
        );
        this.ctx.moveTo(cx + dx, cy + dy);
        this.ctx.lineTo(
          cx + dx - head * Math.cos(ang + 0.4),
          cy + dy - head * Math.sin(ang + 0.4),
        );
      }
      this.ctx.stroke();
    };

    // Outline pass (white/translucent) for contrast on both dark and bright bg.
    this.ctx.strokeStyle = "rgba(255, 255, 255, 0.55)";
    this.ctx.lineWidth = 4.5;
    this.ctx.lineCap = "round";
    drawArrows();

    // Fill pass (dark) on top of outline.
    this.ctx.strokeStyle = "rgba(15, 15, 15, 0.92)";
    this.ctx.lineWidth = 2.4;
    drawArrows();
  }

  private drawCharges(sim: Simulator): void {
    const { Ny } = sim.params;
    const px = this.opts.pixelsPerCell;
    for (const c of sim.charges) {
      const radius = Math.max(4, c.sigma * px * 0.8);
      // B holds +q (red when q>0); A holds −q (blue when q>0).
      this.drawMarker(c.xB * px, (Ny - c.yB) * px, radius, c.q >= 0);
      this.drawMarker(c.xA * px, (Ny - c.yA) * px, radius, c.q < 0);
    }
  }

  private drawMarker(x: number, y: number, radius: number, positive: boolean): void {
    this.ctx.beginPath();
    this.ctx.arc(x, y, radius, 0, 2 * Math.PI);
    this.ctx.fillStyle = positive ? "rgba(220, 40, 40, 0.85)" : "rgba(40, 80, 220, 0.85)";
    this.ctx.fill();
    this.ctx.strokeStyle = "rgba(0, 0, 0, 0.6)";
    this.ctx.lineWidth = 1;
    this.ctx.stroke();
    this.ctx.fillStyle = "white";
    this.ctx.font = `${Math.round(radius * 1.3)}px sans-serif`;
    this.ctx.textAlign = "center";
    this.ctx.textBaseline = "middle";
    this.ctx.fillText(positive ? "+" : "−", x, y);
  }
}
