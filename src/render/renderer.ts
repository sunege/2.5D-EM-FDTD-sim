import type { Charge } from "../sim/charges";
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
    const arrowLen = stride * px * 0.45;

    this.ctx.strokeStyle = "rgba(20, 20, 20, 0.85)";
    this.ctx.lineWidth = 1.2;
    this.ctx.lineCap = "round";
    this.ctx.beginPath();
    for (let j = stride / 2; j < Ny; j += stride) {
      for (let i = stride / 2; i < Nx; i += stride) {
        // Center the arrow on the cell center; flip y for canvas coords.
        const cx = (i + 0.5) * px;
        const cy = (Ny - j - 0.5) * px;
        const ii = Math.floor(i);
        const jj = Math.floor(j);
        // Average to cell center: Ex at (i+0.5, j) → use Ex[i, j]; Ey at (i, j+0.5) → use Ey[i, j]
        const exV = Ex[idx(ii, jj, kMid, Nx, Ny)];
        const eyV = Ey[idx(ii, jj, kMid, Nx, Ny)];
        const m = Math.hypot(exV, eyV) * inv;
        if (m < 0.05) continue;
        const sat = Math.tanh(m);
        // Flip the y component to match the flipped canvas coordinates.
        const dx = (exV / Math.hypot(exV, eyV)) * arrowLen * sat;
        const dy = -(eyV / Math.hypot(exV, eyV)) * arrowLen * sat;
        this.ctx.moveTo(cx, cy);
        this.ctx.lineTo(cx + dx, cy + dy);
        // simple arrowhead
        const ang = Math.atan2(dy, dx);
        const head = arrowLen * 0.25 * sat;
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
    }
    this.ctx.stroke();
  }

  private drawCharges(sim: Simulator): void {
    const { Ny } = sim.params;
    const px = this.opts.pixelsPerCell;
    for (const c of sim.charges) {
      const x = c.x * px;
      const y = (Ny - c.y) * px;
      const radius = Math.max(4, c.sigma * px * 0.8);
      this.ctx.beginPath();
      this.ctx.arc(x, y, radius, 0, 2 * Math.PI);
      this.ctx.fillStyle = c.q >= 0 ? "rgba(220, 40, 40, 0.85)" : "rgba(40, 80, 220, 0.85)";
      this.ctx.fill();
      this.ctx.strokeStyle = "rgba(0, 0, 0, 0.6)";
      this.ctx.lineWidth = 1;
      this.ctx.stroke();
      // sign label
      this.ctx.fillStyle = "white";
      this.ctx.font = `${Math.round(radius * 1.3)}px sans-serif`;
      this.ctx.textAlign = "center";
      this.ctx.textBaseline = "middle";
      this.ctx.fillText(c.q >= 0 ? "+" : "−", x, y);
    }
  }
}
