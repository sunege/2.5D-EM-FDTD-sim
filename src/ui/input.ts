export interface ClickPos {
  gx: number;
  gy: number;
}

// Convert a pointer event over `canvas` into continuous grid coordinates.
// y is flipped so that +y is up on screen.
export function eventToGrid(
  canvas: HTMLCanvasElement,
  evt: PointerEvent | MouseEvent,
  Nx: number,
  Ny: number,
  pixelsPerCell: number,
): ClickPos {
  const rect = canvas.getBoundingClientRect();
  const sx = canvas.width / rect.width;
  const sy = canvas.height / rect.height;
  const px = (evt.clientX - rect.left) * sx;
  const py = (evt.clientY - rect.top) * sy;
  const gx = px / pixelsPerCell;
  const gy = Ny - py / pixelsPerCell;
  return {
    gx: Math.max(0, Math.min(Nx - 1, gx)),
    gy: Math.max(0, Math.min(Ny - 1, gy)),
  };
}
