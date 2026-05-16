// Bipolar (blue → white → red) colormap. Maps v ∈ [-1, 1] to (r, g, b).
export function bipolar(v: number): [number, number, number] {
  const t = Math.max(-1, Math.min(1, v));
  if (t >= 0) {
    // 0 → white, +1 → red
    const r = 255;
    const g = Math.round(255 * (1 - t));
    const b = Math.round(255 * (1 - t));
    return [r, g, b];
  }
  const r = Math.round(255 * (1 + t));
  const g = Math.round(255 * (1 + t));
  const b = 255;
  return [r, g, b];
}
