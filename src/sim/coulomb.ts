// Analytic electric field of a 3D Gaussian charge distribution
//   ρ(r) = (q / (π^(3/2) σ³)) · exp(-r²/σ²)   (total charge = q)
//
// By Gauss's law (in normalized units ε₀ = 1):
//   E_r(r) = q · [erf(r/σ) - (2 r / (σ √π)) · exp(-r²/σ²)] / (4 π r²)
//
// For r → 0 the field is linear in r:
//   E_r(r) ≈ q · r / (3 π^(3/2) σ³)
//
// The function `radialField` returns the magnitude of E along r̂; multiply
// by (Δx_i / r) per component to get a vector value.

const SQRT_PI = Math.sqrt(Math.PI);
const FOUR_PI = 4 * Math.PI;

// Abramowitz & Stegun 7.1.26 — max error ~1.5e-7.
export function erf(x: number): number {
  const sign = x < 0 ? -1 : 1;
  const ax = Math.abs(x);
  const t = 1 / (1 + 0.3275911 * ax);
  const y =
    1 -
    (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t +
      0.254829592) *
      t *
      Math.exp(-ax * ax);
  return sign * y;
}

// For small u = r/σ the closed form erf(u) − (2u/√π)·exp(−u²) suffers
// catastrophic cancellation. The series expansion is
//   (4/√π)·u³·(1/3 − u²/5 + u⁴/14 − u⁶/54 + …)
// Use the series for u < 0.1.
export function radialField(r: number, sigma: number, q: number): number {
  const u = r / sigma;
  if (u < 0.1) {
    const u2 = u * u;
    const series = 1 / 3 - u2 / 5 + (u2 * u2) / 14 - (u2 * u2 * u2) / 54;
    // qDiff = q · (4/√π) · u³ · series
    const qDiff = (q * 4 * u * u * u * series) / SQRT_PI;
    if (r === 0) return 0;
    return qDiff / (FOUR_PI * r * r);
  }
  const expTerm = ((2 * r) / (sigma * SQRT_PI)) * Math.exp(-u * u);
  return (q * (erf(u) - expTerm)) / (FOUR_PI * r * r);
}

// 3D Gaussian volume charge density:
//   ρ(r) = q · exp(-r²/σ²) / (π^(3/2) σ³)
// Integral over all space equals q. Used by the motion-current deposit
// (J = ρ·v) which preserves continuity ∇·J = -∂ρ/∂t exactly for a
// translating profile.
const PI_3_2 = Math.pow(Math.PI, 1.5);
export function gaussianDensity(r: number, sigma: number, q: number): number {
  const u = r / sigma;
  return (q * Math.exp(-u * u)) / (PI_3_2 * sigma * sigma * sigma);
}

// Volume-integrated charge density inside radius r (= total charge for r→∞):
export function enclosedCharge(r: number, sigma: number, q: number): number {
  const u = r / sigma;
  if (u < 0.1) {
    const u2 = u * u;
    const series = 1 / 3 - u2 / 5 + (u2 * u2) / 14 - (u2 * u2 * u2) / 54;
    return (q * 4 * u * u * u * series) / SQRT_PI;
  }
  const expTerm = ((2 * r) / (sigma * SQRT_PI)) * Math.exp(-u * u);
  return q * (erf(u) - expTerm);
}
