import { describe, expect, it } from "vitest";
import { enclosedCharge, erf, gaussianDensity, radialField } from "../src/sim/coulomb";

describe("erf", () => {
  it("erf(0) = 0", () => {
    expect(erf(0)).toBeCloseTo(0, 6);
  });
  it("erf(inf) → 1", () => {
    expect(erf(10)).toBeCloseTo(1, 6);
    expect(erf(-10)).toBeCloseTo(-1, 6);
  });
  it("matches reference values", () => {
    // erf(1) ≈ 0.8427007929
    expect(erf(1)).toBeCloseTo(0.8427007929, 5);
    // erf(0.5) ≈ 0.5204998778
    expect(erf(0.5)).toBeCloseTo(0.5204998778, 5);
    // erf(2) ≈ 0.9953222650
    expect(erf(2)).toBeCloseTo(0.9953222650, 5);
  });
});

describe("radialField for Gaussian charge", () => {
  const sigma = 2;
  const q = 1;

  it("is zero at the origin", () => {
    expect(radialField(0, sigma, q)).toBeCloseTo(0, 8);
  });

  it("is linear in r for small r", () => {
    const r1 = 1e-3 * sigma;
    const r2 = 2e-3 * sigma;
    const e1 = radialField(r1, sigma, q);
    const e2 = radialField(r2, sigma, q);
    expect(e2 / e1).toBeCloseTo(2, 3);
  });

  it("approaches q/(4π r²) far from the charge", () => {
    const r = 50 * sigma;
    const E = radialField(r, sigma, q);
    const expected = q / (4 * Math.PI * r * r);
    expect(E).toBeCloseTo(expected, 6);
  });

  it("matches the small-r linear approximation", () => {
    const r = 0.001;
    const got = radialField(r, sigma, q);
    const expected =
      (q * r) / (3 * Math.sqrt(Math.PI) * Math.PI * sigma * sigma * sigma);
    expect(got).toBeCloseTo(expected, 8);
  });

  it("scales linearly with q", () => {
    const r = sigma;
    expect(radialField(r, sigma, 2)).toBeCloseTo(2 * radialField(r, sigma, 1), 6);
  });
});

describe("gaussianDensity", () => {
  const sigma = 2;
  const q = 1;

  it("peak value matches ρ(0) = q / (π^(3/2) σ³)", () => {
    const expected = q / (Math.pow(Math.PI, 1.5) * sigma * sigma * sigma);
    expect(gaussianDensity(0, sigma, q)).toBeCloseTo(expected, 8);
  });

  it("decays as exp(-r²/σ²)", () => {
    const at0 = gaussianDensity(0, sigma, q);
    const atSigma = gaussianDensity(sigma, sigma, q);
    expect(atSigma / at0).toBeCloseTo(Math.exp(-1), 6);
  });

  it("scales linearly with q", () => {
    expect(gaussianDensity(1.3, sigma, 3)).toBeCloseTo(3 * gaussianDensity(1.3, sigma, 1), 8);
  });

  it("3D volume integral equals q", () => {
    // Midpoint rule on a 6σ-radius cube. Tail outside is ≲ exp(-36) ≈ 0.
    const N = 60;
    const halfW = 6 * sigma;
    const dx = (2 * halfW) / N;
    let sum = 0;
    for (let k = 0; k < N; k++) {
      const zk = -halfW + (k + 0.5) * dx;
      for (let j = 0; j < N; j++) {
        const yj = -halfW + (j + 0.5) * dx;
        for (let i = 0; i < N; i++) {
          const xi = -halfW + (i + 0.5) * dx;
          const r = Math.sqrt(xi * xi + yj * yj + zk * zk);
          sum += gaussianDensity(r, sigma, q);
        }
      }
    }
    sum *= dx * dx * dx;
    expect(sum).toBeCloseTo(q, 2);
  });
});

describe("enclosedCharge", () => {
  const sigma = 2;
  const q = 1;

  it("approaches q for large r", () => {
    expect(enclosedCharge(50 * sigma, sigma, q)).toBeCloseTo(q, 4);
  });

  it("is zero at origin", () => {
    expect(enclosedCharge(0, sigma, q)).toBeCloseTo(0, 8);
  });

  it("∇·E = ρ identity: E_r(r) · 4π r² = Q_enc(r)", () => {
    for (const r of [0.5, 1, 2, 4, 8]) {
      const E = radialField(r, sigma, q);
      const flux = E * 4 * Math.PI * r * r;
      expect(flux).toBeCloseTo(enclosedCharge(r, sigma, q), 6);
    }
  });
});
