import { applyMurAbc, createMurState, snapshotBoundary, type MurState } from "./boundary";
import { advancePositions, depositCurrent, depositMotion, makeDipole, type DipoleCharge } from "./charges";
import { createFields, type Fields, zeroFields } from "./fields";
import { totalEnergy, updateB, updateE } from "./fdtd";
import { defaultParams, type SimParams, assertStable } from "./params";

export class Simulator {
  readonly params: SimParams;
  readonly fields: Fields;
  readonly charges: DipoleCharge[] = [];
  readonly zMid: number;
  private mur: MurState;
  stepCount = 0;

  constructor(params: SimParams = defaultParams) {
    assertStable(params);
    this.params = params;
    this.fields = createFields(params);
    this.mur = createMurState(params);
    this.zMid = (params.Nz - 1) / 2;
  }

  addCharge(cx: number, cy: number, q: number, sigma = 2, injectTotal = 20): void {
    this.charges.push(makeDipole(cx, cy, q, sigma, injectTotal));
  }

  step(): void {
    snapshotBoundary(this.fields, this.params, this.mur);
    updateB(this.fields, this.params);
    updateE(this.fields, this.params);
    depositCurrent(this.fields, this.params, this.charges, this.zMid);
    depositMotion(this.fields, this.params, this.charges, this.zMid);
    applyMurAbc(this.fields, this.params, this.mur);
    advancePositions(this.charges, this.params.dt);
    this.stepCount++;
  }

  reset(): void {
    zeroFields(this.fields);
    this.charges.length = 0;
    this.stepCount = 0;
  }

  energy(): number {
    return totalEnergy(this.fields);
  }
}
