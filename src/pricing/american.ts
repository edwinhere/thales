module Thales.Pricing {
    var app = getModule();
    /**
     * 
     * @author Edwin Jose Palathinkal <edwinhere@gmail.com>
     */
    export class American {
        public CND(X: number): number {
            let functionReturnValue: number = 0.0;
            let L: number = 0.0;
            let K: number = 0.0;
            let a1: number = 0.31938153;
            let a2: number = -0.356563782;
            let a3: number = 1.781477937;
            let a4: number = -1.821255978;
            let a5: number = 1.330274429;
            L = Math.abs(X);
            K = 1.0 / (1.0 + 0.2316419 * L);
            functionReturnValue = 1.0 - 1.0 / Math.sqrt(2.0 * Math.PI) * Math.exp(-Math.pow(L, 2.0) / 2.0) * (a1 * K + a2 * Math.pow(K, 2.0) + a3 * Math.pow(K, 3.0) + a4 * Math.pow(K, 4.0) + a5 * Math.pow(K, 5.0));
            if (X < 0.0) {
                functionReturnValue = 1.0 - functionReturnValue;
            }
            return functionReturnValue;
        }

        public CBND(x: number, y: number, rho: number): number {
            let i: number = 0;
            let ISs: number = 0;
            let LG: number = 0;
            let NG: number = 0;
            let XX: number[][] = <any>(function (dims) { let allocate = function (dims) { if (dims.length == 0) { return 0; } else { let array = []; for (let i = 0; i < dims[0]; i++) { array.push(allocate(dims.slice(1))); } return array; } }; return allocate(dims); })([11, 4]);
            let W: number[][] = <any>(function (dims) { let allocate = function (dims) { if (dims.length == 0) { return 0; } else { let array = []; for (let i = 0; i < dims[0]; i++) { array.push(allocate(dims.slice(1))); } return array; } }; return allocate(dims); })([11, 4]);
            let h: number = 0.0;
            let k: number = 0.0;
            let hk: number = 0.0;
            let hs: number = 0.0;
            let BVN: number = 0.0;
            let Ass: number = 0.0;
            let asr: number = 0.0;
            let sn: number = 0.0;
            let A: number = 0.0;
            let b: number = 0.0;
            let bs: number = 0.0;
            let c: number = 0.0;
            let d: number = 0.0;
            let xs: number = 0.0;
            let rs: number = 0.0;
            W[1][1] = 0.17132449237917;
            XX[1][1] = -0.932469514203152;
            W[2][1] = 0.360761573048138;
            XX[2][1] = -0.661209386466265;
            W[3][1] = 0.46791393457269;
            XX[3][1] = -0.238619186083197;
            W[1][2] = 0.0471753363865118;
            XX[1][2] = -0.981560634246719;
            W[2][2] = 0.106939325995318;
            XX[2][2] = -0.904117256370475;
            W[3][2] = 0.160078328543346;
            XX[3][2] = -0.769902674194305;
            W[4][2] = 0.203167426723066;
            XX[4][2] = -0.587317954286617;
            W[5][2] = 0.233492536538355;
            XX[5][2] = -0.36783149899818;
            W[6][2] = 0.249147045813403;
            XX[6][2] = -0.125233408511469;
            W[1][3] = 0.0176140071391521;
            XX[1][3] = -0.993128599185095;
            W[2][3] = 0.0406014298003869;
            XX[2][3] = -0.963971927277914;
            W[3][3] = 0.0626720483341091;
            XX[3][3] = -0.912234428251326;
            W[4][3] = 0.0832767415767048;
            XX[4][3] = -0.839116971822219;
            W[5][3] = 0.10193011981724;
            XX[5][3] = -0.746331906460151;
            W[6][3] = 0.118194531961518;
            XX[6][3] = -0.636053680726515;
            W[7][3] = 0.131688638449177;
            XX[7][3] = -0.510867001950827;
            W[8][3] = 0.142096109318382;
            XX[8][3] = -0.37370608871542;
            W[9][3] = 0.149172986472604;
            XX[9][3] = -0.227785851141645;
            W[10][3] = 0.152753387130726;
            XX[10][3] = -0.0765265211334973;
            if (Math.abs(rho) < 0.3) {
                NG = 1;
                LG = 3;
            } else if (Math.abs(rho) < 0.75) {
                NG = 2;
                LG = 6;
            } else {
                NG = 3;
                LG = 10;
            }
            h = -x;
            k = -y;
            hk = h * k;
            BVN = 0.0;
            if (Math.abs(rho) < 0.925) {
                if (Math.abs(rho) > 0.0) {
                    hs = (h * h + k * k) / 2.0;
                    asr = Math.asin(rho);
                    for (i = 1; i <= LG; i++) {
                        for (ISs = -1; ISs <= 1; ISs += 2) {
                            sn = Math.sin(asr * (ISs * XX[i][NG] + 1.0) / 2.0);
                            BVN = BVN + W[i][NG] * Math.exp((sn * hk - hs) / (1.0 - sn * sn));
                        };
                    };
                    BVN = BVN * asr / (4.0 * Math.PI);
                }
                BVN = BVN + this.CND(-h) * this.CND(-k);
            } else {
                if (rho < 0.0) {
                    k = -k;
                    hk = -hk;
                }
                if (Math.abs(rho) < 1.0) {
                    Ass = (1.0 - rho) * (1.0 + rho);
                    A = Math.sqrt(Ass);
                    bs = Math.pow((h - k), 2.0);
                    c = (4.0 - hk) / 8.0;
                    d = (12.0 - hk) / 16.0;
                    asr = -(bs / Ass + hk) / 2.0;
                    if (asr > -100.0) {
                        BVN = A * Math.exp(asr) * (1.0 - c * (bs - Ass) * (1.0 - d * bs / 5.0) / 3.0 + c * d * Ass * Ass / 5.0);
                    }
                    if (-hk < 100.0) {
                        b = Math.sqrt(bs);
                        BVN = BVN - Math.exp(-hk / 2.0) * Math.sqrt(2.0 * Math.PI) * this.CND(-b / A) * b * (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
                    }
                    A = A / 2.0;
                    for (i = 1; i <= LG; i++) {
                        for (ISs = -1; ISs <= 1; ISs += 2) {
                            xs = Math.pow((A * (ISs * XX[i][NG] + 1.0)), 2.0);
                            rs = Math.sqrt(1.0 - xs);
                            asr = -(bs / xs + hk) / 2.0;
                            if (asr > -100.0) {
                                BVN = BVN + A * W[i][NG] * Math.exp(asr) * (Math.exp(-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) / rs - (1.0 + c * xs * (1.0 + d * xs)));
                            }
                        };
                    };
                    BVN = -BVN / (2.0 * Math.PI);
                }
                if (rho > 0.0) {
                    BVN = BVN + this.CND(-Math.max(h, k));
                } else {
                    BVN = -BVN;
                    if (k > h) {
                        BVN = BVN + this.CND(k) - this.CND(h);
                    }
                }
            }
            return BVN;
        }

        public GBlackScholes(callPutFlag: American.Option, S: number, X: number, T: number, r: number, b: number, v: number): number {
            let functionReturnValue: number = 0.0;
            let d1: number = 0.0;
            let d2: number = 0.0;
            d1 = (Math.log(S / X) + (b + Math.pow(v, 2.0) / 2.0) * T) / (v * Math.sqrt(T));
            d2 = d1 - v * Math.sqrt(T);
            if (American.Option.CALL === callPutFlag) {
                functionReturnValue = S * Math.exp((b - r) * T) * this.CND(d1) - X * Math.exp(-r * T) * this.CND(d2);
            } else if (American.Option.PUT === callPutFlag) {
                functionReturnValue = X * Math.exp(-r * T) * this.CND(-d2) - S * Math.exp((b - r) * T) * this.CND(-d1);
            }
            return functionReturnValue;
        }

        public BSAmericanApprox2002(callPutFlag: American.Option, S: number, X: number, T: number, r: number, b: number, v: number): number {
            let functionReturnValue: number = 0.0;
            if (American.Option.CALL === callPutFlag) {
                functionReturnValue = this.BSAmericanCallApprox2002(S, X, T, r, b, v);
            } else if (American.Option.PUT === callPutFlag) {
                functionReturnValue = this.BSAmericanCallApprox2002(X, S, T, r - b, -b, v);
            }
            return functionReturnValue;
        }

        public BSAmericanCallApprox2002(S: number, X: number, T: number, r: number, b: number, v: number): number {
            let functionReturnValue: number = 0.0;
            let BInfinity: number = 0.0;
            let B0: number = 0.0;
            let ht1: number = 0.0;
            let ht2: number = 0.0;
            let I1: number = 0.0;
            let I2: number = 0.0;
            let alfa1: number = 0.0;
            let alfa2: number = 0.0;
            let Beta: number = 0.0;
            let t1: number = 0.0;
            t1 = 1.0 / 2.0 * (Math.sqrt(5.0) - 1.0) * T;
            if (b >= r) {
                functionReturnValue = this.GBlackScholes(American.Option.CALL, S, X, T, r, b, v);
            } else {
                Beta = (1.0 / 2.0 - b / Math.pow(v, 2.0)) + Math.sqrt(Math.pow((b / Math.pow(v, 2.0) - 1.0 / 2.0), 2.0) + 2.0 * r / Math.pow(v, 2.0));
                BInfinity = Beta / (Beta - 1.0) * X;
                B0 = Math.max(X, r / (r - b) * X);
                ht1 = -(b * t1 + 2.0 * v * Math.sqrt(t1)) * Math.pow(X, 2.0) / ((BInfinity - B0) * B0);
                ht2 = -(b * T + 2.0 * v * Math.sqrt(T)) * Math.pow(X, 2.0) / ((BInfinity - B0) * B0);
                I1 = B0 + (BInfinity - B0) * (1.0 - Math.exp(ht1));
                I2 = B0 + (BInfinity - B0) * (1.0 - Math.exp(ht2));
                alfa1 = (I1 - X) * Math.pow(I1, (-Beta));
                alfa2 = (I2 - X) * Math.pow(I2, (-Beta));
                if (S >= I2) {
                    functionReturnValue = S - X;
                } else {
                    functionReturnValue = alfa2 * Math.pow(S, Beta) - alfa2 * this.phi(S, t1, Beta, I2, I2, r, b, v) + this.phi(S, t1, 1.0, I2, I2, r, b, v) - this.phi(S, t1, 1.0, I1, I2, r, b, v) - X * this.phi(S, t1, 0.0, I2, I2, r, b, v) + X * this.phi(S, t1, 0.0, I1, I2, r, b, v) + alfa1 * this.phi(S, t1, Beta, I1, I2, r, b, v) - alfa1 * this.ksi(S, T, Beta, I1, I2, I1, t1, r, b, v) + this.ksi(S, T, 1.0, I1, I2, I1, t1, r, b, v) - this.ksi(S, T, 1.0, X, I2, I1, t1, r, b, v) - X * this.ksi(S, T, 0.0, I1, I2, I1, t1, r, b, v) + X * this.ksi(S, T, 0.0, X, I2, I1, t1, r, b, v);
                }
            }
            return functionReturnValue;
        }

        phi(S: number, T: number, gamma: number, h: number, i: number, r: number, b: number, v: number): number {
            let lambda: number = 0.0;
            let kappa: number = 0.0;
            let d: number = 0.0;
            lambda = (-r + gamma * b + 0.5 * gamma * (gamma - 1.0) * Math.pow(v, 2.0)) * T;
            d = -(Math.log(S / h) + (b + (gamma - 0.5) * Math.pow(v, 2.0)) * T) / (v * Math.sqrt(T));
            kappa = 2.0 * b / Math.pow(v, 2.0) + 2.0 * gamma - 1.0;
            return Math.exp(lambda) * Math.pow(S, gamma) * (this.CND(d) - Math.pow((i / S), kappa) * this.CND(d - 2.0 * Math.log(i / S) / (v * Math.sqrt(T))));
        }

        public ksi(S: number, T2: number, gamma: number, h: number, I2: number, I1: number, t1: number, r: number, b: number, v: number): number {
            let e1: number = 0.0;
            let e2: number = 0.0;
            let e3: number = 0.0;
            let e4: number = 0.0;
            let f1: number = 0.0;
            let f2: number = 0.0;
            let f3: number = 0.0;
            let f4: number = 0.0;
            let rho: number = 0.0;
            let kappa: number = 0.0;
            let lambda: number = 0.0;
            e1 = (Math.log(S / I1) + (b + (gamma - 0.5) * Math.pow(v, 2.0)) * t1) / (v * Math.sqrt(t1));
            e2 = (Math.log(Math.pow(I2, 2.0) / (S * I1)) + (b + (gamma - 0.5) * Math.pow(v, 2.0)) * t1) / (v * Math.sqrt(t1));
            e3 = (Math.log(S / I1) - (b + (gamma - 0.5) * Math.pow(v, 2.0)) * t1) / (v * Math.sqrt(t1));
            e4 = (Math.log(Math.pow(I2, 2.0) / (S * I1)) - (b + (gamma - 0.5) * Math.pow(v, 2.0)) * t1) / (v * Math.sqrt(t1));
            f1 = (Math.log(S / h) + (b + (gamma - 0.5) * Math.pow(v, 2.0)) * T2) / (v * Math.sqrt(T2));
            f2 = (Math.log(Math.pow(I2, 2.0) / (S * h)) + (b + (gamma - 0.5) * Math.pow(v, 2.0)) * T2) / (v * Math.sqrt(T2));
            f3 = (Math.log(Math.pow(I1, 2.0) / (S * h)) + (b + (gamma - 0.5) * Math.pow(v, 2.0)) * T2) / (v * Math.sqrt(T2));
            f4 = (Math.log(S * Math.pow(I1, 2.0) / (h * Math.pow(I2, 2.0))) + (b + (gamma - 0.5) * Math.pow(v, 2.0)) * T2) / (v * Math.sqrt(T2));
            rho = Math.sqrt(t1 / T2);
            lambda = -r + gamma * b + 0.5 * gamma * (gamma - 1.0) * Math.pow(v, 2.0);
            kappa = 2.0 * b / (Math.pow(v, 2.0)) + (2.0 * gamma - 1.0);
            return Math.exp(lambda * T2) * Math.pow(S, gamma) * (this.CBND(-e1, -f1, rho) - Math.pow((I2 / S), kappa) * this.CBND(-e2, -f2, rho) - Math.pow((I1 / S), kappa) * this.CBND(-e3, -f3, -rho) + Math.pow((I1 / I2), kappa) * this.CBND(-e4, -f4, -rho));
        }

        public GBlackScholesImpVolBisection(callPutFlag: American.Option, S: number, x: number, T: number, r: number, b: number, cm: number): number {
            let vLow: number = 0.0;
            let vHigh: number = 0.0;
            let vi: number = 0.0;
            let cLow: number = 0.0;
            let cHigh: number = 0.0;
            let epsilon: number = 0.0;
            let counter: number = 0;
            vLow = 0.005;
            vHigh = 4.0;
            epsilon = 1.0E-4;
            cLow = this.GBlackScholes(callPutFlag, S, x, T, r, b, vLow);
            cHigh = this.GBlackScholes(callPutFlag, S, x, T, r, b, vHigh);
            counter = 0;
            vi = vLow + (cm - cLow) * (vHigh - vLow) / (cHigh - cLow);
            while ((Math.abs(cm - this.GBlackScholes(callPutFlag, S, x, T, r, b, vi)) > epsilon)) {
                counter = counter + 1;
                if (counter === 1000) {
                    return NaN;
                }
                if (this.GBlackScholes(callPutFlag, S, x, T, r, b, vi) < cm) {
                    vLow = vi;
                } else {
                    vHigh = vi;
                }
                cLow = this.GBlackScholes(callPutFlag, S, x, T, r, b, vLow);
                cHigh = this.GBlackScholes(callPutFlag, S, x, T, r, b, vHigh);
                vi = vLow + (cm - cLow) * (vHigh - vLow) / (cHigh - cLow);
            };
            return vi;
        }

        public BSAmericanApprox2002ImpVolBisection(callPutFlag: American.Option, S: number, x: number, T: number, r: number, b: number, cm: number): number {
            let vLow: number = 0.0;
            let vHigh: number = 0.0;
            let vi: number = 0.0;
            let cLow: number = 0.0;
            let cHigh: number = 0.0;
            let epsilon: number = 0.0;
            let counter: number = 0;
            vLow = 0.005;
            vHigh = 4.0;
            epsilon = 1.0E-4;
            cLow = this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, vLow);
            cHigh = this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, vHigh);
            counter = 0;
            vi = vLow + (cm - cLow) * (vHigh - vLow) / (cHigh - cLow);
            while ((Math.abs(cm - this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, vi)) > epsilon)) {
                counter = counter + 1;
                if (counter === 1000) {
                    return NaN;
                }
                if (this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, vi) < cm) {
                    vLow = vi;
                } else {
                    vHigh = vi;
                }
                cLow = this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, vLow);
                cHigh = this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, vHigh);
                vi = vLow + (cm - cLow) * (vHigh - vLow) / (cHigh - cLow);
            };
            return vi;
        }

        public GBlackScholesNGreeks(outputFlag: American.Greek, callPutFlag: American.Option, S: number, x: number, T: number, r: number, b: number, v: number, dS: number): number {
            if (/* isInfinite */((value) => Number.NEGATIVE_INFINITY === value || Number.POSITIVE_INFINITY === value)(dS) || /* isNaN */isNaN(dS) || dS <= 0.0) {
                dS = 0.01;
            }
            if (null != outputFlag) {
                switch ((outputFlag)) {
                    case Thales.Pricing.American.Greek.DELTA:
                        return (this.GBlackScholes(callPutFlag, S + dS, x, T, r, b, v) - this.GBlackScholes(callPutFlag, S - dS, x, T, r, b, v)) / (2.0 * dS);
                    case Thales.Pricing.American.Greek.ELASTICITY:
                        return (this.GBlackScholes(callPutFlag, S + dS, x, T, r, b, v) - this.GBlackScholes(callPutFlag, S - dS, x, T, r, b, v)) / (2.0 * dS) * S / this.GBlackScholes(callPutFlag, S, x, T, r, b, v);
                    case Thales.Pricing.American.Greek.GAMMA:
                        return (this.GBlackScholes(callPutFlag, S + dS, x, T, r, b, v) - 2.0 * this.GBlackScholes(callPutFlag, S, x, T, r, b, v) + this.GBlackScholes(callPutFlag, S - dS, x, T, r, b, v)) / Math.pow(dS, 2.0);
                    case Thales.Pricing.American.Greek.DGAMMADVOL:
                        return (this.GBlackScholes(callPutFlag, S + dS, x, T, r, b, v + 0.01) - 2.0 * this.GBlackScholes(callPutFlag, S, x, T, r, b, v + 0.01) + this.GBlackScholes(callPutFlag, S - dS, x, T, r, b, v + 0.01) - this.GBlackScholes(callPutFlag, S + dS, x, T, r, b, v - 0.01) + 2.0 * this.GBlackScholes(callPutFlag, S, x, T, r, b, v - 0.01) - this.GBlackScholes(callPutFlag, S - dS, x, T, r, b, v - 0.01)) / (2.0 * 0.01 * Math.pow(dS, 2.0)) / 100.0;
                    case Thales.Pricing.American.Greek.GAMMAP:
                        return S / 100.0 * (this.GBlackScholes(callPutFlag, S + dS, x, T, r, b, v) - 2.0 * this.GBlackScholes(callPutFlag, S, x, T, r, b, v) + this.GBlackScholes(callPutFlag, S - dS, x, T, r, b, v)) / Math.pow(dS, 2.0);
                    case Thales.Pricing.American.Greek.VEGA:
                        return (this.GBlackScholes(callPutFlag, S, x, T, r, b, v + 0.01) - this.GBlackScholes(callPutFlag, S, x, T, r, b, v - 0.01)) / 2.0;
                    case Thales.Pricing.American.Greek.DVEGADVOL:
                        return (this.GBlackScholes(callPutFlag, S, x, T, r, b, v + 0.01) - 2.0 * this.GBlackScholes(callPutFlag, S, x, T, r, b, v) + this.GBlackScholes(callPutFlag, S, x, T, r, b, v - 0.01));
                    case Thales.Pricing.American.Greek.VEGAP:
                        return v / 0.1 * (this.GBlackScholes(callPutFlag, S, x, T, r, b, v + 0.01) - this.GBlackScholes(callPutFlag, S, x, T, r, b, v - 0.01)) / 2.0;
                    case Thales.Pricing.American.Greek.THETA:
                        if (T <= 1.0 / 365.0) {
                            return this.GBlackScholes(callPutFlag, S, x, 1.0E-5, r, b, v) - this.GBlackScholes(callPutFlag, S, x, T, r, b, v);
                        } else {
                            return this.GBlackScholes(callPutFlag, S, x, T - 1.0 / 365.0, r, b, v) - this.GBlackScholes(callPutFlag, S, x, T, r, b, v);
                        }
                    case Thales.Pricing.American.Greek.RHO:
                        return (this.GBlackScholes(callPutFlag, S, x, T, r + 0.01, b + 0.01, v) - this.GBlackScholes(callPutFlag, S, x, T, r - 0.01, b - 0.01, v)) / (2.0);
                    case Thales.Pricing.American.Greek.RHO_FUTURES_OPTION:
                        return (this.GBlackScholes(callPutFlag, S, x, T, r + 0.01, 0.0, v) - this.GBlackScholes(callPutFlag, S, x, T, r - 0.01, 0.0, v)) / (2.0);
                    case Thales.Pricing.American.Greek.PHI_RHO2:
                        return (this.GBlackScholes(callPutFlag, S, x, T, r, b - 0.01, v) - this.GBlackScholes(callPutFlag, S, x, T, r, b + 0.01, v)) / (2.0);
                    case Thales.Pricing.American.Greek.CARRY_RHO:
                        return (this.GBlackScholes(callPutFlag, S, x, T, r, b + 0.01, v) - this.GBlackScholes(callPutFlag, S, x, T, r, b - 0.01, v)) / (2.0);
                    case Thales.Pricing.American.Greek.DDELTADVOL:
                        return 1.0 / (4.0 * dS * 0.01) * (this.GBlackScholes(callPutFlag, S + dS, x, T, r, b, v + 0.01) - this.GBlackScholes(callPutFlag, S + dS, x, T, r, b, v - 0.01) - this.GBlackScholes(callPutFlag, S - dS, x, T, r, b, v + 0.01) + this.GBlackScholes(callPutFlag, S - dS, x, T, r, b, v - 0.01)) / 100.0;
                    case Thales.Pricing.American.Greek.STRIKE_DELTA:
                        return (this.GBlackScholes(callPutFlag, S, x + dS, T, r, b, v) - this.GBlackScholes(callPutFlag, S, x - dS, T, r, b, v)) / (2.0 * dS);
                    case Thales.Pricing.American.Greek.SPEED:
                        return 1.0 / Math.pow(dS, 3.0) * (this.GBlackScholes(callPutFlag, S + 2.0 * dS, x, T, r, b, v) - 3.0 * this.GBlackScholes(callPutFlag, S + dS, x, T, r, b, v) + 3.0 * this.GBlackScholes(callPutFlag, S, x, T, r, b, v) - this.GBlackScholes(callPutFlag, S - dS, x, T, r, b, v));
                    default:
                        break;
                }
            }
            return NaN;
        }

        public BSAmericanApprox2002NGreeks(outputFlag: American.Greek, callPutFlag: American.Option, S: number, x: number, T: number, r: number, b: number, v: number, dS: number): number {
            if (/* isInfinite */((value) => Number.NEGATIVE_INFINITY === value || Number.POSITIVE_INFINITY === value)(dS) || /* isNaN */isNaN(dS) || dS <= 0.0) {
                dS = 0.01;
            }
            if (null != outputFlag) {
                switch ((outputFlag)) {
                    case Thales.Pricing.American.Greek.DELTA:
                        return (this.BSAmericanApprox2002(callPutFlag, S + dS, x, T, r, b, v) - this.BSAmericanApprox2002(callPutFlag, S - dS, x, T, r, b, v)) / (2.0 * dS);
                    case Thales.Pricing.American.Greek.ELASTICITY:
                        return (this.BSAmericanApprox2002(callPutFlag, S + dS, x, T, r, b, v) - this.BSAmericanApprox2002(callPutFlag, S - dS, x, T, r, b, v)) / (2.0 * dS) * S / this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v);
                    case Thales.Pricing.American.Greek.GAMMA:
                        return (this.BSAmericanApprox2002(callPutFlag, S + dS, x, T, r, b, v) - 2.0 * this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v) + this.BSAmericanApprox2002(callPutFlag, S - dS, x, T, r, b, v)) / Math.pow(dS, 2.0);
                    case Thales.Pricing.American.Greek.DGAMMADVOL:
                        return (this.BSAmericanApprox2002(callPutFlag, S + dS, x, T, r, b, v + 0.01) - 2.0 * this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v + 0.01) + this.BSAmericanApprox2002(callPutFlag, S - dS, x, T, r, b, v + 0.01) - this.BSAmericanApprox2002(callPutFlag, S + dS, x, T, r, b, v - 0.01) + 2.0 * this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v - 0.01) - this.BSAmericanApprox2002(callPutFlag, S - dS, x, T, r, b, v - 0.01)) / (2.0 * 0.01 * Math.pow(dS, 2.0)) / 100.0;
                    case Thales.Pricing.American.Greek.GAMMAP:
                        return S / 100.0 * (this.BSAmericanApprox2002(callPutFlag, S + dS, x, T, r, b, v) - 2.0 * this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v) + this.BSAmericanApprox2002(callPutFlag, S - dS, x, T, r, b, v)) / Math.pow(dS, 2.0);
                    case Thales.Pricing.American.Greek.VEGA:
                        return (this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v + 0.01) - this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v - 0.01)) / 2.0;
                    case Thales.Pricing.American.Greek.DVEGADVOL:
                        return (this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v + 0.01) - 2.0 * this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v) + this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v - 0.01));
                    case Thales.Pricing.American.Greek.VEGAP:
                        return v / 0.1 * (this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v + 0.01) - this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v - 0.01)) / 2.0;
                    case Thales.Pricing.American.Greek.THETA:
                        if (T <= 1.0 / 365.0) {
                            return this.BSAmericanApprox2002(callPutFlag, S, x, 1.0E-5, r, b, v) - this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v);
                        } else {
                            return this.BSAmericanApprox2002(callPutFlag, S, x, T - 1.0 / 365.0, r, b, v) - this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v);
                        }
                    case Thales.Pricing.American.Greek.RHO:
                        return (this.BSAmericanApprox2002(callPutFlag, S, x, T, r + 0.01, b + 0.01, v) - this.BSAmericanApprox2002(callPutFlag, S, x, T, r - 0.01, b - 0.01, v)) / (2.0);
                    case Thales.Pricing.American.Greek.RHO_FUTURES_OPTION:
                        return (this.BSAmericanApprox2002(callPutFlag, S, x, T, r + 0.01, 0.0, v) - this.BSAmericanApprox2002(callPutFlag, S, x, T, r - 0.01, 0.0, v)) / (2.0);
                    case Thales.Pricing.American.Greek.PHI_RHO2:
                        return (this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b - 0.01, v) - this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b + 0.01, v)) / (2.0);
                    case Thales.Pricing.American.Greek.CARRY_RHO:
                        return (this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b + 0.01, v) - this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b - 0.01, v)) / (2.0);
                    case Thales.Pricing.American.Greek.DDELTADVOL:
                        return 1.0 / (4.0 * dS * 0.01) * (this.BSAmericanApprox2002(callPutFlag, S + dS, x, T, r, b, v + 0.01) - this.BSAmericanApprox2002(callPutFlag, S + dS, x, T, r, b, v - 0.01) - this.BSAmericanApprox2002(callPutFlag, S - dS, x, T, r, b, v + 0.01) + this.BSAmericanApprox2002(callPutFlag, S - dS, x, T, r, b, v - 0.01)) / 100.0;
                    case Thales.Pricing.American.Greek.STRIKE_DELTA:
                        return (this.BSAmericanApprox2002(callPutFlag, S, x + dS, T, r, b, v) - this.BSAmericanApprox2002(callPutFlag, S, x - dS, T, r, b, v)) / (2.0 * dS);
                    case Thales.Pricing.American.Greek.SPEED:
                        return 1.0 / Math.pow(dS, 3.0) * (this.BSAmericanApprox2002(callPutFlag, S + 2.0 * dS, x, T, r, b, v) - 3.0 * this.BSAmericanApprox2002(callPutFlag, S + dS, x, T, r, b, v) + 3.0 * this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, v) - this.BSAmericanApprox2002(callPutFlag, S - dS, x, T, r, b, v));
                    default:
                        break;
                }
            }
            return NaN;
        }

        public testBSAmericanApprox2002(): boolean {
            let expectations: number[] = [0.0205, 1.8757, 10.0, 0.3151, 3.1256, 10.3725, 0.9479, 4.3746, 11.1578, 0.8099, 4.0628, 10.7898, 2.718, 6.7661, 12.9814, 4.9665, 9.4608, 15.5137, 10.0, 1.8757, 0.0408, 10.228, 3.1256, 0.4552, 10.8663, 4.3746, 1.2383, 10.54, 4.0628, 1.0689, 12.4097, 6.7661, 3.2932, 14.6445, 9.4608, 5.8374];
            let callPutFlags: American.Option[] = [American.Option.CALL, American.Option.PUT];
            let X: number = 100.0;
            let r: number = 0.1;
            let b: number = 0.0;
            let Ss: number[] = [90.0, 100.0, 110.0];
            let Ts: number[] = [0.1, 0.5];
            let vs: number[] = [0.15, 0.25, 0.35];
            let i: number = 0;
            for (let index612 = 0; index612 < callPutFlags.length; index612++) {
                let callPutFlag = callPutFlags[index612];
                {
                    for (let index613 = 0; index613 < Ts.length; index613++) {
                        let T = Ts[index613];
                        {
                            for (let index614 = 0; index614 < vs.length; index614++) {
                                let v = vs[index614];
                                {
                                    for (let index615 = 0; index615 < Ss.length; index615++) {
                                        let S = Ss[index615];
                                        {
                                            let result: number = this.BSAmericanApprox2002(callPutFlag, S, X, T, r, b, v);
                                            if (expectations[i] - result > 9.0E-5) {
                                                return false;
                                            }
                                            i++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return true;
        }

        public testGBlackScholesImpVolBisection(): boolean {
            let callPutFlags: American.Option[] = [American.Option.CALL, American.Option.PUT];
            let X: number = 100.0;
            let r: number = 0.1;
            let b: number = 0.0;
            let Ss: number[] = [90.0, 100.0, 110.0];
            let Ts: number[] = [0.1, 0.5];
            let vs: number[] = [0.15, 0.25, 0.35];
            for (let index616 = 0; index616 < callPutFlags.length; index616++) {
                let callPutFlag = callPutFlags[index616];
                {
                    for (let index617 = 0; index617 < Ts.length; index617++) {
                        let T = Ts[index617];
                        {
                            for (let index618 = 0; index618 < vs.length; index618++) {
                                let v = vs[index618];
                                {
                                    for (let index619 = 0; index619 < Ss.length; index619++) {
                                        let S = Ss[index619];
                                        {
                                            let price: number = this.GBlackScholes(callPutFlag, S, X, T, r, b, v);
                                            let impliedVolatility: number = this.GBlackScholesImpVolBisection(callPutFlag, S, X, T, r, b, price);
                                            if (v - impliedVolatility > 9.0E-4) {
                                                return false;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return true;
        }

        public testBSAmericanApprox2002ImpVolBisection(): boolean {
            let callPutFlags: American.Option[] = [American.Option.CALL, American.Option.PUT];
            let X: number = 100.0;
            let r: number = 0.1;
            let b: number = 0.0;
            let Ss: number[] = [90.0, 100.0, 110.0];
            let Ts: number[] = [0.1, 0.5];
            let vs: number[] = [0.15, 0.25, 0.35];
            for (let index616 = 0; index616 < callPutFlags.length; index616++) {
                let callPutFlag = callPutFlags[index616];
                {
                    for (let index617 = 0; index617 < Ts.length; index617++) {
                        let T = Ts[index617];
                        {
                            for (let index618 = 0; index618 < vs.length; index618++) {
                                let v = vs[index618];
                                {
                                    for (let index619 = 0; index619 < Ss.length; index619++) {
                                        let S = Ss[index619];
                                        {
                                            let price: number = this.BSAmericanApprox2002(callPutFlag, S, X, T, r, b, v);
                                            let impliedVolatility: number = this.BSAmericanApprox2002ImpVolBisection(callPutFlag, S, X, T, r, b, price);
                                            if (v - impliedVolatility > 9.0E-4 && Math.abs(S - X) !== price) {
                                                return false;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return true;
        }

        public testGBlackScholesNGreeks(): boolean {
            let greeks: American.Greek[] = [American.Greek.DELTA, American.Greek.ELASTICITY, American.Greek.GAMMA, American.Greek.DGAMMADVOL, American.Greek.GAMMAP, American.Greek.VEGA, American.Greek.DVEGADVOL, American.Greek.VEGAP, American.Greek.THETA, American.Greek.RHO, American.Greek.RHO_FUTURES_OPTION, American.Greek.PHI_RHO2, American.Greek.CARRY_RHO, American.Greek.DDELTADVOL, American.Greek.STRIKE_DELTA, American.Greek.SPEED];
            let expectations: number[] = [0.503105, 9.059951, 0.026794, -8.97E-4, 0.026258, 0.192999, -1.9E-5, 0.578997, -0.03708, 0.109655, -0.012124, -0.123262, 0.123262, 0.00166, -0.438623, -3.17E-4];
            let i: number = 0;
            for (let index620 = 0; index620 < greeks.length; index620++) {
                let greek = greeks[index620];
                {
                    let value: number = this.GBlackScholesNGreeks(greek, American.Option.CALL, 98.0, 100.0, 0.25, 0.1, 0.05, 0.3, -1);
                    if (expectations[i] - value > 9.0E-4) {
                        return false;
                    }
                    i++;
                }
            }
            return true;
        }

        public testBSAmericanApprox2002NGreeks(): boolean {
            let greeks: American.Greek[] = [American.Greek.DELTA, American.Greek.ELASTICITY, American.Greek.GAMMA, American.Greek.DGAMMADVOL, American.Greek.GAMMAP, American.Greek.VEGA, American.Greek.DVEGADVOL, American.Greek.VEGAP, American.Greek.THETA, American.Greek.RHO, American.Greek.RHO_FUTURES_OPTION, American.Greek.PHI_RHO2, American.Greek.CARRY_RHO, American.Greek.DDELTADVOL, American.Greek.STRIKE_DELTA, American.Greek.SPEED];
            let expectations: number[] = [0.503105, 9.059951, 0.026794, -8.97E-4, 0.026258, 0.192999, -1.9E-5, 0.578997, -0.03708, 0.109655, -0.012124, -0.123262, 0.123262, 0.00166, -0.438623, -3.17E-4];
            let i: number = 0;
            for (let index620 = 0; index620 < greeks.length; index620++) {
                let greek = greeks[index620];
                {
                    let value: number = this.BSAmericanApprox2002NGreeks(greek, American.Option.CALL, 98.0, 100.0, 0.25, 0.1, 0.05, 0.3, -1);
                    if (expectations[i] - value > 9.0E-4) {
                        return false;
                    }
                    i++;
                }
            }
            return true;
        }

        public static $inject: string[] = [];
    }
    app.service("Thales.Pricing.American", American);
    
    American["__class"] = "Thales.Pricing.American";
    export namespace American {

        export enum Option {
            CALL, PUT
        }

        export enum Greek {
            DELTA, ELASTICITY, GAMMA, DGAMMADVOL, GAMMAP, VEGA, DVEGADVOL, VEGAP, THETA, RHO, RHO_FUTURES_OPTION, PHI_RHO2, CARRY_RHO, DDELTADVOL, STRIKE_DELTA, SPEED
        }
    }
}

