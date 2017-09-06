"use strict";
var Thales;
(function (Thales) {
    angular.module("thales", ["ngRoute"]);
    Thales.getModule = function () {
        return angular.module("thales");
    };
})(Thales || (Thales = {}));
var Thales;
(function (Thales) {
    var Pricing;
    (function (Pricing) {
        var app = Thales.getModule();
        var American = (function () {
            function American() {
            }
            American.prototype.CND = function (X) {
                var functionReturnValue = 0.0;
                var L = 0.0;
                var K = 0.0;
                var a1 = 0.31938153;
                var a2 = -0.356563782;
                var a3 = 1.781477937;
                var a4 = -1.821255978;
                var a5 = 1.330274429;
                L = Math.abs(X);
                K = 1.0 / (1.0 + 0.2316419 * L);
                functionReturnValue = 1.0 - 1.0 / Math.sqrt(2.0 * Math.PI) * Math.exp(-Math.pow(L, 2.0) / 2.0) * (a1 * K + a2 * Math.pow(K, 2.0) + a3 * Math.pow(K, 3.0) + a4 * Math.pow(K, 4.0) + a5 * Math.pow(K, 5.0));
                if (X < 0.0) {
                    functionReturnValue = 1.0 - functionReturnValue;
                }
                return functionReturnValue;
            };
            American.prototype.CBND = function (x, y, rho) {
                var i = 0;
                var ISs = 0;
                var LG = 0;
                var NG = 0;
                var XX = (function (dims) { var allocate = function (dims) { if (dims.length == 0) {
                    return 0;
                }
                else {
                    var array = [];
                    for (var i_1 = 0; i_1 < dims[0]; i_1++) {
                        array.push(allocate(dims.slice(1)));
                    }
                    return array;
                } }; return allocate(dims); })([11, 4]);
                var W = (function (dims) { var allocate = function (dims) { if (dims.length == 0) {
                    return 0;
                }
                else {
                    var array = [];
                    for (var i_2 = 0; i_2 < dims[0]; i_2++) {
                        array.push(allocate(dims.slice(1)));
                    }
                    return array;
                } }; return allocate(dims); })([11, 4]);
                var h = 0.0;
                var k = 0.0;
                var hk = 0.0;
                var hs = 0.0;
                var BVN = 0.0;
                var Ass = 0.0;
                var asr = 0.0;
                var sn = 0.0;
                var A = 0.0;
                var b = 0.0;
                var bs = 0.0;
                var c = 0.0;
                var d = 0.0;
                var xs = 0.0;
                var rs = 0.0;
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
                }
                else if (Math.abs(rho) < 0.75) {
                    NG = 2;
                    LG = 6;
                }
                else {
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
                            }
                            ;
                        }
                        ;
                        BVN = BVN * asr / (4.0 * Math.PI);
                    }
                    BVN = BVN + this.CND(-h) * this.CND(-k);
                }
                else {
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
                            }
                            ;
                        }
                        ;
                        BVN = -BVN / (2.0 * Math.PI);
                    }
                    if (rho > 0.0) {
                        BVN = BVN + this.CND(-Math.max(h, k));
                    }
                    else {
                        BVN = -BVN;
                        if (k > h) {
                            BVN = BVN + this.CND(k) - this.CND(h);
                        }
                    }
                }
                return BVN;
            };
            American.prototype.GBlackScholes = function (callPutFlag, S, X, T, r, b, v) {
                var functionReturnValue = 0.0;
                var d1 = 0.0;
                var d2 = 0.0;
                d1 = (Math.log(S / X) + (b + Math.pow(v, 2.0) / 2.0) * T) / (v * Math.sqrt(T));
                d2 = d1 - v * Math.sqrt(T);
                if (American.Option.CALL === callPutFlag) {
                    functionReturnValue = S * Math.exp((b - r) * T) * this.CND(d1) - X * Math.exp(-r * T) * this.CND(d2);
                }
                else if (American.Option.PUT === callPutFlag) {
                    functionReturnValue = X * Math.exp(-r * T) * this.CND(-d2) - S * Math.exp((b - r) * T) * this.CND(-d1);
                }
                return functionReturnValue;
            };
            American.prototype.BSAmericanApprox2002 = function (callPutFlag, S, X, T, r, b, v) {
                var functionReturnValue = 0.0;
                if (American.Option.CALL === callPutFlag) {
                    functionReturnValue = this.BSAmericanCallApprox2002(S, X, T, r, b, v);
                }
                else if (American.Option.PUT === callPutFlag) {
                    functionReturnValue = this.BSAmericanCallApprox2002(X, S, T, r - b, -b, v);
                }
                return functionReturnValue;
            };
            American.prototype.BSAmericanCallApprox2002 = function (S, X, T, r, b, v) {
                var functionReturnValue = 0.0;
                var BInfinity = 0.0;
                var B0 = 0.0;
                var ht1 = 0.0;
                var ht2 = 0.0;
                var I1 = 0.0;
                var I2 = 0.0;
                var alfa1 = 0.0;
                var alfa2 = 0.0;
                var Beta = 0.0;
                var t1 = 0.0;
                t1 = 1.0 / 2.0 * (Math.sqrt(5.0) - 1.0) * T;
                if (b >= r) {
                    functionReturnValue = this.GBlackScholes(American.Option.CALL, S, X, T, r, b, v);
                }
                else {
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
                    }
                    else {
                        functionReturnValue = alfa2 * Math.pow(S, Beta) - alfa2 * this.phi(S, t1, Beta, I2, I2, r, b, v) + this.phi(S, t1, 1.0, I2, I2, r, b, v) - this.phi(S, t1, 1.0, I1, I2, r, b, v) - X * this.phi(S, t1, 0.0, I2, I2, r, b, v) + X * this.phi(S, t1, 0.0, I1, I2, r, b, v) + alfa1 * this.phi(S, t1, Beta, I1, I2, r, b, v) - alfa1 * this.ksi(S, T, Beta, I1, I2, I1, t1, r, b, v) + this.ksi(S, T, 1.0, I1, I2, I1, t1, r, b, v) - this.ksi(S, T, 1.0, X, I2, I1, t1, r, b, v) - X * this.ksi(S, T, 0.0, I1, I2, I1, t1, r, b, v) + X * this.ksi(S, T, 0.0, X, I2, I1, t1, r, b, v);
                    }
                }
                return functionReturnValue;
            };
            American.prototype.phi = function (S, T, gamma, h, i, r, b, v) {
                var lambda = 0.0;
                var kappa = 0.0;
                var d = 0.0;
                lambda = (-r + gamma * b + 0.5 * gamma * (gamma - 1.0) * Math.pow(v, 2.0)) * T;
                d = -(Math.log(S / h) + (b + (gamma - 0.5) * Math.pow(v, 2.0)) * T) / (v * Math.sqrt(T));
                kappa = 2.0 * b / Math.pow(v, 2.0) + 2.0 * gamma - 1.0;
                return Math.exp(lambda) * Math.pow(S, gamma) * (this.CND(d) - Math.pow((i / S), kappa) * this.CND(d - 2.0 * Math.log(i / S) / (v * Math.sqrt(T))));
            };
            American.prototype.ksi = function (S, T2, gamma, h, I2, I1, t1, r, b, v) {
                var e1 = 0.0;
                var e2 = 0.0;
                var e3 = 0.0;
                var e4 = 0.0;
                var f1 = 0.0;
                var f2 = 0.0;
                var f3 = 0.0;
                var f4 = 0.0;
                var rho = 0.0;
                var kappa = 0.0;
                var lambda = 0.0;
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
            };
            American.prototype.GBlackScholesImpVolBisection = function (callPutFlag, S, x, T, r, b, cm) {
                var vLow = 0.0;
                var vHigh = 0.0;
                var vi = 0.0;
                var cLow = 0.0;
                var cHigh = 0.0;
                var epsilon = 0.0;
                var counter = 0;
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
                    }
                    else {
                        vHigh = vi;
                    }
                    cLow = this.GBlackScholes(callPutFlag, S, x, T, r, b, vLow);
                    cHigh = this.GBlackScholes(callPutFlag, S, x, T, r, b, vHigh);
                    vi = vLow + (cm - cLow) * (vHigh - vLow) / (cHigh - cLow);
                }
                ;
                return vi;
            };
            American.prototype.BSAmericanApprox2002ImpVolBisection = function (callPutFlag, S, x, T, r, b, cm) {
                var vLow = 0.0;
                var vHigh = 0.0;
                var vi = 0.0;
                var cLow = 0.0;
                var cHigh = 0.0;
                var epsilon = 0.0;
                var counter = 0;
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
                    }
                    else {
                        vHigh = vi;
                    }
                    cLow = this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, vLow);
                    cHigh = this.BSAmericanApprox2002(callPutFlag, S, x, T, r, b, vHigh);
                    vi = vLow + (cm - cLow) * (vHigh - vLow) / (cHigh - cLow);
                }
                ;
                return vi;
            };
            American.prototype.GBlackScholesNGreeks = function (outputFlag, callPutFlag, S, x, T, r, b, v, dS) {
                if ((function (value) { return Number.NEGATIVE_INFINITY === value || Number.POSITIVE_INFINITY === value; })(dS) || isNaN(dS) || dS <= 0.0) {
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
                            }
                            else {
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
            };
            American.prototype.BSAmericanApprox2002NGreeks = function (outputFlag, callPutFlag, S, x, T, r, b, v, dS) {
                if ((function (value) { return Number.NEGATIVE_INFINITY === value || Number.POSITIVE_INFINITY === value; })(dS) || isNaN(dS) || dS <= 0.0) {
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
                            }
                            else {
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
            };
            American.prototype.testBSAmericanApprox2002 = function () {
                var expectations = [0.0205, 1.8757, 10.0, 0.3151, 3.1256, 10.3725, 0.9479, 4.3746, 11.1578, 0.8099, 4.0628, 10.7898, 2.718, 6.7661, 12.9814, 4.9665, 9.4608, 15.5137, 10.0, 1.8757, 0.0408, 10.228, 3.1256, 0.4552, 10.8663, 4.3746, 1.2383, 10.54, 4.0628, 1.0689, 12.4097, 6.7661, 3.2932, 14.6445, 9.4608, 5.8374];
                var callPutFlags = [American.Option.CALL, American.Option.PUT];
                var X = 100.0;
                var r = 0.1;
                var b = 0.0;
                var Ss = [90.0, 100.0, 110.0];
                var Ts = [0.1, 0.5];
                var vs = [0.15, 0.25, 0.35];
                var i = 0;
                for (var index612 = 0; index612 < callPutFlags.length; index612++) {
                    var callPutFlag = callPutFlags[index612];
                    {
                        for (var index613 = 0; index613 < Ts.length; index613++) {
                            var T = Ts[index613];
                            {
                                for (var index614 = 0; index614 < vs.length; index614++) {
                                    var v = vs[index614];
                                    {
                                        for (var index615 = 0; index615 < Ss.length; index615++) {
                                            var S = Ss[index615];
                                            {
                                                var result = this.BSAmericanApprox2002(callPutFlag, S, X, T, r, b, v);
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
            };
            American.prototype.testGBlackScholesImpVolBisection = function () {
                var callPutFlags = [American.Option.CALL, American.Option.PUT];
                var X = 100.0;
                var r = 0.1;
                var b = 0.0;
                var Ss = [90.0, 100.0, 110.0];
                var Ts = [0.1, 0.5];
                var vs = [0.15, 0.25, 0.35];
                for (var index616 = 0; index616 < callPutFlags.length; index616++) {
                    var callPutFlag = callPutFlags[index616];
                    {
                        for (var index617 = 0; index617 < Ts.length; index617++) {
                            var T = Ts[index617];
                            {
                                for (var index618 = 0; index618 < vs.length; index618++) {
                                    var v = vs[index618];
                                    {
                                        for (var index619 = 0; index619 < Ss.length; index619++) {
                                            var S = Ss[index619];
                                            {
                                                var price = this.GBlackScholes(callPutFlag, S, X, T, r, b, v);
                                                var impliedVolatility = this.GBlackScholesImpVolBisection(callPutFlag, S, X, T, r, b, price);
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
            };
            American.prototype.testBSAmericanApprox2002ImpVolBisection = function () {
                var callPutFlags = [American.Option.CALL, American.Option.PUT];
                var X = 100.0;
                var r = 0.1;
                var b = 0.0;
                var Ss = [90.0, 100.0, 110.0];
                var Ts = [0.1, 0.5];
                var vs = [0.15, 0.25, 0.35];
                for (var index616 = 0; index616 < callPutFlags.length; index616++) {
                    var callPutFlag = callPutFlags[index616];
                    {
                        for (var index617 = 0; index617 < Ts.length; index617++) {
                            var T = Ts[index617];
                            {
                                for (var index618 = 0; index618 < vs.length; index618++) {
                                    var v = vs[index618];
                                    {
                                        for (var index619 = 0; index619 < Ss.length; index619++) {
                                            var S = Ss[index619];
                                            {
                                                var price = this.BSAmericanApprox2002(callPutFlag, S, X, T, r, b, v);
                                                var impliedVolatility = this.BSAmericanApprox2002ImpVolBisection(callPutFlag, S, X, T, r, b, price);
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
            };
            American.prototype.testGBlackScholesNGreeks = function () {
                var greeks = [American.Greek.DELTA, American.Greek.ELASTICITY, American.Greek.GAMMA, American.Greek.DGAMMADVOL, American.Greek.GAMMAP, American.Greek.VEGA, American.Greek.DVEGADVOL, American.Greek.VEGAP, American.Greek.THETA, American.Greek.RHO, American.Greek.RHO_FUTURES_OPTION, American.Greek.PHI_RHO2, American.Greek.CARRY_RHO, American.Greek.DDELTADVOL, American.Greek.STRIKE_DELTA, American.Greek.SPEED];
                var expectations = [0.503105, 9.059951, 0.026794, -8.97E-4, 0.026258, 0.192999, -1.9E-5, 0.578997, -0.03708, 0.109655, -0.012124, -0.123262, 0.123262, 0.00166, -0.438623, -3.17E-4];
                var i = 0;
                for (var index620 = 0; index620 < greeks.length; index620++) {
                    var greek = greeks[index620];
                    {
                        var value = this.GBlackScholesNGreeks(greek, American.Option.CALL, 98.0, 100.0, 0.25, 0.1, 0.05, 0.3, -1);
                        if (expectations[i] - value > 9.0E-4) {
                            return false;
                        }
                        i++;
                    }
                }
                return true;
            };
            American.prototype.testBSAmericanApprox2002NGreeks = function () {
                var greeks = [American.Greek.DELTA, American.Greek.ELASTICITY, American.Greek.GAMMA, American.Greek.DGAMMADVOL, American.Greek.GAMMAP, American.Greek.VEGA, American.Greek.DVEGADVOL, American.Greek.VEGAP, American.Greek.THETA, American.Greek.RHO, American.Greek.RHO_FUTURES_OPTION, American.Greek.PHI_RHO2, American.Greek.CARRY_RHO, American.Greek.DDELTADVOL, American.Greek.STRIKE_DELTA, American.Greek.SPEED];
                var expectations = [0.503105, 9.059951, 0.026794, -8.97E-4, 0.026258, 0.192999, -1.9E-5, 0.578997, -0.03708, 0.109655, -0.012124, -0.123262, 0.123262, 0.00166, -0.438623, -3.17E-4];
                var i = 0;
                for (var index620 = 0; index620 < greeks.length; index620++) {
                    var greek = greeks[index620];
                    {
                        var value = this.BSAmericanApprox2002NGreeks(greek, American.Option.CALL, 98.0, 100.0, 0.25, 0.1, 0.05, 0.3, -1);
                        if (expectations[i] - value > 9.0E-4) {
                            return false;
                        }
                        i++;
                    }
                }
                return true;
            };
            American.$inject = [];
            return American;
        }());
        Pricing.American = American;
        app.service("Thales.Pricing.American", American);
        American["__class"] = "Thales.Pricing.American";
        (function (American) {
            var Option;
            (function (Option) {
                Option[Option["CALL"] = 0] = "CALL";
                Option[Option["PUT"] = 1] = "PUT";
            })(Option = American.Option || (American.Option = {}));
            var Greek;
            (function (Greek) {
                Greek[Greek["DELTA"] = 0] = "DELTA";
                Greek[Greek["ELASTICITY"] = 1] = "ELASTICITY";
                Greek[Greek["GAMMA"] = 2] = "GAMMA";
                Greek[Greek["DGAMMADVOL"] = 3] = "DGAMMADVOL";
                Greek[Greek["GAMMAP"] = 4] = "GAMMAP";
                Greek[Greek["VEGA"] = 5] = "VEGA";
                Greek[Greek["DVEGADVOL"] = 6] = "DVEGADVOL";
                Greek[Greek["VEGAP"] = 7] = "VEGAP";
                Greek[Greek["THETA"] = 8] = "THETA";
                Greek[Greek["RHO"] = 9] = "RHO";
                Greek[Greek["RHO_FUTURES_OPTION"] = 10] = "RHO_FUTURES_OPTION";
                Greek[Greek["PHI_RHO2"] = 11] = "PHI_RHO2";
                Greek[Greek["CARRY_RHO"] = 12] = "CARRY_RHO";
                Greek[Greek["DDELTADVOL"] = 13] = "DDELTADVOL";
                Greek[Greek["STRIKE_DELTA"] = 14] = "STRIKE_DELTA";
                Greek[Greek["SPEED"] = 15] = "SPEED";
            })(Greek = American.Greek || (American.Greek = {}));
        })(American = Pricing.American || (Pricing.American = {}));
    })(Pricing = Thales.Pricing || (Thales.Pricing = {}));
})(Thales || (Thales = {}));
var Thales;
(function (Thales) {
    var Yahoo;
    (function (Yahoo) {
        var app = Thales.getModule();
        var Download = (function () {
            function Download($http, $q, $sce) {
                this.$http = $http;
                this.$q = $q;
                this.$sce = $sce;
            }
            Download.prototype.download = function (symbol) {
                var _this = this;
                var results = this.$http.jsonp(this.$sce.trustAsResourceUrl("https://query1.finance.yahoo.com/v7/finance/options/" + symbol))
                    .then(function (promiseValue) {
                    var rootObject = promiseValue.data;
                    if (rootObject.optionChain.error) {
                        rootObject.optionChain.error = "Error A ~ " + JSON.stringify(rootObject.optionChain.error);
                        return [rootObject.optionChain];
                    }
                    else {
                        var expirations = rootObject.optionChain.result[0].expirationDates;
                        var resultPromises = [];
                        angular.forEach(expirations, function (expiration, key) {
                            resultPromises.push(_this.downloadByExpiration(symbol, expiration));
                        });
                        var onePromise = _this.$q.all(resultPromises);
                        return onePromise;
                    }
                }).catch(function (reason) {
                    return [{
                            error: "Error C ~ " + JSON.stringify(reason),
                            result: undefined
                        }];
                });
                return results;
            };
            Download.prototype.downloadByExpiration = function (symbol, expiration) {
                var results = this.$http.jsonp(this.$sce.trustAsResourceUrl("https://query1.finance.yahoo.com/v7/finance/options/" + symbol + "?date=" + expiration))
                    .then(function (promiseValue) {
                    return promiseValue.data.optionChain;
                }).catch(function (reason) {
                    return {
                        error: "Error D ~ " + JSON.stringify(reason),
                        result: undefined
                    };
                });
                return results;
            };
            Download.prototype.testDownload = function () {
                var _this = this;
                var promise = this.download("AAPL");
                promise.then(function (promiseValue) {
                    var result = promiseValue.every(function (value, index) {
                        if (value.error) {
                            console.error(value.error);
                        }
                        return !value.error;
                    });
                    _this.testDownloadStatus = result;
                });
                this.testDownloadStatus = undefined;
            };
            Download.$inject = ["$http", "$q", "$sce"];
            return Download;
        }());
        Yahoo.Download = Download;
        app.service("Thales.Yahoo.Download", Download);
    })(Yahoo = Thales.Yahoo || (Thales.Yahoo = {}));
})(Thales || (Thales = {}));
var Thales;
(function (Thales) {
    var Controllers;
    (function (Controllers) {
        var app = Thales.getModule();
        var Test = (function () {
            function Test(yahoo, pricing) {
                this.yahoo = yahoo;
                this.pricing = pricing;
                this.testBSAmericanApprox2002 = pricing.testBSAmericanApprox2002();
                this.testGBlackScholesImpVolBisection = pricing.testGBlackScholesImpVolBisection();
                this.testGBlackScholesNGreeks = pricing.testGBlackScholesNGreeks();
                this.testBSAmericanApprox2002ImpVolBisection = pricing.testBSAmericanApprox2002ImpVolBisection();
                this.testBSAmericanApprox2002NGreeks = pricing.testBSAmericanApprox2002NGreeks();
                yahoo.testDownload();
            }
            Test.$inject = [
                "Thales.Yahoo.Download",
                "Thales.Pricing.American"
            ];
            return Test;
        }());
        Controllers.Test = Test;
        app.controller("Thales.Controllers.Test", Test);
    })(Controllers = Thales.Controllers || (Thales.Controllers = {}));
})(Thales || (Thales = {}));
var Thales;
(function (Thales) {
    var app = Thales.getModule();
    app.config([
        "$routeProvider", function ($routeProvider) {
            $routeProvider.when("/test", {
                templateUrl: "views/test.html",
                controller: "Thales.Controllers.Test",
                controllerAs: "test"
            });
            $routeProvider.otherwise({
                templateUrl: "views/test.html",
                controller: "Thales.Controllers.Test",
                controllerAs: "test"
            });
        }
    ]);
})(Thales || (Thales = {}));
