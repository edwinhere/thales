module Thales.Controllers {
    var app = getModule();
    export class Test {
        private testBSAmericanApprox2002: boolean;
        private testGBlackScholesImpVolBisection: boolean;
        private testGBlackScholesNGreeks: boolean;
        private testBSAmericanApprox2002ImpVolBisection: boolean;
        private testBSAmericanApprox2002NGreeks: boolean;
        constructor(
            private yahoo: Thales.Yahoo.Download,
            private pricing: Thales.Pricing.American
        ) {
            this.testBSAmericanApprox2002 = pricing.testBSAmericanApprox2002();
            this.testGBlackScholesImpVolBisection = pricing.testGBlackScholesImpVolBisection();
            this.testGBlackScholesNGreeks = pricing.testGBlackScholesNGreeks();
            this.testBSAmericanApprox2002ImpVolBisection = pricing.testBSAmericanApprox2002ImpVolBisection();
            this.testBSAmericanApprox2002NGreeks = pricing.testBSAmericanApprox2002NGreeks();
            yahoo.testDownload();
        }
        public static $inject: string[] = [
            "Thales.Yahoo.Download",
            "Thales.Pricing.American"
        ];
    }
    app.controller("Thales.Controllers.Test", Test);
}