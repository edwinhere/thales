/// <reference path="model.d.ts"/>
/// <reference path="../main.ts"/>
module Thales.Yahoo {
    var app = getModule();
    export class Download {
        constructor(
            private $http: ng.IHttpService,
            private $q: ng.IQService,
            private $sce: ng.ISCEService
        ) {
        }
        public download(symbol: string): angular.IPromise<Pricing.Model.Underlying> {
            var result = this.downloadAll(symbol).then((promiseValue: Model.OptionChain[]) => {
                var price: number = promiseValue[0].result[0].quote.regularMarketPrice;
                var dividend: number = promiseValue[0].result[0].quote.trailingAnnualDividendYield;
                var expirations: Pricing.Model.Expiration[] = <Pricing.Model.Expiration[]>[];
                angular.forEach(promiseValue, (optionChain: Model.OptionChain, index: number) => {
                    optionChain.result[0].options;
                });
                return <Pricing.Model.Underlying>{
                    symbol: symbol,
                    price: price,
                    dividend: dividend,
                    expirations: expirations
                };
            });
            return result;
        }
        public downloadAll(symbol: string): angular.IPromise<Model.OptionChain[]> {
            var results = this.$http.jsonp(this.$sce.trustAsResourceUrl("https://query1.finance.yahoo.com/v7/finance/options/" + symbol))
                .then((promiseValue: angular.IHttpResponse<Model.RootObject>) => {
                    var rootObject: Model.RootObject = promiseValue.data;
                    if (rootObject.optionChain.error) {
                        rootObject.optionChain.error = "Error A ~ " + JSON.stringify(rootObject.optionChain.error);
                        return [rootObject.optionChain] as Model.OptionChain[];
                    } else {
                        var expirations: number[] = rootObject.optionChain.result[0].expirationDates;
                        var resultPromises: angular.IPromise<Model.OptionChain>[] = [] as angular.IPromise<Model.OptionChain>[];
                        angular.forEach(expirations, (expiration: number, key: number): void => {
                            resultPromises.push(this.downloadByExpiration(symbol, expiration));
                        });
                        var onePromise = this.$q.all(resultPromises);
                        return onePromise;
                    }
                }).catch((reason: any) => {
                    return [<Model.OptionChain>{
                        error: "Error C ~ " + JSON.stringify(reason),
                        result: undefined
                    }] as Model.OptionChain[];
                });
            return results;
        }

        private downloadByExpiration(symbol: string, expiration: number): angular.IPromise<Model.OptionChain> {
            var results: angular.IPromise<Model.OptionChain> = this.$http.jsonp(this.$sce.trustAsResourceUrl("https://query1.finance.yahoo.com/v7/finance/options/" + symbol + "?date=" + expiration))
                .then((promiseValue: angular.IHttpResponse<Model.RootObject>) => {
                    return promiseValue.data.optionChain;
                }).catch((reason: any) => {
                    return <Model.OptionChain>{
                        error: "Error D ~ " + JSON.stringify(reason),
                        result: undefined
                    };
                });
            return results;
        }

        public testDownload(): void {
            var promise = this.downloadAll("AAPL");
            promise.then((promiseValue: Model.OptionChain[]) => {
                var result = promiseValue.every((value: Model.OptionChain, index: number): boolean => {
                    if (value.error) {
                        console.error(value.error);
                    }
                    return !value.error;
                });
                this.testDownloadStatus = result;
            });
            this.testDownloadStatus = undefined;
        }

        public testDownloadStatus: boolean;
        public static $inject: string[] = ["$http", "$q", "$sce"]
    }
    app.service("Thales.Yahoo.Download", Download);
}