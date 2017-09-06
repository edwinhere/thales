/// <reference path="main.ts"/>
/// <reference path="pricing/american.ts"/>
/// <reference path="yahoo/download.ts"/>
/// <reference path="controllers/test.ts"/>
module Thales {
    var app = getModule();
    app.config([
        "$routeProvider", ($routeProvider: angular.route.IRouteProvider) => {
            $routeProvider.when("/test", {
                templateUrl: "views/test.html",
                controller: "Thales.Controllers.Test",
                controllerAs: "test"
            });
            $routeProvider.otherwise({
                // templateUrl: "views/main.html",
                templateUrl: "views/test.html",
                controller: "Thales.Controllers.Test",
                controllerAs: "test"
            });
        }
    ]);
}