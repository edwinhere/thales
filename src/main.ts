module Thales {
    angular.module("thales", ["ngRoute"]);

    export var getModule: () => ng.IModule = () => {
        return angular.module("thales");
    }
}