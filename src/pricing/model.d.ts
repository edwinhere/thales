declare module Thales.Pricing.Model {
    export interface OptionChain {
        symbol: string,
        price: string
        dividend: number,
        calls: Strike[],
        puts : Strike[]
    }

    export interface Strike {
        [strike: number] : Option
    }

    export interface Option {
        expiryInYears: number,
        expiryInDays: number
        expiryAsDate: Date,
        volatility: number
    }
}