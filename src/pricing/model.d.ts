declare module Thales.Pricing.Model {
    export interface Underlying {
        symbol: string,
        price: number,
        dividend: number,
        expirations: Expiration[]
    }

    export interface Expiration {
        [expiryInDays: number] : OptionChain 
    }

    export interface OptionChain {
        calls: Strike[],
        puts : Strike[]
    }

    export interface Strike {
        [strike: number] : Option
    }

    export interface Option {
        expiryInYears: number,
        expiryAsDate: Date,
        volatility: number,
        theoPrice: number
    }
}