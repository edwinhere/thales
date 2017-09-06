declare module Thales.Yahoo.Model {

    export interface Quote {
        language: string;
        quoteType: string;
        quoteSourceName: string;
        currency: string;
        twoHundredDayAverage: number;
        twoHundredDayAverageChange: number;
        tradeable: boolean;
        exchange: string;
        regularMarketPrice: number;
        regularMarketTime: number;
        regularMarketChange: number;
        regularMarketOpen: number;
        regularMarketDayHigh: number;
        regularMarketDayLow: number;
        regularMarketVolume: number;
        priceHint: number;
        earningsTimestampEnd: number;
        trailingAnnualDividendRate: number;
        trailingPE: number;
        trailingAnnualDividendYield: number;
        epsTrailingTwelveMonths: number;
        epsForward: number;
        sharesOutstanding: number;
        bookValue: number;
        fiftyDayAverage: number;
        fiftyDayAverageChange: number;
        twoHundredDayAverageChangePercent: number;
        forwardPE: number;
        priceToBook: number;
        sourceInterval: number;
        exchangeTimezoneName: string;
        exchangeTimezoneShortName: string;
        gmtOffSetMilliseconds: number;
        longName: string;
        market: string;
        shortName: string;
        postMarketPrice: number;
        postMarketChange: number;
        regularMarketChangePercent: number;
        regularMarketPreviousClose: number;
        bid: number;
        ask: number;
        bidSize: number;
        askSize: number;
        marketState: string;
        fiftyDayAverageChangePercent: number;
        marketCap: number;
        postMarketChangePercent: number;
        postMarketTime: number;
        messageBoardId: string;
        fullExchangeName: string;
        financialCurrency: string;
        averageDailyVolume3Month: number;
        averageDailyVolume10Day: number;
        fiftyTwoWeekLowChange: number;
        fiftyTwoWeekLowChangePercent: number;
        fiftyTwoWeekHighChange: number;
        fiftyTwoWeekHighChangePercent: number;
        fiftyTwoWeekLow: number;
        fiftyTwoWeekHigh: number;
        dividendDate: number;
        earningsTimestamp: number;
        earningsTimestampStart: number;
        symbol: string;
    }

    export interface Call {
        contractSymbol: string;
        strike: number;
        currency: string;
        lastPrice: number;
        change: number;
        percentChange: number;
        volume: number;
        openInterest: number;
        bid: number;
        ask: number;
        contractSize: string;
        expiration: number;
        lastTradeDate: number;
        impliedVolatility: number;
        inTheMoney: boolean;
    }

    export interface Put {
        contractSymbol: string;
        strike: number;
        currency: string;
        lastPrice: number;
        change: number;
        percentChange: number;
        volume: number;
        openInterest: number;
        bid: number;
        ask: number;
        contractSize: string;
        expiration: number;
        lastTradeDate: number;
        impliedVolatility: number;
        inTheMoney: boolean;
    }

    export interface Option {
        expirationDate: number;
        hasMiniOptions: boolean;
        calls: Call[];
        puts: Put[];
    }

    export interface Result {
        underlyingSymbol: string;
        expirationDates: number[];
        strikes: number[];
        hasMiniOptions: boolean;
        quote: Quote;
        options: Option[];
    }

    export interface OptionChain {
        result: Result[];
        error?: any;
    }

    export interface RootObject {
        optionChain: OptionChain;
    }

}

