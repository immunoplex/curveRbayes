# Compile a curveRbayes Stan Model (Cached)

Compiles the Stan model via cmdstanr. Compilation is cached by cmdstanr
so subsequent calls are instant.

## Usage

``` r
compile_stan_model(model_family = "logistic4")
```

## Arguments

- model_family:

  Character. One of the curveRcore model names.

## Value

A `CmdStanModel` object.
