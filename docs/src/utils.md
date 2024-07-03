```@meta
CurrentModule = BioFindr
```

# Utilities

## Data preparation

The [`findr`](@ref) function expects that input DataFrames use [scientific types](https://juliaai.github.io/ScientificTypes.jl/dev/), mainly to distinguish between count-based expression data and categorical genotype data which could otherwise both be presented as integer-valued data. A utility function [`coerce_scitypes!`](@ref) is provided to convert data to the right [scientific type](https://juliaai.github.io/ScientificTypes.jl/dev/):

```@docs
coerce_scitypes!
```

## Postprocessing functions

Several utility functions are used when [`findr`](@ref) is called with DataFrame inputs, some of which may be useful when manually post-processing the output of [`findr_matrix`](@ref) calls with matrix-based inputs.

```@docs
getpairs
symprobs
combineprobs
stackprobs
globalfdr!
globalfdr
qvalue
```

# Generating simulated data

BioFindr includes a function `generate_test_data` for generating simple simulated data for testing the package:

```@docs
generate_test_data
```