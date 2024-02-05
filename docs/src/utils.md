```@meta
CurrentModule = BioFindr
```

# Utilities

Utility functions which are used when [`findr`](@ref) is called with DataFrame inputs, some of which may be useful when manually post-processing the output of [`findr`](@ref) calls with matrix-based inputs.

```@docs
getpairs
symprobs
combineprobs
stackprobs
globalfdr!
globalfdr
qvalue
```