```@meta
CurrentModule = Findr
```

# Findr.jl Documentation

This is the documentation for [Findr.jl](https://github.com/tmichoel/Findr.jl), an implementation of the [Findr software](https://github.com/lingfeiwang/findr) in [Julia](https://julialang.org/). 

The methods implemented in Findr were developed by [Lingfei Wang](https://github.com/lingfeiwang) and were first described in the paper ["Efficient and accurate causal inference with hidden confounders from genome-transcriptome variation data"](https://doi.org/10.1371/journal.pcbi.1005703).

This documentation is also an exploration of what an "executable methods paper" in computational biology could look like, by interweaving the mathematical description of the method (structured as and copied from the original paper) with the documentation of its implementation.

[Findr.jl](https://github.com/tmichoel/Findr.jl) is not a literal translation of the [original Findr software](https://github.com/lingfeiwang/findr). For those familiar with the original version, the following differences are to be paid attention to:

- The main [`findr(X,G,pairGX)`](@ref) interface function takes as input expression and genotype matrices, and a list of pairs to match eQTLs with a subset of genes. This avoids the need to reshape the gene expression data every time [`findr`](@ref) is called with a different set of eQTLs.

- Input and output are structured by **columns**, that is, in the gene expression and genotype data, columns are genes or SNPs and rows are samples, and in the posterior probability matrices, each column contains the probabilities of a causal relation from the gene corresponding to that column to all other genes. This is the opposite of the [original software](https://github.com/lingfeiwang/findr) where variables corresponded to rows. This is to boost performance as Julia stores arrays in [column-major](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-column-major) format.

- LLRs ... 

## Table of contents

```@contents
```