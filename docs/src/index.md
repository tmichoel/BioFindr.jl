```@meta
CurrentModule = Findr
```

# Findr.jl Documentation

This is the documentation for [Findr.jl](https://github.com/tmichoel/Findr.jl), an implementation of the [Findr software](https://github.com/lingfeiwang/findr) in [Julia](https://julialang.org/). 

The methods implemented in Findr were developed by [Lingfei Wang](https://github.com/lingfeiwang) and [Tom Michoel](https://github.com/tmichoel), and were first described in the paper ["Efficient and accurate causal inference with hidden confounders from genome-transcriptome variation data"](https://doi.org/10.1371/journal.pcbi.1005703).[^Wang2017]

This documentation copies both the structure and (most of) the contents of the [Materials and methods section](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005703#sec009) of the original paper, directly linking the mathematical description of the method to the documentation and source code of its implementation. Maybe this is how methods papers in computational biology should be written in the first place?

If you haven't used Findr before, the [FindrTutorials website](https://tmichoel.github.io/FindrTutorials) may be useful **(under development, link not active yet)**.

If you are familiar with the [original Findr software](https://github.com/lingfeiwang/findr), you should be aware that
[Findr.jl](https://github.com/tmichoel/Findr.jl) is not a literal translation. In particular, pay attention to the following differences:

- The main `findr` interface function for [Causal inference](@ref) takes as input expression and genotype matrices or [DataFrames](https://dataframes.juliadata.org/stable/), and a list or [DataFrame](https://dataframes.juliadata.org/stable/) of pairs to match eQTLs with a subset of genes. This avoids the need to reshape the gene expression data every time `findr` is called with a different set of eQTLs.

- Input and output are structured by **columns**, that is, in the gene expression and genotype data, columns are genes or SNPs and rows are samples, and in the posterior probability matrices, each column contains the probabilities of a causal relation from the gene corresponding to that column to all other genes. This is the **opposite** of the [original software](https://github.com/lingfeiwang/findr) where variables corresponded to rows. This is to boost performance as Julia stores arrays in [column-major](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-column-major) format.

-  If you call `findr` with [DataFrame](https://dataframes.juliadata.org/stable/) inputs (where columns naturally correspond to variables, that is, genes or SNPs), you need not worry about remembering the role of rows and columns in the output, as the output is returned in the form of a [DataFrame](https://dataframes.juliadata.org/stable/) with `Source`, `Target`, `Probability`, and `q-value` columns.

- You can pass a desired global FDR value for filtering inferred associations as a parameter when calling `findr`, no more need for manual post-processing of the output.  

- Estimation of the observed distribution of log-likelihood ratios uses either a [new, parametric method of moments](@ref mom_postprobs), or a [kernel density estimation method](@ref kde_postprobs), replacing the previous histogram-based method. 

## Table of contents

```@contents
```

[^Wang2017]: Wang L, Michoel T (2017) [Efficient and accurate causal inference with hidden confounders from genome-transcriptome variation data](https://doi.org/10.1371/journal.pcbi.1005703). PLoS Comput Biol 13(8): e1005703.
