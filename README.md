# Findr.jl

Findr.jl is an implementation of the [Findr software](https://github.com/lingfeiwang/findr) in [Julia](https://julialang.org/). See the paper ["Efficient and accurate causal inference with hidden confounders from genome-transcriptome variation data"](https://julialang.org/) for  details of the methods implemented in Findr.

If you haven't used Findr before, see the [documentation](https://tmichoel.github.io/Findr.jl).

If you have used Findr before, be aware that Findr.jl is not a literal translation of the original code. Major differences are:

- The main Findr interface function takes as input a matrix of gene expression data, a matrix of eQTL genotype data, and a list of pairs to match eQTLs with a subset of genes. This avoids the need to reshape the gene expression data every time Findr is called with a different set of eQTLs.

- Input and output are structured by **columns**, that is, in the input data, columns must correspond to variables (genes or eQTLs) and rows to samples, and in the posterior probability matrices, each column contains the probabilities of a causal relation from the gene corresponding to that column to all other genes. This is the *opposite* of the [original software](https://github.com/lingfeiwang/findr). This is to boost performance as Julia stores arrays in [column-major](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-column-major) format, and to make it easier for the software to interact with [DataFrames](https://dataframes.juliadata.org/), where columns tend to correspond to variables.

- Posterior probabilities are estimated by fitting a 2-parameter mixture distribution to the log-likelihood ratio (LLR) values. This avoids the need for the complicated adaptive histogram fitting in the [original software](https://github.com/lingfeiwang/findr), and allows to compute a goodness-of-fit in the form of testing the uniformity of p-values for observing the real LLRs under the fitted real distribution. See the [documentation](https://tmichoel.github.io/Findr.jl) for details.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tmichoel.github.io/Findr.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tmichoel.github.io/Findr.jl/dev/)
[![Build Status](https://github.com/tmichoel/Findr.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tmichoel/Findr.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/tmichoel/Findr.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tmichoel/Findr.jl)
