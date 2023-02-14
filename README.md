# Findr.jl

Findr.jl is an implementation of the [Findr software](https://github.com/lingfeiwang/findr) in [Julia](https://julialang.org/). See the paper ["Efficient and accurate causal inference with hidden confounders from genome-transcriptome variation data"](https://julialang.org/) for  details of the methods implemented in Findr.

Findr.jl is not a literal translation of the original Findr code. Major differences are:

- The main Findr interface function takes as input a matrix of gene expression data, a matrix of eQTL genotype data, and a list of pairs to match eQTLs with a subset of genes. This avoids the need to reshape the gene expression data every time Findr is called with a different set of eQTLs.

- Input and output are structured by **columns**, that is, in the gene expression and genotype data, columns are genes or SNPs and rows are samples, and in the posterior probability matrix, each column contains the probabilities of a causal relation from the gene corresponding to that column to all other genes. This is the opposite of the [original software](https://github.com/lingfeiwang/findr). This is to boost performance as Julia stores arrays in [column-major](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-column-major) format.

- LLRs ... 

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tmichoel.github.io/Findr.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tmichoel.github.io/Findr.jl/dev/)
[![Build Status](https://github.com/tmichoel/Findr.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tmichoel/Findr.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/tmichoel/Findr.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tmichoel/Findr.jl)
