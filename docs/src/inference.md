```@meta
CurrentModule = Findr
```


# General inference algorithm

Consider a set of observations sampled from a mixture distribution of a null and an alternative hypothesis. For instance in gene regulation, every observation can correspond to expression levels of a pair of genes wich are sampled from a bivariate normal distribution with zero (null hypothesis) or non-zero (alternative hypothesis) correlation coefficient. In Findr, we predict the probability that any sample follows the alternative hypothesis with the following algorithm (based on and modified from [^Chen2007]):

1.  For robustness against outliers, we convert every continuous variable into standard normally distributed $N(0,1)$ values using a rank-based inverse normal transformation across all samples. We name this step as *supernormalization*.
2.  We propose a null and an alternative hypothesis for every likelihood ratio test (LRT) of interest where, by definition, the null hypothesis space is a subset of the alternative hypothesis. Model parameters are replaced with their maximum likelihood estimators (MLEs) to obtain the log likelihood ratio (LLR) between the alternative and null hypotheses.
3.  We derive the analytical expression for the probablity density function (PDF) of the LLR when samples follow the null hypothesis.
4.  We convert LLRs into posterior probabilities of the hypothesis of interest with the empirical estimation of local FDR.

While primarily developed for causal inference, several common tasks in genome-wide studies can be implemented using this algorithm and are supported in Findr. 

[^Chen2007]: Chen L, Emmert-Streib F, Storey J. [Harnessing naturally randomized transcription to infer regulatory relationships among genes](https://doi.org/10.1186/gb-2007-8-10-r219). Genome Biol 8, R219 (2007).

## Coexpression analysis

```@docs
findr(X::Matrix{T}) where T<:AbstractFloat
findr(dX::T) where T<:AbstractDataFrame
```

## Differential expression analysis

```@docs
findr(X::Matrix{T}, G::Array{S}) where {T<:AbstractFloat, S<:Integer}
findr(dX::T, dG::T) where T<:AbstractDataFrame
```

## Causal inference

```@docs
findr(X::Matrix{T},G::Matrix{S},pairGX::Matrix{R}; combination="none") where {T<:AbstractFloat, S<:Integer, R<:Integer}
findr(dX::T, dG::T, dE::T; colG=1, colX=2, combination="IV") where T<:AbstractDataFrame
```

## Bipartite causal inference

```@docs
findr(X1::Matrix{T}, X2::Array{T}, G::Array{S}, pairGX::Matrix{R}; combination="none")  where {T<:AbstractFloat, S<:Integer, R<:Integer}
findr(dX1::T, dX2::T, dG::T, dE::T; colG=1, colX=2, combination="IV") where T<:AbstractDataFrame
```