```@meta
CurrentModule = BioFindr
```


# General inference algorithm

Consider a set of observations sampled from a mixture distribution of a null and an alternative hypothesis. For instance in gene regulation, every observation can correspond to expression levels of a pair of genes wich are sampled from a bivariate normal distribution with zero (null hypothesis) or non-zero (alternative hypothesis) correlation coefficient. In BioFindr, we predict the probability that any sample follows the alternative hypothesis with the following algorithm (based on and modified from [^Chen2007]):

**1\.**  For robustness against outliers, we convert every continuous variable into standard normally distributed $N(0,1)$ values using a rank-based inverse normal transformation across all samples. We name this step as *supernormalization*.

```@docs
supernormalize
```

**2\.**  We propose a null and an alternative hypothesis for all [Likelihood ratio tests](@ref) of interest where, by definition, the null hypothesis space is a subset of the alternative hypothesis. Model parameters are replaced with their maximum likelihood estimators (MLEs) to obtain the log likelihood ratio (LLR) between the alternative and null hypotheses.

**3\.**  We derive the analytical expression for the probablity density function (PDF) of the [Null distributions of the log-likelihood ratios](@ref) when samples follow the null hypothesis.

**4\.**  We convert LLRs using [Bayesian inference of posterior probabilities](@ref) of the hypothesis of interest with empirical estimation of local false discovery rate.

**5\.** We consider multiple [Tests to evaluate](@ref), consisting of combinations of the basic [Likelihood ratio tests](@ref), for common tasks in genome-wide studies: 

- [Coexpression analysis](@ref) 
- [Differential expression analysis](@ref)
- [Causal inference](@ref)
- [Bipartite causal inference](@ref).



[^Chen2007]: Chen L, Emmert-Streib F, Storey J. [Harnessing naturally randomized transcription to infer regulatory relationships among genes](https://doi.org/10.1186/gb-2007-8-10-r219). Genome Biol 8, R219 (2007).
