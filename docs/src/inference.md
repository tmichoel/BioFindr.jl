```@meta
CurrentModule = Findr
```


# General inference algorithm

Consider a set of observations sampled from a mixture distribution of a null and an alternative hypothesis. For instance in gene regulation, every observation can correspond to expression levels of a pair of genes wich are sampled from a bivariate normal distribution with zero (null hypothesis) or non-zero (alternative hypothesis) correlation coefficient. In Findr, we predict the probability that any sample follows the alternative hypothesis with the following algorithm (based on and modified from [@Chen:2007]):

1.  For robustness against outliers, we convert every continuous variable into standard normally distributed $N(0,1)$ values using a rank-based inverse normal transformation across all samples. We name this step as *supernormalization*.
2.  We propose a null and an alternative hypothesis for every likelihood ratio test (LRT) of interest where, by definition, the null hypothesis space is a subset of the alternative hypothesis. Model parameters are replaced with their maximum likelihood estimators (MLEs) to obtain the log likelihood ratio (LLR) between the alternative and null hypotheses.
3.  We derive the analytical expression for the probablity density function (PDF) of the LLR when samples follow the null hypothesis.
4.  We convert LLRs into posterior probabilities of the hypothesis of interest with the empirical estimation of local FDR.

```@docs
findr
```