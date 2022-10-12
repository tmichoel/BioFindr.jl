"""
    LBeta(α,β)
    
The *LBeta distribution* with parameters ``\\alpha`` and ``\\beta`` is defined as the distribution of a random variable 

``X=-\\frac{1}{2}(\\ln(1-Y))``

where ``Y\\sim\\operatorname{Beta}(\\alpha, \\beta)``
"""
struct LBeta{T<:Real} <: ContinuousUnivariateDistribution 
    α :: T
    β :: T
end

"""
    params(d)

Get the parameters of an LBeta distribution.
"""
Distributions.params(d::LBeta) = (d.α, d.β)

"""
    pdf(d,x)

Evaluate the probability density function of an LBeta distribution with support on ``x\\geq 0``. 
"""
function Distributions.pdf(d::LBeta, x::Real)
    α = params(d)[1]
    β = params(d)[2]
    return x >= 0 ? 2*(1-exp(-2x))^(0.5α-1)*exp(-β*x)/beta(0.5α,0.5β) : 0.
end

"""
    cdf(d,x)

Evaluate the cumulative distribution function of an LBeta distribution using its relation to the Beta distribution. 
"""
function Distributions.cdf(d::LBeta, x::Real)
    bd = Beta(params(d)...) 
    return cdf(bd, 1-exp(-2x))
end

"""
    randomLLRcorr(ns)

Return an LBeta distributed random variable for the null distribution of the log-likelihood ratio for Findr test 0 (**correlation test**) with sample size `ns`. 
"""
function randomLLRcorr(ns)
    return LBeta(1,ns-2)
end

"""
    randomLLRlink(ns,ng)

Return an LBeta distributed random variable for the null distribution of the log-likelihood ratio for Findr test 1/2 (**linkage test**) with sample size `ns` and number of genotype groups `ng`. 
"""
function randomLLRlink(ns,ng)
    return LBeta(ng-1,ns-ng)
end

"""
    randomLLRmed(ns,ng)

Return an LBeta distributed random variable for the null distribution of the log-likelihood ratio for Findr test 3 (**mediation test**) with sample size `ns` and number of genotype groups `ng`. 
"""
function randomLLRmed(ns,ng)
    return LBeta(ng-1,ns-ng-1)
end

"""
    randomLLRRrelev(ns,ng)

Return an LBeta distributed random variable for the null distribution of the log-likelihood ratio for Findr test 4 (**relevance test**) with sample size `ns` and number of genotype groups `ng`. 
"""
function randomLLRRrelev(ns,ng)
    return LBeta(ng,ns-ng-1)
end

"""
    randomLLRpleio(ns,ng)

Return an LBeta distributed random variable for the null distribution of the log-likelihood ratio for Findr test 5 (**pleiotropy test**) with sample size `ns` and number of genotype groups `ng`. 
"""
function randomLLRpleio(ns,ng)
    return LBeta(1,ns-ng-1)
end