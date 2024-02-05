"""
    LBeta(α,β)
    
The *LBeta distribution* with parameters ``\\alpha`` and ``\\beta`` is defined as the distribution of a random variable 

``X=-\\frac{1}{2}(\\ln(1-Y))``

where ``Y\\sim\\operatorname{Beta}(\\alpha/2, \\beta/2)``
"""
struct LBeta <: ContinuousUnivariateDistribution
    α :: Float64
    β :: Float64
end

# struct LBeta{T<:Real} <: ContinuousUnivariateDistribution 
#     α :: T
#     β :: T
# end

"""
    params(d)

Get the parameters of an LBeta distribution.
"""
#Distributions.params(d::LBeta) = (d.α, d.β)
function params(d::LBeta)  
    (d.α, d.β)
end


"""
    pdf(d,x)

Evaluate the probability density function of an LBeta distribution with support on ``x\\geq 0``. 
"""
function Distributions.pdf(d::LBeta, x::Real)
    α,β = params(d)
    return x >= 0 ? 2*(1-exp(-2x))^(0.5α-1)*exp(-β*x)/beta(0.5α,0.5β) : 0.
end

"""
    logpdf(d,x)

Evaluate the probability density function of an LBeta distribution.  
"""
function Distributions.logpdf(d::LBeta, x::Real)
    α,β = params(d)
    return x >= 0 ? log(2) - logbeta(0.5α,0.5β) + (0.5α-1)*log(1-exp(-2x)) - β*x : 0.
end

"""
    cdf(d,x)

Evaluate the cumulative distribution function of an LBeta distribution using its relation to the Beta distribution. 
"""
function Distributions.cdf(d::LBeta, x::Real)
    bp = 0.5 .* params(d)
    bd = Beta(bp...) 
    return cdf(bd, 1-exp(-2x))
end

"""
    ccdf(d,x)

Evaluate the complementary cumulative distribution function of an LBeta distribution using its relation to the Beta distribution. 
"""
function Distributions.ccdf(d::LBeta, x::Real)
    bp = 0.5 .* params(d)
    bd = Beta(bp...)
    return ccdf(bd, 1-exp(-2x))
end

"""
    logccdf(d,x)

Evaluate the log complementary cumulative distribution function of an LBeta distribution using its relation to the Beta distribution. 
"""
function Distributions.logccdf(d::LBeta, x::Real)
    bp = 0.5 .* params(d)
    bd = Beta(bp...)
    return logccdf(bd, 1-exp(-2x))
end

"""
    compute_moments(d::LBeta)

Compute the first and second moments `m1` and `m2` of the Beta distribution corresponding to the given `LBeta` distribution. 
"""
function compute_moments(d::LBeta)
    α,β = params(d)
    m1 = α / (α + β)
    m2 = α * (α + 1) / ((α + β) * (α + β + 1))
    return m1, m2
end

"""
    test_moments(m1,m2)

Test if the given first and second moments `m1` and `m2` satisfy the conditions for a valid Beta distribution. (see also the [Beta distribution wiki](https://en.wikipedia.org/wiki/Beta_distribution#Two_unknown_parameters)):

    ``m_1>m_2 \\wedge m_2>m_1^2``
"""
function test_moments(m1,m2)
    return m1>m2 && m2>m1^2
end

"""
    fit_mom(LBeta, m1, m2)

Fit an `LBeta` distribution to given first and second moments `m1` and `m2` of the corresponding Beta distribution. This requires that `m1` and `m2` satisfy the following relations (see also the [Beta distribution wiki](https://en.wikipedia.org/wiki/Beta_distribution#Two_unknown_parameters)):

``m_1>m_2 \\wedge m_2>m_1^2``

An [AssertionError](https://docs.julialang.org/en/v1/base/base/#Core.AssertionError) is thrown if the condition evaluates to `false`.
"""
function fit_mom(LBeta, m1, m2)
    @assert test_moments(m1,m2) "Invalid Beta distribution moments."
    α = 2 * m1 * (m1 - m2) / (m2 - m1^2)
    β = 2 * (1 - m1) * (m1 - m2) / (m2 - m1^2)
    return LBeta(α,β)
end


"""
    fit(LBeta, x)

Fit an `LBeta` distribution to data `x` using the method of moments by exploiting its relationship to the Beta distribution.

See https://github.com/JuliaStats/Distributions.jl/blob/master/src/univariate/continuous/beta.jl
"""
function fit(::Type{<:LBeta}, x::AbstractArray{T}) where T<:Real
    z = 1 .-  exp.(-2 .* x) # if `x` is LBeta distributed, then `z` is Beta distributed
    bd = fit(Beta,z)
    lbp = 2 .* params(bd) # multiply fitted Beta param to obtain LBeta parameters
    return LBeta(lbp...)
end

"""
    fit_mle(LBeta, x)

Fit an `LBeta` distribution to data `x` using maximum-likelihood estimation by exploiting its relationship to the Beta distribution.

See https://github.com/JuliaStats/Distributions.jl/blob/master/src/univariate/continuous/beta.jl
"""
function Distributions.fit_mle(::Type{<:LBeta}, x::AbstractArray{T}) where T<:Real
    z = 1 .-  exp.(-2 .* x) # if `x` is LBeta distributed, then `z` is Beta distributed
    bd = Distributions.fit_mle(Beta,z)
    lbp = 2 .* params(bd) # multiply fitted Beta param to obtain LBeta parameters
    return LBeta(lbp...)
end