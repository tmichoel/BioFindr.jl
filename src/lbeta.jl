"""
    LBeta(α,β)
    
The *LBeta distribution* with parameters ``\\alpha`` and ``\\beta`` is defined as the distribution of a random variable 

``X=-\\frac{1}{2}(\\ln(1-Y))``

where ``Y\\sim\\operatorname{Beta}(\\alpha/2, \\beta/2)``
"""
# struct LBeta{T<:Real} <: ContinuousUnivariateDistribution 
#     α :: T
#     β :: T
# end

struct LBeta <: ContinuousUnivariateDistribution
    α :: Float64
    β :: Float64
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
    fit(LBeta, x)

Fit an `LBeta` distribution to data `x` using the method of moments by exploiting its relationship to the Beta distribution.

See https://github.com/JuliaStats/Distributions.jl/blob/master/src/univariate/continuous/beta.jl
"""
function Distributions.fit(::Type{<:LBeta}, x::AbstractArray{T}) where T<:Real
    z = 1 .-  exp.(-2 .* x) # if `x` is LBeta distributed, then `z` is Beta distributed
    bd = Distributions.fit(Beta,z)
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

"""
    fit_weighted(LBeta, x, w)

Fit an `LBeta` distribution to data `x` where each observation in `x` has a weight `w` between 0 and 1, using the method of moments. The function exploits the relationship to the Beta distribution and is a simple adaptation of the `fit` function in

https://github.com/JuliaStats/Distributions.jl/blob/master/src/univariate/continuous/beta.jl
"""
function fit_weighted(::Type{<:LBeta}, x::AbstractArray{T}, w::AbstractWeights) where T<:Real
    z = 1 .-  exp.(-2 .* x) # if `x` is LBeta distributed, then `z` is Beta distributed

    # the following lines adapated from beta.jl fit function by replacing calls to `mean` and `var` by their weighted versions
    z_bar = mean(z, w)
    v_bar = var(z, w; mean=z_bar, corrected=true)
    if v_bar < z_bar*(1-z_bar)
        temp = ((z_bar * (1. - z_bar)) / v_bar) - 1.
        α = z_bar * temp
        β = (1. - z_bar) * temp
    else
        α = 0.
        β = 0.
    end
    
    lbp = 2 .* [α,β] # multiply fitted Beta param to obtain LBeta parameters
    return LBeta(lbp...)
end

# """
#     beta_mixture_model(z)
# """
# @model function beta_mixture_model(z,α0,β0,π0)
#     if z === missing
#         z = Vector{Float64}(undef,1)
#     end

#     # Draw the parameters of the alternative distribution from an uninformative prior with constraints imposed by the null distribution
#     α ~ FlatPos(α0)
#     β ~ Uniform(0.,β0)

#     # Latent variable distribution
#     latent = Categorical([π0; 1-π0])

#     # Construct Beta distributions for the null and alternative
#     beta_components = [Beta(0.5α0, 0.5β0); Beta(0.5α, 0.5β)]

#     # Draw assignments for each observation
#     N = length(z)
#     u = Vector{Int}(undef, N)
#     for i = eachindex(z)
#         u[i] ~ latent
#         z[i] ~ beta_components[u[i]]
#     end
# end