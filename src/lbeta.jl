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
    Sampling
"""
function sampler(d::LBeta)
    bp = 0.5 .* params(d)
    bd = Beta(bp...)
    return sampler(bd)
end

function rand(rng::AbstractRNG, d::LBeta)
    bp = 0.5 .* params(d)
    bd = Beta(bp...)
    x = rand(rng,bd)
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
#     fit_mle_weighted(LBeta, x, w)

# Fit an `LBeta` distribution to data `x` where each observation in `x` has a weight `w` between 0 and 1, using maximum likelihood. The function exploits the relationship to the Beta distribution and is a simple adaptation of the `fit_mle` function in

# https://github.com/JuliaStats/Distributions.jl/blob/master/src/univariate/continuous/beta.jl
# """
# function fit_mle_weighted(::Type{<:LBeta}, x::AbstractArray{T}, w::AbstractWeights; 
#     maxiter::Int=1000, tol::Float64=1e-14) where T<:Real

#     z = 1 .-  exp.(-2 .* x) # if `x` is LBeta distributed, then `z` is Beta distributed

#     # the following lines adapated from beta.jl fit function by replacing calls to `mean` and `var` by their weighted versions
    
#     α₀,β₀ = 0.5*params(fit_weighted(LBeta,x,w)) #initial guess of parameters
#     g₁ = mean(log.(z),w)
#     g₂ = mean(log.(one(T) .- z),w)
#     θ= [α₀ ; β₀ ]

#     converged = false
#     t=0
#     while !converged && t < maxiter #newton method
#         t+=1
#         temp1 = digamma(θ[1]+θ[2])
#         temp2 = trigamma(θ[1]+θ[2])
#         grad = [g₁+temp1-digamma(θ[1])
#                temp1+g₂-digamma(θ[2])]
#         hess = [temp2-trigamma(θ[1]) temp2
#                 temp2 temp2-trigamma(θ[2])]
#         Δθ = hess\grad #newton step
#         θ .-= Δθ
#         converged = dot(Δθ,Δθ) < 2*tol #stopping criterion
#     end


#     lbp = 2 .* θ # multiply fitted Beta param to obtain LBeta parameters
#     return LBeta(lbp...)
# end

