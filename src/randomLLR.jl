struct LBeta{T<:Real} <: ContinuousUnivariateDistribution 
    α :: T
    β :: T
end

Distributions.params(d::LBeta) = (d.α, d.β)

function Distributions.pdf(d::LBeta, x::Real)
    α = params(d)[1]
    β = params(d)[2]
    return 2*(1-exp(-2x))^(0.5α-1)*exp(-β*x)/beta(0.5α,0.5β)
end