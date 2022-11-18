"""
    pi0est(pval)

Estimate the proportion π0 of truly null features in a vector `pval` of p-values using Storey's method
"""
function pi0est(pval)
    λ = 0:0.05:0.95
    pval = sort(pval)
    m = length(pval)
    π0 = map(x -> sum(pval.>=x)/(m*(1-x)),λ)
    minπ0 = quantile(π0,0.1)
    W = map(x -> sum(pval.>x),λ)
    mse = (W ./ (m^2 * (1 .- λ).^2)) .* (1 .- W/m) .+ (π0 .- minπ0).^2
    π0[argmin(mse)]
end

"""
    pi0est2(pval)

Estimate the proportion π0 of truly null features in a vector `llr` of log-likelihood ratios
"""
function pi0est2(llr)
    λ = 0:0.05:1.
    llr = sort(llr)
    m = length(llr)
    π0 = zeros(size(λ))
    for k=eachindex(λ)
        π0[k] = sum(llr.<λ[k])/(m*λ[k])
    end
    π0
end
