"""
    pi0est(pval)

Estimate the proportion π0 of truly null features in a vector `pval` of p-values using Storey's method

See also http://varianceexplained.org/files/pi0boot.pdf
"""
function pi0est(pval)
    λ = 0:0.05:0.95
    pval = sort(pval)
    m = length(pval)
    W = map(x -> sum(pval.>=x),λ)
    π0 = W ./ (m*(1 .- λ)) #map(x -> sum(pval.>=x)/(m*(1-x)),λ)
    minπ0 = minimum(π0) # quantile(π0,0.1) # 
    mse = (W ./ (m^2 * (1 .- λ).^2)) .* (1 .- W/m) .+ (π0 .- minπ0).^2
    π0[argmin(mse)]
end

"""
    fitdist(llr,ns,[ng,test])

Fit an LBeta distribution to a vector of log-likelihood ratio values `llr` for a given Findr `test` with sample size `ns` and number of genotype groups `ng`. The input variable `test` can take the values:

    - ':corr' - **correlation test** (test 0)
    - ':link' - **linkage test** (test 1/2)
    - ':med' - **mediation test** (test 3)
    - ':relev' - **relevance test** (test 4)
    - ':pleio' - **pleiotropy test** (test 5)
    
With two input arguments, the correlation test with `ns` samples is used. With three input arguments, or with four arguments and `test` equal to ":corr", the correlation test with `ns` samples is used and the third argument is ignored.
"""
function fitdist(llr,ns,ng=1,test=:corr)
    dnull = nulldist(ns,ng,test)
    α = Float64(params(dnull)[1])
    β0 = Float64(params(dnull)[2])
    z = log.(exp.(2 .* llr[llr.<Inf]) .- 1)
    β = min(β0 , 2 * invdigamma( digamma(0.5α) - mean(z) ) )
    dfit = LBeta(α,β)
    π0 = beta(0.5α, 0.5β0) / beta(0.5α, 0.5β)
    return dfit, π0
end

"""
    postprob(llr,ns,[ng,test])

Return posterior probabilities for a vector of log-likelihood ratio values `llr` for a given Findr `test` with sample size `ns` and number of genotype groups `ng`. The input variable `test` can take the values:

- ':corr' - **correlation test** (test 0)
- ':link' - **linkage test** (test 1/2)
- ':med' - **mediation test** (test 3)
- ':relev' - **relevance test** (test 4)
- ':pleio' - **pleiotropy test** (test 5)

With two input arguments, the correlation test with `ns` samples is used. With three input arguments, or with four arguments and `test` equal to ":corr", the correlation test with `ns` samples is used and the third argument is ignored.
"""
function postprob(llr,ns,ng=1,test=:corr)
    # Set the null distribution
    dnull = nulldist(ns, ng, test)
    # Evaluate the null distribution p.d.f. on the log-likelihood ratios
    pnull = pdf.(dnull, llr)
    # Fit an LBeta distribution with same α to the log-likelihood ratios
    dreal = fitdist(llr, ns, ng, test)
    # Evaluate the real distribution p.d.f. on the log-likelihood ratios
    preal = pdf.(dreal, llr)
    # Estimate the proportion of truly null features among the data
    π0 = beta(0.5.*params(dnull)...) / beta(0.5.*params(dreal)...) 
    # Compute the posterior probabilities. Take into account potential 
    1 .- π0 * pnull ./ preal
end