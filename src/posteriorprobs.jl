"""
    pprob_corr_row(Y,row)

Compute the posterior probabilities for Findr test 0 (**correlation test**) for a given `row` (gene) of gene expression matrix `Y` against all other rows of `Y`.

`Y` is assumed to have undergone supernormalization with each row having mean zero and variance one. The LLRs are scaled by the number of samples.
"""
function pprob_corr_row(Y,row)
    # number of samples
    ns = size(Y,2) 
    # log-likelihood ratios
    llr = realLLRcorr_row(Y,row) 
    # remove diagonal element?
    
    # posterior probabilities
    pp, dalt, dnull = fit_mixdist_EM(llr,ns)
    return pp, π0, dalt, dnull
end

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
    π0 = W ./ (m*(1 .- λ))
    minπ0 = minimum(π0) # quantile(π0,0.1) # 
    mse = (W ./ (m^2 * (1 .- λ).^2)) .* (1 .- W/m) .+ (π0 .- minπ0).^2
    π0[argmin(mse)]
end

"""
    fit_mixdist_EM(llr,pi0,ns,ng=1,test=:corr; maxiter=1000, tol=1e-14)

Fit a two-component mixture distribution of two LBeta distribution to a vector of log-likelihood ratios `llr` using an EM algorithm. The first component is the true null distribution for a given Findr `test` with sample size `ns` and number of genotype groups `ng`. The second component is the alternative distribution, assumed to follow an LBeta distribution. The prior probability `pi0` of an observation belonging to the null component is fixed and determined by the `pi0est` function. Hence only the parameters of the alternative component need to be estimated.

The EM algorithm outputs posterior probabilities of the alternative hypothesis being true, in the form of the estimated recognition distribution. The optional parameters `maxiter` (default value 1000) and `tol` (default value 1e-14) control the convergence of the EM algorithm.

The input variable `test` can take the values:

- ':corr' - **correlation test** (test 0)
- ':link' - **linkage test** (test 1/2)
- ':med' - **mediation test** (test 3)
- ':relev' - **relevance test** (test 4)
- ':pleio' - **pleiotropy test** (test 5)

With two input arguments, the correlation test with `ns` samples is used. With three input arguments, or with four arguments and `test` equal to ":corr", the correlation test with `ns` samples is used and the third argument is ignored.
"""
function fit_mixdist_EM(llr,ns,ng=1,test=:corr; maxiter::Int=1000, tol::Float64=1e-5)

    dnull = nulldist(ns, ng, test) # null distribution
    pnull = pdf.(dnull, llr) # pvalues under the null hypothesis
    π0 = pi0est(pnull) # estimated proportion of true nulls

    # Initial guess for α and β: take the parameters of the null distribution and multiply α0 by a factor greater than one to ensure that in the limit LLR -> 0, all observations come from the null
    α, β = params(dnull)
    α = α/π0

    # Set the current alternative distribution
    dalt = LBeta(α,β)

    # Set the current recognition / posterior probability values
    palt = pdf.(dalt, llr) # pvalues under the alternative distribution
    pp = (1-π0) .*  palt./ (π0 .* pnull .+ (1-π0) .* palt) # posterior probabilities

    # EM until convergence
    converged = false
    it = 0
    while !converged && it < maxiter
        it += 1
        # Update α, β using current pp
        w = pweights(pp)
        dalt = fit_weighted(LBeta, llr, w)
        if dalt.α < dnull.α
            dalt = LBeta(dnull.α,dalt.β)
        end
        # Update pp using new params
        palt = pdf.(dalt, llr)
       # π0 = 1 - sum(pp)/length(pp)
        pp = (1-π0) .*  palt./ (π0 .* pnull .+ (1-π0) .* palt)
        
        # Check convergence
        converged = norm(pp.-w,Inf) < tol 
        #println(norm(pp.-w,Inf))
    end
    println(it)
    return pp, dalt, dnull
end

"""
    fit_mixdist_KDE(llr,ns,[ng,test])

Return posterior probabilities for a vector of log-likelihood ratio values `llr` for a given Findr `test` with sample size `ns` and number of genotype groups `ng`. The input variable `test` can take the values:

- ':corr' - **correlation test** (test 0)
- ':link' - **linkage test** (test 1/2)
- ':med' - **mediation test** (test 3)
- ':relev' - **relevance test** (test 4)
- ':pleio' - **pleiotropy test** (test 5)

With two input arguments, the correlation test with `ns` samples is used. With three input arguments, or with four arguments and `test` equal to ":corr", the correlation test with `ns` samples is used and the third argument is ignored.
"""
function fit_mixdist_KDE(llr,π0,ns,ng=1,test=:corr)
    # Set the null distribution
    dnull = nulldist(ns, ng, test)
    
    # Evaluate the null distribution p.d.f. on the log-likelihood ratios
    pnull = pdf.(dnull, llr)
    
    # Fit and evaluate the real distribution p.d.f. on the log-likelihood ratios
    preal = fit_kde(llr)

    # Compute the posterior probabilities of the alternative hypothesis being true
    pp = 1 .- π0 * pnull ./ preal

    # Smoothen posterior probs and make monotonically increasing

    # Fit an LBeta distribution with same α to the log-likelihood ratios
    # dreal,π0 = fitdist(llr, ns, ng, test)
    # # Evaluate the real distribution p.d.f. on the log-likelihood ratios
    # preal = pdf.(dreal, llr)
    # # Estimate the proportion of truly null features among the data
    # π0 = beta(0.5.*params(dnull)...) / beta(0.5.*params(dreal)...) 
    # Compute the posterior probabilities. 
    # 1 .- exp.(-(dnull.β - dreal.β)*llr)
end


"""
    fit_kde(llr)

Fit a distribution function to a vector of log-likelihood ratio values `llr` using kernel density estimation. To avoid boundary effects in the KDE, the log-likelihoods (which take values in ``[0,\\infty)``) are first transformed to a vector of `z`,
    
``
z = \\log \\left( e^{2 LLR} - 1 \\right)
``

which takes values in ``(-\\infty,\\infty)``. KDE is applied to `z`, and the probability density function (pdf) for `llr` is obtained from the pdf for `z` by the usual transformation rule for functions of random variables.
"""
function fit_kde(llr)
    # Transform llr values to avoid edge effects in density estimation
    z = log.(exp.(2 .* llr[llr.<Inf]) .- 1)
    
    # Apply kernel density estimation to transformed values
    dfit = kde(z)

    # Evaluate pdf at llr values using rule for transformations of random variables
    pd = zeros(size(llr))
    pd[llr.<Inf] = 2 * pdf(dfit, z) .* (1 .+ exp.(-z))
    return pd
end
