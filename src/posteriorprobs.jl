"""
    pprob_corr_col(Y,col)

Compute the posterior probabilities for Findr test 0 (**correlation test**) for a given column `col` (gene) of gene expression matrix `Y` against all other columns of `Y`.

`Y` is assumed to have undergone supernormalization with each column having mean zero and variance one. The LLRs are scaled by the number of rows (samples).
"""
function pprob_corr_col(Y,col)
    # number of samples
    ns = size(Y,1) 
    # log-likelihood ratios
    llr = realLLRcorr_col(Y,col) 
    # vector to hold the output
    pp = ones(size(llr))
    # posterior probabilities, diagonal element set to 1 by default
    pp[Not(col)], dreal = fit_mixdist_EM(llr[Not(col)],ns)
    #pp[Not(col)] = fit_mixdist_KDE(llr[Not(col)],ns)
    return pp
end

"""
    pprob_causal_col(Y,col)

Compute the posterior probabilities for a given column `col` (gene) of gene expression matrix `Y` with categorical instrument `E` against all other columns of `Y` for Findr causal tests: 

    - Test 2 (**Linkage test**) 
    - Test 3 (**Mediation test**)
    - Test 4 (**Relevance test**)
    - Test 5 (**Pleiotropy test**)
    
`Y` is assumed to have undergone supernormalization with each column having mean zero and variance one.

For test 2, 4, and 5 the posterior probabilities are the probabilities of the alternative hypothesis being true. For test 3 they are the probabilities of the null hypothesis being true.
"""
function pprob_causal_col(Y,E,col)
    # number of samples and groups
    ns = size(Y,1) 
    ng = length(unique(E))
    # log-likelihood ratios
    llr2, llr3, llr4, llr5 = realLLRcausal_col(Y,E,col)
    # allocate output array
    pp = ones(length(llr2),4)
    # posterior probabilities for test 2, here we keep the true diagonal element
    pp[:,1], _ = fit_mixdist_EM(llr2,ns,ng,:link)
    # posterior probabilities for test 3, here we keep the default diagonal element and swap the role of null and alternative
    pp[Not(col),2], _ = fit_mixdist_EM(llr3[Not(col)],ns,ng,:med)
    # posterior probabilities for test 4 and 5, here we keep the default diagonal element
    pp[Not(col),3], _ = fit_mixdist_EM(llr4[Not(col)],ns,ng,:relev)
    pp[Not(col),4], _ = fit_mixdist_EM(llr5[Not(col)],ns,ng,:pleio)

    return pp
end



"""
    fit_mixdist_EM(llr,ns,ng=1,test=:corr; maxiter=1000, tol=1e-14)

Fit a two-component mixture distribution of two LBeta distributions to a vector of log-likelihood ratios `llr` using an EM algorithm. The first component is the true null distribution for a given Findr `test` with sample size `ns` and number of genotype groups `ng`. The second component is the alternative distribution, assumed to follow an LBeta distribution. The prior probability `pi0` of an observation belonging to the null component is fixed and determined by the `pi0est` function. Hence only the parameters of the alternative component need to be estimated.

The EM algorithm outputs posterior probabilities of the alternative hypothesis being true, in the form of the estimated recognition distribution. The optional parameters `maxiter` (default value 1000) and `tol` (default value 1e-3) control the convergence of the EM algorithm.

The input variable `test` can take the values:

- ':corr' - **correlation test** (test 0)
- ':link' - **linkage test** (test 1/2)
- ':med' - **mediation test** (test 3)
- ':relev' - **relevance test** (test 4)
- ':pleio' - **pleiotropy test** (test 5)

With two input arguments, the correlation test with `ns` samples is used. With three input arguments, or with four arguments and `test` equal to ":corr", the correlation test with `ns` samples is used and the third argument is ignored.
"""
function fit_mixdist_EM(llr,ns,ng=1,test=:corr; maxiter::Int=1000, tol::Float64=1e-3)

    # set null distribution and estimate proportion of true nulls
    dnull = nulldist(ns, ng, test) 
    π0 = pi0est( nullpval(llr, ns, ng, test) )

    # Initial guess for the alternative distribution: fit an LBeta distribution to the 50% largest llr values
    dalt = fit(LBeta,llr[llr .> quantile(llr,0.5)])
#    println([dalt.α,dalt.β])

    # Set the current recognition / posterior probability values
    pnull = pdf.(dnull, llr)
    palt = pdf.(dalt, llr) 
    pp = (1-π0) .*  palt./ (π0 .* pnull .+ (1-π0) .* palt)

    # EM until convergence
    converged = false
    it = 0
    while !converged && it < maxiter
        it += 1
        # Update the alternative distribution using current pp
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
    # println(it)
    # Set mixture distribution
    dreal = MixtureModel(LBeta[dnull, dalt],[π0, 1-π0])

    # Return posterior probabilities and estimated mixture distribution
    return pp, dreal
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
function fit_mixdist_KDE(llr,ns,ng=1,test=:corr)
    # Set the null distribution
    dnull = nulldist(ns, ng, test)
    π0 = pi0est( nullpval(llr, ns, ng, test) )
    
    # Evaluate the null distribution p.d.f. on the log-likelihood ratios
    pnull = pdf.(dnull, llr)
    
    # Fit and evaluate the real distribution p.d.f. on the log-likelihood ratios
    preal = fit_kde(llr)

    # Compute the posterior probabilities of the alternative hypothesis being true
    pp = 1 .- π0 * pnull ./ preal

    # Smoothen posterior probs and make monotonically increasing

    return pp
    
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
    lfdr(pval)

Estimate local false discovery rates for a vector `pval` of p-values using Storey's method.
"""
function lfdr(pval)
    
end