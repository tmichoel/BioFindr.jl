
"""
    nulldist(ns,[ng,test])

Return an LBeta distribution that is the null distribution of the log-likelihood ratio for a given Findr test with sample size `ns` and number of genotype groups `ng`. The input variable `test` can take the values:

- :corr - **correlation test** (test 0)
- :link - **linkage test** (test 1/2)
- :med - **mediation test** (test 3)
- :relev - **relevance test** (test 4)
- :pleio - **pleiotropy test** (test 5)

With only one input argument, the null distribution for the correlation test with `ns` samples is returned. With two input arguments, or  with three arguments and `test` equal to `:corr`, the null distribution for the correlation test with `ns` samples is returned and the second argument is ignored.
"""
function nulldist(ns,ng=1,test=:corr)
    if test==:corr
        return LBeta(1,ns-2)
    elseif test==:link
        return LBeta(ng-1,ns-ng)
    elseif test==:med
        return LBeta(ng-1,ns-ng-1)
    elseif test==:relev
        return LBeta(ng,ns-ng-1)
    elseif test==:pleio
        return LBeta(1,ns-ng-1)
    else
        error("the third argument must be a symbol from the set {:corr,:link,:med,:relev,:pleio}")
    end
end

"""
    nullpval(llr,ns,[ng,test])

Return p-values for a vector of log-likelihood ratio values `llr` under the null distribution of the log-likelihood ratio for a given Findr test with sample size `ns` and number of genotype groups `ng`. The input variable `test` can take the values:

- :corr - **correlation test** (test 0)
- :link - **linkage test** (test 1/2)
- :med - **mediation test** (test 3)
- :relev - **relevance test** (test 4)
- :pleio - **pleiotropy test** (test 5)

With two input arguments, the correlation test with `ns` samples is used. With three input arguments, or with four arguments and `test` equal to `:corr`, the correlation test with `ns` samples is used and the third argument is ignored.
"""
function nullpval(llr,ns,ng=1,test=:corr)
    # create null distribution
    nulld = nulldist(ns,ng,test)
    # check for Inf and NaN
    #tf = llr.==Inf | isnan(llr) 
    # compute p-values
    ccdf(nulld,llr)
end

"""
    nulllog10pval(llr,ns,[ng,test])

Return negative log10 p-values for a vector of log-likelihood ratio values `llr` under the null distribution of the log-likelihood ratio for a given Findr test with sample size `ns` and number of genotype groups `ng`. The input variable `test` can take the values:

- :corr - **correlation test** (test 0)
- :link - **linkage test** (test 1/2)
- :med - **mediation test** (test 3)
- :relev - **relevance test** (test 4)
- :pleio - **pleiotropy test** (test 5)

With two input arguments, the correlation test with `ns` samples is used. With three input arguments, or with four arguments and `test` equal to `:corr`, the correlation test with `ns` samples is used and the third argument is ignored.
"""
function nulllog10pval(llr,ns,ng=1,test=:corr)
    # create null distribution
    nulld = nulldist(ns,ng,test)
    # compute negative log10 p-values
    -logccdf(nulld,llr)/log(10) 
end

"""
    nullpdf(llr,ns,[ng,test])

Return probability distribution function evaluations for a vector of log-likelihood ratio values `llr` under the null distribution of the log-likelihood ratio for a given Findr test with sample size `ns` and number of genotype groups `ng`. The input variable `test` can take the values:

- :corr - **correlation test** (test 0)
- :link - **linkage test** (test 1/2)
- :med - **mediation test** (test 3)
- :relev - **relevance test** (test 4)
- :pleio - **pleiotropy test** (test 5)

With two input arguments, the correlation test with `ns` samples is used. With three input arguments, or with four arguments and `test` equal to `:corr`, the correlation test with `ns` samples is used and the third argument is ignored.
"""
function nullpdf(llr,ns,ng=1,test=:corr)
    # create null distribution
    nulld = nulldist(ns,ng,test)
    # check for Inf and NaN
    #tf = llr.==Inf | isnan(llr) 
    # compute p-values
    pdf.(nulld,llr)
end