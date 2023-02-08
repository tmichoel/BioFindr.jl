"""
    realLLRcorr_row(Y,row)

Compute the log-likelihood ratios for Findr test 0 (**correlation test**) for a given `row` (gene) of gene expression matrix `Y` against all other rows of `Y`.

`Y` is assumed to have undergone supernormalization with each row having mean zero and variance one. The LLRs are scaled by the number of samples.
"""
function realLLRcorr_row(Y,row)
    ρ = cov(Y[row,:],Y,dims=2,corrected=false)'
    ρ[row] = 1. # set self to exact values
    -0.5*log.(1 .- ρ.^2)
end

"""
    realLLRcausal_row(Y,E,row)

Compute for a given `row` (gene) of gene expression matrix `Y` with categorical instrument `E` against all other rows of `Y` the log-likelihood ratios for Findr causal tests: 

- Test 2 (**Linkage test**) 
- Test 3 (**Mediation test**)
- Test 4 (**Relevance test**)
- Test 5 (**Pleiotropy test**)

`Y` is assumed to have undergone supernormalization with each row having mean zero and variance one. The LLRs are scaled by the number of samples.
"""
function realLLRcausal_row(Y,E,row)
    # compute the sufficient statistics
    ρ, σ1, σ2 = llrstats_row(Y,E,row)

    # test 2
    llr2 = -0.5*log.(σ1)

    # test 4
    # the abs in the argument of the log is to have a positive argument when applied to the "row"th entry, where the exact value is 0
    llr4 = -0.5*log.(abs.(σ1[row].*σ1 .- (ρ .+ σ2 .- 1).^2)) .+ 0.5*log.(σ1[row])
    llr4[row] = Inf # set self to Inf

    # test 3
    llr3 = llr4 .+ 0.5*log.(1 .- ρ.^2)
    llr3[row] = Inf

    # test 5
    llr5 = llr4 .+ 0.5*log.(σ1)
    llr5[row] = Inf

    # output
    llr2, llr3, llr4, llr5
end

"""
    llrstats_row(Y,E,row)

Compute the sufficient statistics to compute the log-likelihood ratios for Findr tests 2-5  for a given `row` (gene) of gene expression matrix `Y` with categorical instrument `E` against all other rows of `Y`.

`Y` is assumed to have undergone supernormalization with each row having mean zero and variance one. The LLRs are scaled by the number of samples.

The sufficient statistics are:

- the covariance `ρ` between the given `row` of matrix `Y` and all other rows of `Y`
- the weighted average variances `σ1` of each row of matrix `Y` over the groups (unique values) in `E`
- the weighted average covariance `σ2` between the given `row` of `Y` and all other rows of `Y` over the groups (unique values) of `E`
"""
function llrstats_row(Y,E,row)
    ρ = cov(Y[row,:],Y,dims=2,corrected=false)'
    ρ[row] = 1. # set self to exact values

    gs, μ = groupmeans(Y,E)
    
    ns = length(E) # number of samples
    w = pweights(gs/ns) # probability weights

    σ1 = vec(1 .- sum(μ.^2,w,dims=2))
    σ2 = vec(1 .- sum(μ.*μ[row,:]',w,dims=2))

    ρ, σ1, σ2
end

"""
    groupmeans(Y,E)

Compute the size and mean of each row of matrix `Y` for each of the groups (unique values) in categorical vector `E`.
"""
function groupmeans(Y,E)
    uE = unique(E)
    gs = zeros(Int,length(uE))
    μ = zeros(size(Y)[1],length(uE))
    Threads.@threads for i = eachindex(uE)
        gs[i] = sum(E.==uE[i])
        μ[:,i] = mean(Y[:,E.==uE[i]], dims=2)
    end
    gs, μ
end