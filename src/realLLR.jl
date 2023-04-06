"""
    realLLRcorr_col(Y,col)

Compute the log-likelihood ratios for Findr test 0 (**correlation test**) for a given column `col` (gene) of gene expression matrix `Y` against all other columns of `Y`.

`Y` is assumed to have undergone supernormalization with each column having mean zero and variance one. The LLRs are scaled by the number of rows (samples).
"""
function realLLRcorr_col(Y,col)
    ρ = vec(cov(Y,Y[:,col],corrected=false))
    ρ[col] = 1. # set self to exact value
    -0.5*log.(abs.(1 .- ρ.^2))
end

"""
    realLLRcausal_col(Y,E,col)

Compute for a given column `col` (gene) of gene expression matrix `Y` with categorical instrument `E` against all other columns of `Y` the log-likelihood ratios for Findr causal tests: 

- Test 2 (**Linkage test**) 
- Test 3 (**Mediation test**)
- Test 4 (**Relevance test**)
- Test 5 (**Pleiotropy test**)

`Y` is assumed to have undergone supernormalization with each column having mean zero and variance one. The LLRs are scaled by the number of rows (samples).
"""
function realLLRcausal_col(Y,E,col)
    # compute the sufficient statistics
    ρ, σ1, σ2 = llrstats_col(Y,E,col)

    # test 2
    llr2 = -0.5*log.(σ1)

    # test 4
    # the abs in the argument of the log is to have a positive argument when applied to the "col"th entry, where the exact value is 0
    llr4 = -0.5*log.(abs.(σ1[col].*σ1 .- (ρ .+ σ2 .- 1).^2)) .+ 0.5*log.(σ1[col])
    llr4[col] = Inf # set self to Inf

    # test 3
    llr3 = llr4 .+ 0.5*log.(abs.(1 .- ρ.^2))
    llr3[col] = Inf

    # test 5
    llr5 = llr4 .+ 0.5*log.(σ1)
    llr5[col] = Inf

    # output
    llr2, llr3, llr4, llr5
end

"""
    llrstats_col(Y,E,col)

Compute the sufficient statistics to compute the log-likelihood ratios for Findr tests 2-5  for a given column `col` (gene) of gene expression matrix `Y` with categorical instrument `E` against all other columns of `Y`.

`Y` is assumed to have undergone supernormalization with each col having mean zero and variance one. The LLRs are scaled by the number of rows (samples).

The sufficient statistics are:

- the covariance `ρ` between the given column `col` of matrix `Y` and all other columns of `Y`
- the weighted average variances `σ1` of each column of matrix `Y` over the groups (unique values) in `E`
- the weighted average covariance `σ2` between the given column `col` of `Y` and all other columns of `Y` over the groups (unique values) of `E`
"""
function llrstats_col(Y,E,col)
    ρ = vec(cov(Y,Y[:,col],corrected=false))
    ρ[col] = 1. # set self to exact values

    gs, μ = groupmeans(Y,E)
    
    ns = length(E) # number of samples
    w = pweights(gs/ns) # probability weights

    σ1 = vec(1 .- sum(μ.^2,w,dims=2))
    σ2 = vec(1 .- sum(μ.*μ[col,:]',w,dims=2))

    ρ, σ1, σ2
end

"""
    groupmeans(Y,E)

Compute the size and mean of each column of matrix `Y` for each of the groups (unique values) in categorical vector `E`.
"""
function groupmeans(Y,E)
    uE = unique(E)
    gs = zeros(Int,length(uE))
    μ = zeros(size(Y,2),length(uE))
    Threads.@threads for i = eachindex(uE)
        gs[i] = sum(E.==uE[i])
        μ[:,i] = mean(Y[E.==uE[i],:], dims=1)
    end
    gs, μ
end