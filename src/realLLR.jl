"""
    realLLRcorr_row(Y,row)

Compute the log-likelihood ratios for Findr test 0 (**correlation test**) for a given `row` (gene) of gene expression matrix `Y` against all other rows of `Y`.

`Y` is assumed to have undergone supernormalization with each row having mean zero and variance one. The LLRs are scaled by the number of samples.
"""
function realLLRcorr_row(Y,row)
    ρ = cov_row(Y,row)
    return -0.5*log.(1 .- ρ.^2)
end

"""
    realLLRlink_row(Y,E)

Compute the log-likelihood ratios for Findr test 2 (**linkage test**) for categorical vector `E` against all rows of matrix `Y`.

`Y` is assumed to have undergone supernormalization with each row having mean zero and variance one. The LLRs are scaled by the number of samples.
"""
function realLLRlink_row(Y,E)
    σ = groupvar(Y,E)
    return -0.5*log.(σ)
end

"""
    realLLRmed_row(Y,E,row)

Compute the log-likelihood ratios for Findr test 3 (**mediation test**) for a given `row` (gene) of gene expression matrix `Y` with categorical instrument `E` against all other rows of `Y`.

`Y` is assumed to have undergone supernormalization with each row having mean zero and variance one. The LLRs are scaled by the number of samples.
"""
function realLLRmed_row(Y,E,row)
    ρ = cov_row(Y,row)
    σ1 = groupvar(Y,E)
    σ2 = groupcov(Y,E,row)
    return -0.5*log(σ1[row].*σ1 .- (ρ .+ σ2 .- 1).^2) .+ 0.5*log(σ1[row]) + 0.5*log(1 .- ρ.^2)
end

"""
    cov_row(Y,row)

Compute the covariance between a given `row` of matrix `Y` and all other rows of `Y`.
"""
function cov_row(Y,row)
    return cov(Y[row,:],Y,dims=2,corrected=false)'
end
"""
    groupmeans(Y,E)

Compute the mean of each row of matrix `Y` for each of the groups (unique values) in categorical vector `E`.
"""
function groupmeans(Y,E)
    uE = unique(E)
    μ = zeros(size(Y)[1],length(uE))
    for i = eachindex(uE)
        μ[:,i] = mean(Y[:,E.==uE[i]], dims=2)
    end
    return μ
end

"""
    groupsizes(E)

Compute the group sizes (number of elements) for each group (unique value) in categorical vector `E`.
"""
function groupsizes(E)
    uE = unique(E)
    gs = zeros(Int,length(uE))
    for i = eachindex(uE)
        gs[i] = sum(E.==uE[i])
    end
    return gs
end

"""
    groupvar(Y,E)

Compute the weighted average variance of each row of matrix `Y` over the groups defined by the unique values of vector `E`.
"""
function groupvar(Y,E)
    ns = length(E) # number of samples
    w = pweights(groupsizes(E)/ns) # probability weights
    μ = groupmeans(Y,E)
    return 1 .- sum(μ.^2,w,dims=2)
end

"""
    groupcov_row(Y,E,row)

Compute the weighted average covariance between a given `row` of matrix `Y` and all other rows of `Y` over the groups defined by the unique values of vector `E`.
"""
function groupcov(Y,E,row)
    ns = length(E) # number of samples
    w = pweights(groupsizes(E)/ns) # probability weights
    μ = groupmeans(Y,E)
    return 1 .- sum(μ.*μ[row,:]',w,dims=2)
end