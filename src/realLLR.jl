"""
    realLLRcorr_row(Y,row)

Compute the log-likelihood ratios for Findr test 0 (**correlation test**) for a given `row` (gene) of gene expression matrix `Y` against all other rows of `Y`.

`Y` is assumed to have undergone supernormalization with each row having mean zero and variance one. The LLRs are scaled by the number of samples.
"""
function realLLRcorr_row(Y,row)
    ρ = cov_row(Y,row)
    -0.5*log.(1 .- ρ.^2)
end

"""
    realLLRcorr_row(ρ)

Compute the log-likelihood ratios for Findr test 0 (**correlation test**) with precomputed correlation values `̢ρ`
"""
realLLRcorr_row(ρ) = -0.5*log.(1 .- ρ.^2)

"""
    realLLRlink_row(Y,E)

Compute the log-likelihood ratios for Findr test 2 (**linkage test**) for categorical vector `E` against all rows of matrix `Y`.

`Y` is assumed to have undergone supernormalization with each row having mean zero and variance one. The LLRs are scaled by the number of samples.
"""
function realLLRlink_row(Y,E)
    σ = groupvar(Y,E)
    -0.5*log.(σ)
end

"""
    realLLRlink_row(σ)

Compute the log-likelihood ratios for Findr test 2 (**linkage test**) with precomputed weighted average group variances `σ`.
"""
realLLRlink_row(σ) = -0.5*log.(σ)

"""
    realLLRmed_row(Y,E,row)

Compute the log-likelihood ratios for Findr test 3 (**mediation test**) for a given `row` (gene) of gene expression matrix `Y` with categorical instrument `E` against all other rows of `Y`.

`Y` is assumed to have undergone supernormalization with each row having mean zero and variance one. The LLRs are scaled by the number of samples.
"""
function realLLRmed_row(Y,E,row)
    ρ = cov_row(Y,row)
    σ1 = groupvar(Y,E)
    σ2 = groupcov(Y,E,row)
    -0.5*log(σ1[row].*σ1 .- (ρ .+ σ2 .- 1).^2) .+ 0.5*log(σ1[row]) + 0.5*log(1 .- ρ.^2)
end

"""
    realLLRmed_row(ρ,σ1,σ2)

Compute the log-likelihood ratios for Findr test 3 (**mediation test**) with precomputed correlation values `ρ` and weighted average group variances `σ1` and covariances `σ2`
"""
realLLRmed_row(ρ,σ1,σ2) = -0.5*log(σ1[row].*σ1 .- (ρ .+ σ2 .- 1).^2) .+ 0.5*log(σ1[row]) + 0.5*log(1 .- ρ.^2)

"""
    realLLRrelev_row(Y,E,row)

Compute the log-likelihood ratios for Findr test 4 (**relevance test**) for a given `row` (gene) of gene expression matrix `Y` with categorical instrument `E` against all other rows of `Y`.

`Y` is assumed to have undergone supernormalization with each row having mean zero and variance one. The LLRs are scaled by the number of samples.
"""
function realLLRrelev_row(Y,E,row)
    ρ = cov_row(Y,row)
    σ1 = groupvar(Y,E)
    σ2 = groupcov(Y,E,row)
    -0.5*log(σ1[row].*σ1 .- (ρ .+ σ2 .- 1).^2) .+ 0.5*log(σ1[row])
end

"""
    realLLRrelev_row(ρ,σ1,σ2)

Compute the log-likelihood ratios for Findr test 3 (**mediation test**) with precomputed correlation values `ρ` and weighted average group variances `σ1` and covariances `σ2`
"""
realLLRrelev_row(ρ,σ1,σ2) = -0.5*log(σ1[row].*σ1 .- (ρ .+ σ2 .- 1).^2) .+ 0.5*log(σ1[row])


"""
    realLLRpleio_row(Y,E,row)

Compute the log-likelihood ratios for Findr test 5 (**pleiotropy test**) for a given `row` (gene) of gene expression matrix `Y` with categorical instrument `E` against all other rows of `Y`.

`Y` is assumed to have undergone supernormalization with each row having mean zero and variance one. The LLRs are scaled by the number of samples.
"""
function realLLRpleio_row(Y,E,row)
    ρ = cov_row(Y,row)
    σ1 = groupvar(Y,E)
    σ2 = groupcov(Y,E,row)
    -0.5*log(σ1[row].*σ1 .- (ρ .+ σ2 .- 1).^2) .+ 0.5*log(σ1[row].*σ1)
end

"""
    realLLRpleio_row(ρ,σ1,σ2)

Compute the log-likelihood ratios for Findr test 3 (**mediation test**) with precomputed correlation values `ρ` and weighted average group variances `σ1` and covariances `σ2`
"""
realLLRpleio_row(ρ,σ1,σ2) = -0.5*log(σ1[row].*σ1 .- (ρ .+ σ2 .- 1).^2) .+ 0.5*log(σ1[row].*σ1)


"""
    llrstats_row(Y,E,row)

Compute the sufficient statistics to compute the log-likelihood ratios for Findr tests 2-5  for a given `row` (gene) of gene expression matrix `Y` with categorical instrument `E` against all other rows of `Y`.

`Y` is assumed to have undergone supernormalization with each row having mean zero and variance one. The LLRs are scaled by the number of samples.

The sufficient statistics are:

- the group sizes (number of elements) for each group (unique value) in `E`
- the mean of each row of matrix `Y` for each of the groups (unique values) in `E`
- the weighted average variance of each row of matrix `Y` over the groups (unique values) in `E`
- the weighted average covariance between the given `row` of `Y` and all other rows of `Y` over the groups (unique values) of `E`
"""
function llrstats_row(Y,E,row)
    ρ = cov_row(Y,row)
    σ1 = groupvar(Y,E)
    σ2 = groupcov(Y,E,row)
    ρ, σ1, σ2
end

"""
    cov_row(Y,row)

Compute the covariance between a given `row` of matrix `Y` and all other rows of `Y`.
"""
function cov_row(Y,row)
    return cov(Y[row,:],Y,dims=2,corrected=false)'
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