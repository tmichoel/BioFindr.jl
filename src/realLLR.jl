"""
    realLLR_col(Y::Matrix{T},Ycol::Vector{T}) where T<:AbstractFloat

Compute the log-likelihood ratios for Findr test 0 (**correlation test**) for a given column vector `Ycol` against all columns of matrix `Y`.

`Y` and `Ycol` are assumed to have undergone supernormalization with each column having mean zero and variance one. The LLRs are scaled by the number of rows (samples).

See also [`supernormalize`](@ref).
"""
function realLLR_col(Y::Matrix{T},Ycol::Vector{T}) where T<:AbstractFloat
    ρ = vec(cov(Y,Ycol,corrected=false))
    # ρ[col] = 1. # set self to exact value
    -0.5*log.(abs.(1 .- ρ.^2))
end

"""
    realLLR_col(Y::Matrix{T},Ycol::Vector{T},E::Vector{S}) where {T<:AbstractFloat, S<:Integer}

Compute the log-likelihood ratios for the Findr causal tests for a given column vector `Ycol` with categorical instrument `E` against all columns of matrix `Y` : 

- Test 2 (**Linkage test**) 
- Test 3 (**Mediation test**)
- Test 4 (**Relevance test**)
- Test 5 (**Pleiotropy test**)

`Y` and `Ycol` are assumed to have undergone supernormalization with each column having mean zero and variance one. The LLRs are scaled by the number of rows (samples).

See also [`supernormalize`](@ref), [`llrstats_col`](@ref).
"""
function realLLR_col(Y::Matrix{T},Ycol::Vector{T},E::Vector{S}) where {T<:AbstractFloat, S<:Integer}
    # compute the sufficient statistics
    ρ, σ, σcol = llrstats_col(Y,Ycol,E)

    # test 2
    llr2 = -0.5*log.(σ[:,1])

    # test 4
    # the abs in the argument of the log is to have a positive argument in cases where the exact value is 0
    llr4 = -0.5*log.(abs.(σcol.*σ[:,1] .- (ρ .+ σ[:,2] .- 1).^2)) .+ 0.5*log.(σcol)
    #llr4[col] = Inf # set self to Inf

    # test 3
    llr3 = llr4 .+ 0.5*log.(abs.(1 .- ρ.^2))
    # llr3[col] = Inf

    # test 5
    llr5 = llr4 .+ 0.5*log.(σ[:,1])
    # llr5[col] = Inf

    # output
    llr2, llr3, llr4, llr5
end


"""
    realLLR_col(Y::Matrix{T},E::Vector{S}) where {T<:AbstractFloat, S<:Integer}

Compute the log-likelihood ratios for Findr test 2 (**Linkage test**)  for a given categorical vector `E` against all columns of matrix `Y`.

`Y` is assumed to have undergone supernormalization with each column having mean zero and variance one. The LLRs are scaled by the number of rows (samples).

See also [`supernormalize`](@ref), [`llrstats_col`](@ref).
"""
function realLLR_col(Y::Matrix{T},E::Vector{S}) where {T<:AbstractFloat, S<:Integer}
    # compute the sufficient statistics
    σ = llrstats_col(Y,E)

    # test 2
    -0.5*log.(σ)
end

"""
    llrstats_col(Y,Ycol,E)

Compute the sufficient statistics to compute the log-likelihood ratios for Findr tests 2-5  for a given column vector `Ycol` (gene) of gene expression values with categorical instrument `E` against all columns of `Y`.

`Ycol` and `Y` are assumed to have undergone supernormalization with each col having mean zero and variance one. The LLRs are scaled by the number of rows (samples).

The sufficient statistics are:

- the covariance `ρ` between the given `Ycol` and all columns of `Y`
- the weighted average variances `σ1` of each column of matrix `Y` over the groups (unique values) in `E`
- the weighted average covariance `σ2` between `Ycol` and all  columns of `Y` over the groups (unique values) of `E`

See also [`groupmeans`](@ref).
"""
function llrstats_col(Y,Ycol,E)
    ρ = vec(cov(Y,Ycol,corrected=false))

    gs, μ, μcol = groupmeans(Y,Ycol,E)
    
    ns = length(E) # number of samples
    w = pweights(gs/ns) # probability weights

    σ = [vec(1 .- sum(μ.^2,w,dims=2)) vec(1 .- sum(μ.*μcol',w,dims=2))]
    σcol = 1 - sum(μcol.^2,w)

    ρ, σ, σcol
end


"""
    llrstats_col(Y,E)

Compute the sufficient statistics to compute the log-likelihood ratios for Findr tests 2  for a given categorical vector `E` against all columns of `Y`.

`Y` is assumed to have undergone supernormalization with each col having mean zero and variance one. The LLRs are scaled by the number of rows (samples).

The sufficient statistics are the weighted average variances `σ1` of each column of matrix `Y` over the groups (unique values) in `E`.

See also [`groupmeans`](@ref).
"""
function llrstats_col(Y,E)
    gs, μ = groupmeans(Y,E)
    
    ns = length(E) # number of samples
    w = pweights(gs/ns) # probability weights

    vec(1 .- sum(μ.^2,w,dims=2))
end

"""
    groupmeans(Y,Ycol,E)

Compute the size and mean of each column of matrix `Y` and of the column vector `Ycol` for each of the groups (unique values) in categorical vector `E`.
"""
function groupmeans(Y,Ycol,E)
    uE = unique(E)
    gs = zeros(Int,length(uE))
    μ = zeros(size(Y,2),length(uE))
    μcol = zeros(length(uE))
    Threads.@threads for i = eachindex(uE)
        gs[i] = sum(E.==uE[i])
        μ[:,i] = mean(Y[E.==uE[i],:], dims=1)
        μcol[i] = mean(Ycol[E.==uE[i]])
    end
    gs, μ, μcol
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