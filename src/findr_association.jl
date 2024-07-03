####################################################################
#   Methods for association and differential expression analysis   #
####################################################################


"""
    findr_matrix(X::Matrix{T},G::Array{S}; method="moments") where {T<:AbstractFloat, S<:Integer}

Compute posterior probabilities for nonzero differential expression of colunns of input matrix `X` across groups defined by one or more categorical variables (columns of `G`).

Return a matrix of size ncols(X) x ncols(G)

The optional parameter `method` determines the LLR mixture distribution fitting method and can be either `moments` (default) for the method of moments, or `kde` for kernel-based density estimation.

See also [`findr(::DataFrame,::DataFrame)`](@ref), [`supernormalize`](@ref), [`pprob_col`](@ref).

!!! note
    `G` is currently assumed to be an array (vector or matrix) of integers. CategoricalArrays will be supported in the future.
"""
function findr_matrix(X::Matrix{T}, G::Array{S}; method="moments") where {T<:AbstractFloat, S<:Integer}
    # Inverse-normal transformation and standardization for each column of X
    Y = supernormalize(X)
    # Matrix to store posterior probabilities
    PP = zeros(size(X,2),size(G,2))
    # Compute posterior probabilities for each column separately
    Threads.@threads for col = axes(G,2)
        PP[:,col] = pprob_col(Y,G[:,col]; method = method)
    end
    return PP
end


"""
    findr(dX::T, dG::T; method="moments", FDR=1.0, sorted=true) where T<:AbstractDataFrame

Wrapper for `findr(Matrix(dX), Matrix(dG))` when the inputs `dX` and `dG` are in the form of a DataFrame. The output is then also wrapped in a DataFrame with `Source`, `Target` (Posterior) `Probability`, and `qvalue` columns.

The optional parameter `method` determines the LLR mixture distribution fitting method and can be either `moments` (default) for the method of moments, or `kde` for kernel-based density estimation.

The optional parameter `FDR` can be used to return only a subset of interactions with a desired expected FDR value (q-value threshold) (default 1.0, no filtering).

The optional parameter `sorted` determines if the output must be sorted by increasing q-value / decreasing posterior probability (`sorted=true`, the default) or by causal factor (column names of `dX`) (`sorted=false`)

Note that depending on the type of `Matrix(dG)`, different matrix-based methods are called. If `Matrix(dG)` consists of Floats, posterior probabilities for nonzero pairwise correlations between the variables in `dG` and variables in `dX` are computed. If `Matrix(dG)` consists of integers, posterior probabilities for nonzero differential expression of variables in `dX` across groups defined by the variables in `dG` are computed

See also [`findr(::Matrix,::Array)`], [`stackprobs`](@ref), [`globalfdr!`](@ref). 
"""
function findr(dX::T, dG::T; method="moments", FDR=1.0, sorted=true) where T<:AbstractDataFrame
    dP = stackprobs(findr(Matrix(dX), Matrix(dG); method = method), names(dG), names(dX))
    globalfdr!(dP, FDR = FDR, sorted = sorted)
    return dP
end
