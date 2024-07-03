#########################################
#   Methods for coexpression analysis   #
#########################################


"""
    findr(dX::T; colnames=[], method="moments", FDR=1.0, sorted=true, combination="none") where T<:AbstractDataFrame

Wrapper for `findr(Matrix(dX))` when the input `dX` is in the form of a DataFrame. The output is then also wrapped in a DataFrame with `Source`, `Target`, (Posterior) `Probability`, and `qvalue` columns.

The optional parameter `colnames` (vector of strings) determines whether we consider all columns of `dX` as source nodes (`colnames=[]`, default), or only a subset of columns determined by the variable names in the vector `colnames`.

The optional parameter `method` determines the LLR mixture distribution fitting method and can be either `moments` (default) for the method of moments, or `kde` for kernel-based density estimation.

The optional parameter `FDR` can be used to return only a subset of interactions with a desired expected FDR value (q-value threshold) (default 1.0, no filtering).

The optional parameter `sorted` determines if the output must be sorted by increasing q-value / decreasing posterior probability (`sorted=true`, the default) or by causal factor (column names of `dX`) (`sorted=false`)

The optional parameter `combination` determines whether the output must be symmetrized. Possible values are `none` (default), `prod`, `mean`, or `anti`. If the optional parameter `colnames` is non-empty, symmetrization makes no sense and an error will be thrown unless `combination="none"`.

See also [`findr(::Matrix)`](@ref), [`symprobs`](@ref), [`stackprobs`](@ref), [`globalfdr!`](@ref).
"""
function findr(dX::T; colnames=[], method="moments", FDR=1.0, sorted=true, combination="none") where T<:AbstractDataFrame
    if !isempty(colnames)
        colnames = intersect(colnames, names(dX))
        cols = indexin(colnames, names(dX))
    else
        colnames = names(dX)
        cols = []
    end
    PP = findr(Matrix(dX); cols = cols, method = method, combination = combination)
    dP = stackprobs(PP, colnames, names(dX))
    globalfdr!(dP, FDR = FDR, sorted = sorted)
    return dP
end

"""
    findr_matrix(X::Matrix{T}; cols=[], method="moments", combination="none") where T<:AbstractFloat

Compute posterior probabilities for nonzero pairwise correlations between columns of input matrix `X`. The probabilities are directed (asymmetric) in the sense that they are estimated from a column-specific background distribution.

The optional parameter `cols` (vector of integers) determines whether we consider all columns of `X` as source nodes (`cols=[]`, default), or only a subset of columns determined by the indices in the vector `cols`.

The optional parameter `method` determines the LLR mixture distribution fitting method and can be either `moments` (default) for the method of moments, or `kde` for kernel-based density estimation.

The optional parameter `combination` determines whether the output must be symmetrized. Possible values are `none` (default), `prod`, `mean`, or `anti`. If the optional parameter `cols` is non-empty, symmetrization makes no sense and an error will be thrown unless `combination="none"`.

See also [`findr(::DataFrame)`](@ref), [`symprobs`](@ref), [`supernormalize`](@ref), [`pprob_col`](@ref).
"""
function findr_matrix(X::Matrix{T}; cols=[], method="moments", combination="none") where T<:AbstractFloat
    # Inverse-normal transformation and standardization for each columns of X
    Y = supernormalize(X)
    # check if we need to use all columns or only a subset as source nodes
    if isempty(cols)
        cols = axes(Y,2) # use all columns
    end
    # Matrix to store posterior probabilities
    nall = size(Y,2)
    ncols = length(cols)
    PP = ones(nall,ncols) # this sets the diagonal elements to one
    # Compute posterior probabilities for each column separately
    Threads.@threads for i in eachindex(cols) 
        PP[Not(cols[i]),i] = pprob_col(Y[:,Not(cols[i])],Y[:,cols[i]]; method = method)
    end
    return symprobs(PP, combination = combination)
end

"""
    findr_matrix(X1::Matrix{T}, X2::Matrix{T}; method="moments") where T<:AbstractFloat

Compute posterior probabilities for nonzero pairwise correlations between columns of input matrices `X1` and `X2`. The probabilities are directed (asymmetric) from the columns of `X2` to the columns of `X1` in the sense that they are estimated from a column-specific background distribution for each column of `X2`.

The optional parameter `method` determines the LLR mixture distribution fitting method and can be either `moments` (default) for the method of moments, or `kde` for kernel-based density estimation.

Only use this method if `X1` and `X2` are distinct (no overlapping columns). For `X2` consisting of a subset of columns with indices `idx`, use `findr(X1; cols=idx)` instead. 

See also [`findr(::DataFrame)`](@ref), [`symprobs`](@ref), [`supernormalize`](@ref), [`pprob_col`](@ref).
"""
function findr_matrix(X1::Matrix{T}, X2::Array{T}; method="moments") where T<:AbstractFloat
    # Inverse-normal transformation and standardization for each columns of X1 and X2
    Y1 = supernormalize(X1)
    Y2 = supernormalize(X2)
    # Matrix to store posterior probabilities
    PP = ones(size(X1,2), size(X2,2))
    # Compute posterior probabilities for each column separately
    Threads.@threads for col = axes(Y2,2)
        PP[:,col] = pprob_col(Y1,Y2[:,col]; method = method)
    end
    return PP
end
  