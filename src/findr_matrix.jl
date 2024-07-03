"""
    findr_matrix(X::Matrix{T}; cols=[], method="moments", combination="none") where T<:AbstractFloat

Compute posterior probabilities for nonzero pairwise correlations between columns of input matrix `X`. The probabilities are directed (asymmetric) in the sense that they are estimated from a column-specific background distribution.

The optional parameter `cols` (vector of integers) determines whether we consider all columns of `X` as source nodes (`cols=[]`, default), or only a subset of columns determined by the indices in the vector `cols`.

The optional parameter `method` determines the LLR mixture distribution fitting method and can be either `moments` (default) for the method of moments, or `kde` for kernel-based density estimation.

The optional parameter `combination` determines whether the output must be symmetrized. Possible values are `none` (default), `prod`, `mean`, or `anti`. If the optional parameter `cols` is non-empty, symmetrization makes no sense and an error will be thrown unless `combination="none"`.

See also [`findr(::DataFrame)`](@ref), [`symprobs`](@ref), [`supernormalize`](@ref), [`pprob_col`](@ref).
"""
function findr_matrix(X::Matrix{T}; cols=[], method="moments", combination="none") where T<:Real
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
    findr_matrix(X::Matrix{T},G::Matrix{S},pairGX::Matrix{S}; method="moments", combination="none") where {T<:AbstractFloat, S<:Integer}

Compute posterior probabilities for nonzero causal relations between columns of input matrix `X`. The probabilities are estimated for relations going from a subset of columns of `X` that have a (discrete) instrumental variable in input matrix `G` to all columns of `X`, while excluding self-interactions (given default value 1). The matching between columns of `X` and columns of `G` is given by `pairGX`, a two-column array where the first column corresponds to a column index in `G` and the second to a column index in `X`.

Posterior probabilities are computed for the following tests

 - Test 2 (**Linkage test**) 
 - Test 3 (**Mediation test**)
 - Test 4 (**Relevance test**)
 - Test 5 (**Pleiotropy test**)

which can be combined into the mediation test (``P_2 P_3``; `combination="mediation"`), the instrumental variable or non-independence test (``P_2 P_5``; `combination="IV"`), or BioFindr's original combination (``\\frac{1}{2}(P_2 P_5 + P_4)``; `combination="orig"`). By default, individual probability matrices for all tests are returned (`combination="none"`).

The optional parameter `method` determines the LLR mixture distribution fitting method and can be either `moments` (default) for the method of moments, or `kde` for kernel-based density estimation.

If `combination="none"`, then the output has size ncols(X) x 4 x ncols(G), where the middle index indexes the tests, and otherwise the output has size ncols(X) x ncols(G). 

See also [`findr(::DataFrame,::DataFrame,::DataFrame)`](@ref), [`supernormalize`](@ref), [`pprob_col`](@ref), [`combineprobs`](@ref).

!!! note
    `G` is currently assumed to be an array (vector or matrix) of integers. I intend to use CategoricalArrays in the future.
"""
function findr_matrix(X::Matrix{T},G::Array{S},pairGX::Matrix{R}; method="moments", combination="none") where {T<:AbstractFloat, S<:Integer, R<:Integer}
    if !(combination in Set(["none","IV","mediation","orig"]))
        error("combination parameter must be one of \"none\", \"IV\", \"mediation\", or \"orig\"")
    end
    # Inverse-normal transformation and standardization for each column of X
    Y = supernormalize(X)
    # Matrix to store posterior probabilities
    PP = ones(size(X,2),4,size(pairGX,1)) # this sets the diagonal to one
    # Compute posterior probabilities for each column separately
    Threads.@threads for col = axes(pairGX,1)
        # println(row)
        colG = pairGX[col,1]
        colX = pairGX[col,2]
        PP[Not(colX),:,col] = pprob_col(Y[:,Not(colX)], Y[:,colX], G[:,colG]; method = method)
    end
    return combineprobs(PP; combination = combination)
end

"""
    findr_matrix(X1::Matrix{T},X2::Array{T},G::Array{S},pairGX::Matrix{R}; method="moments", combination="none")  where {T<:AbstractFloat, S<:Integer}

Compute posterior probabilities for nonzero causal relations from columns of input matrix `X2` to columns of input matrix `X1`. The probabilities are estimated for a subset of columns of `X2` that have a (discrete) instrumental variable in input matrix `G`. The matching between columns of `X2` and columns of `G` is given by `pairGX`, a two-column array where the first column corresponds to a column index in `G` and the second to a column index in `X2`.

Posterior probabilities are computed for the following tests

 - Test 2 (**Linkage test**) 
 - Test 3 (**Mediation test**)
 - Test 4 (**Relevance test**)
 - Test 5 (**Pleiotropy test**)

which can be combined into the mediation test (``P_2 P_3``; `combination="mediation"`), the instrumental variable or non-independence test (``P_2 P_5``; `combination="IV"`), or BioFindr's original combination (``\\frac{1}{2}(P_2 P_5 + P_4)``; `combination="orig"`). By default, individual probability matrices for all tests are returned (`combination="none"`).

The optional parameter `method` determines the LLR mixture distribution fitting method and can be either `moments` (default) for the method of moments, or `kde` for kernel-based density estimation.

If `combination="none"`, then the output has size ncols(X1) x 4 x ncols(X2), where the middle index indexes the tests, and otherwise the output has size ncols(X1) x ncols(X2).

See also [`findr(::DataFrame,::DataFrame,::DataFrame,::DataFrame)`](@ref), [`combineprobs`](@ref), [`supernormalize`](@ref), [`pprob_col`](@ref)

!!! note
    `G` is currently assumed to be an array (vector or matrix) of integers. I intend to use CategoricalArrays in the future.
"""
function findr_matrix(X1::Matrix{T}, X2::Array{T}, G::Array{S}, pairGX::Matrix{R}; method="moments", combination="none")  where {T<:AbstractFloat, S<:Integer, R<:Integer}
    if !(combination in Set(["none","IV","mediation","orig"]))
        error("combination parameter must be one of \"none\", \"IV\", \"mediation\", or \"orig\"")
    end
    # Inverse-normal transformation and standardization for each column of X1 and X2
    Y1 = supernormalize(X1)
    Y2 = supernormalize(X2)
    # Matrix to store posterior probabilities
    PP = ones(size(X1,2),4,size(pairGX,1)) # this sets the diagonal to one
    # Compute posterior probabilities for each column separately
    Threads.@threads for col = axes(pairGX,1)
        # println(row)
        colG = pairGX[col,1]
        colX = pairGX[col,2]
        PP[:,:,col] = pprob_col(Y1, Y2[:,colX], G[:,colG]; method = method)
    end
    return combineprobs(PP; combination = combination)
end