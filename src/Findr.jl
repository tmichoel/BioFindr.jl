module Findr

# External packages

using DataFrames
using Statistics
using StatsBase
using Distributions
using Random
using SpecialFunctions
using KernelDensity
using LinearAlgebra
using InvertedIndices

# Define test names as variables
corr = "corr"
link = "link"
med = "med"
relev = "relev"
pleio = "pleio"

# Load the different code files

include("supernormalization.jl")

include("lbeta.jl")

include("randomLLR.jl")

include("realLLR.jl")

include("posteriorprobs.jl")

include("bayesiannets.jl")

# Main Findr function calls

"""
    findr(X)

Compute posterior probabilities for nonzero pairwise correlations between columns of input matrix `X`. The probabilities are directed (asymmetric) in the sense that they are estimated from a column-specific background distribution. 
"""
function findr(X::Matrix{T}) where T<:AbstractFloat
    # Inverse-normal transformation and standardization for each columns of X
    Y = supernormalize(X)
    # Matrix to store posterior probabilities
    ncols = size(Y,2)
    PP = ones(ncols,ncols) # this sets the diagonal elements to one
    # Compute posterior probabilities for each column separately
    Threads.@threads for col = axes(Y,2)
        PP[Not(col),col] = pprob_col(Y[:,Not(col)],Y[:,col])
    end
    return PP
end

"""
    findr(X,G)

Compute posterior probabilities for nonzero differential expression of colunns of input matrix `X` across groups defined by one or more categorical variables (columns of `G`).

Return a matrix of size ncols(X) x ncols(G)
"""
function findr(X::Matrix{T},G::Array{S}) where {T<:AbstractFloat, S<:Integer}
    # Inverse-normal transformation and standardization for each column of X
    Y = supernormalize(X)
    # Matrix to store posterior probabilities
    PP = zeros(size(X,2),size(G,2))
    # Compute posterior probabilities for each column separately
    Threads.@threads for col = axes(G,2)
        PP[:,col] = pprob_col(Y,G[:,col])
    end
    return PP
end

"""
    findr(X,G,pairGX)

Compute posterior probabilities for nonzero causal relations between columns of input matrix `X`. The probabilities are estimated for a subset of columns of `X` that have a (discrete) instrumental variable in input matrix `G`. The matching between columns of `X` and columns of `G` is given by `pairGX`, a two-column array where the first column corresponds to a column index in `G` and the second to a column index in `X`.

Posterior probabilities are computed for the following tests

 - Test 2 (**Linkage test**) 
 - Test 3 (**Mediation test**)
 - Test 4 (**Relevance test**)
 - Test 5 (**Pleiotropy test**)

which can be combined into the mediation test (``P_2 P_3``), the non-independence test (``P_2 P_5``), or Findr's legacy combination (``\\frac{1}{2}(P_2 P_5 + P_4)``). Alternatively, individual probability matrices for all tests are returned.

All return matrices have size ncols(X) x ncols(G) .
"""
function findr(X::Matrix{T},G::Matrix{S},pairGX::Matrix{S}) where {T<:AbstractFloat, S<:Integer}
    # Inverse-normal transformation and standardization for each column of X
    Y = supernormalize(X)
    # Matrix to store posterior probabilities
    PP = ones(size(X,2),4,size(pairGX,1)) # this sets the diagonal to one
    # Compute posterior probabilities for each column separately
    Threads.@threads for col = axes(pairGX,1)
        # println(row)
        colG = pairGX[col,1]
        colX = pairGX[col,2]
        PP[Not(colX),:,col] = pprob_col(Y[:,Not(colX)], Y[:,colX], G[:,colG])
    end
    return PP
end

"""
    findr(X1,X2,G)

Compute posterior probabilities for nonzero causal relations from columns of input matrix `X2` to columns of input matrix `X1`. The columns of input matrix `G` are (discrete) instrumental variables for the corresponding columns in `X2`.

Posterior probabilities are computed for the following tests

 - Test 2 (**Linkage test**) 
 - Test 3 (**Mediation test**)
 - Test 4 (**Relevance test**)
 - Test 5 (**Pleiotropy test**)

which can be combined into the mediation test (``P_2 P_3``), the non-independence test (``P_2 P_5``), or Findr's legacy combination (``\\frac{1}{2}(P_2 P_5 + P_4)``). Alternatively, individual probability matrices for all tests are returned.

All return matrices have size ncols(X1) x ncols(X2) .
"""
function findr(X1::Matrix{T},X2::Array{T},G::Array{S})  where {T<:AbstractFloat, S<:Integer}
    # Inverse-normal transformation and standardization for each column of X1 and X2
    Y1 = supernormalize(X1)
    Y2 = supernormalize(X2)
    # Matrix to store posterior probabilities
    PP = ones(size(X1,2),4,size(X2,2)) # this sets the diagonal to one
    # Compute posterior probabilities for each column separately
    Threads.@threads for col = axes(Y2,2)
        # println(row)
        PP[:,:,col] = pprob_col(Y1, Y2[:,col], G[:,col])
    end
    return PP
end

end
