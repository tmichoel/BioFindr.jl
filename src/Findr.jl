module Findr

# Write your package code here.

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

"""
    findr_corr(X)

Compute posterior probabilities for nonzero pairwise correlations between rows of input matrix `X`. The probabilities are directed in the sense that they are estimated from a row-specific background distribution. 
"""
function findr_corr(X)
    # Inverse-normal transformation and standardization for each row of X
    Y = supernormalize(X)
    # Matrix to store posterior probabilities
    ncols = size(X,2)
    PP = zeros(ncols,ncols)
    # Compute posterior probabilities for each row separately
    Threads.@threads for col = 1:ncols
        PP[:,col] = pprob_corr_row(Y,col)
    end
    return PP
end

"""
    findr_diffexp(X,E)

Compute posterior probabilities for nonzero differential expression of rows of input matrix `X` across groups defined by one or more categorical variables (rows of `E`).

Return a matrix of size num_rows(E) x num_rows(X)
"""
function findr_diffexp(X,E)

end

"""
    findr_causal(X,E,pairEX)

Compute posterior probabilities for nonzero causal relations between rows of input matrix `X`. The probabilities are estimated for a subset of rows of `X` that have a (discrete) instrumental variable in input matrix `E`. The matching between rows of `X` and rows of `E` is given by `pairEX`, a two-column array where the first column corresponds to a row index in `E` and the second to a row index in `X`.

Posterior probabilities are computed for the following tests

 - Test 2 (**Linkage test**) 
 - Test 3 (**Mediation test**)
 - Test 4 (**Relevance test**)
 - Test 5 (**Pleiotropy test**)

which can be combined into the mediation test (P2*P3), the non-independence test (P2*P5), or Findr's legacy combination (0.5*(P2*P5 + P4)). Alternatively, individual probability matrices for all tests are returned.

All return matrices have size num_rows(E) x num_rows(X).
"""
function findr_causal(X,G,pairGX)
    # Inverse-normal transformation and standardization for each row of X
    Y = supernormalize(X)
    # Matrix to store posterior probabilities
    npairs = size(pairGX,1)
    ncolsX = size(X,1)
    PP2 = zeros(ncolsX,npairs)
    PP3 = zeros(ncolsX,npairs)
    PP4 = zeros(ncolsX,npairs)
    PP5 = zeros(ncolsX.npairs)
    # Compute posterior probabilities for each row separately
    Threads.@threads for col = 1:npairs
        # println(row)
        colG = pairGX[col,1]
        colX = pairGX[col,2]
        PP2[:,col], PP3[:,col], PP4[:,col], PP5[:,col] = pprob_causal_row(Y,G[:,colG],colX)
    end
    return PP2, PP3, PP4, PP5
end


end
