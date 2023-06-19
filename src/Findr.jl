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

include("utils.jl")

# Main Findr function calls

"""
    findr(X::Matrix{T}) where T<:AbstractFloat

Compute posterior probabilities for nonzero pairwise correlations between columns of input matrix `X`. The probabilities are directed (asymmetric) in the sense that they are estimated from a column-specific background distribution. 

See also [`findr(dX::DataFrame)`](@ref)-
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
    findr(dX::T) where T<:AbstractDataFrame

Wrapper for `findr(Matrix(dX))` when the input `dX` is in the form of a DataFrame. The output is then also wrapped in a DataFrame with `Source`, `Target` and `Posterior probability` columns.

See also [`stackprobs`](@ref).
"""
function findr(dX::T) where T<:AbstractDataFrame
    stackprobs(findr(Matrix(dX)), names(dX), names(dX))
end

"""
    findr(X::Matrix{T},G::Array{S}) where {T<:AbstractFloat, S<:Integer}

Compute posterior probabilities for nonzero differential expression of colunns of input matrix `X` across groups defined by one or more categorical variables (columns of `G`).

Return a matrix of size ncols(X) x ncols(G)

See also [`findr(dX::DataFrame,dG::DataFrame)`](@ref).

!!! note
    `G` is currently assumed to be an array (vector or matrix) of integers. CategoricalArrays will be supported in the future.
"""
function findr(X::Matrix{T}, G::Array{S}) where {T<:AbstractFloat, S<:Integer}
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
    findr(dX::T, dG::T) where T<:AbstractDataFrame

Wrapper for `findr(Matrix(dX), Matrix(dG))` when the inputs `dX` and `dG` are in the form of a DataFrame. The output is then also wrapped in a DataFrame with `Source`, `Target` and `Posterior probability` columns.

See also [`stackprobs`](@ref). 
"""
function findr(dX::T, dG::T) where T<:AbstractDataFrame
    stackprobs(findr(Matrix(dX), Matrix(dG)), names(dG), names(dX))
end


"""
    findr(X::Matrix{T},G::Matrix{S},pairGX::Matrix{S}; combination="none") where {T<:AbstractFloat, S<:Integer}

Compute posterior probabilities for nonzero causal relations between columns of input matrix `X`. The probabilities are estimated for a subset of columns of `X` that have a (discrete) instrumental variable in input matrix `G`. The matching between columns of `X` and columns of `G` is given by `pairGX`, a two-column array where the first column corresponds to a column index in `G` and the second to a column index in `X`.

Posterior probabilities are computed for the following tests

 - Test 2 (**Linkage test**) 
 - Test 3 (**Mediation test**)
 - Test 4 (**Relevance test**)
 - Test 5 (**Pleiotropy test**)

which can be combined into the mediation test (``P_2 P_3``; `combination="mediation"`), the instrumental variable or non-independence test (``P_2 P_5``; `combination="IV"`), or Findr's original combination (``\\frac{1}{2}(P_2 P_5 + P_4)``; `combination="orig"`). By default, individual probability matrices for all tests are returned (`combination="none"`).

All return matrices have size ncols(X) x ncols(G).

!!! note
    `G` is currently assumed to be an array (vector or matrix) of integers. I intend to use CategoricalArrays in the future.
"""
function findr(X::Matrix{T},G::Matrix{S},pairGX::Matrix{R}; combination="none") where {T<:AbstractFloat, S<:Integer, R<:Integer}
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
        PP[Not(colX),:,col] = pprob_col(Y[:,Not(colX)], Y[:,colX], G[:,colG])
    end
    return combineprobs(PP; combination = combination)
end

"""
    findr(dX::T, dG::T, dE::T, colX=2, colG=1, combination="IV") where T<:AbstractDataFrame

Wrapper for `findr(X, G, pairGX)` when the inputs are in the form of a DataFrame. The output is then also wrapped in a DataFrame with `Source`, `Target` and `Posterior probability` columns. When DataFrames are used, only combined posterior probabilities can be returned (`combination="IV"` (default), `"mediation"`, or `"orig"`).

The input dataframes are:

- `dX` - DataFrame with expression data, columns are genes
- `dG` - DataFrame with genotype data, columns are variants (SNPs)
- `dE` - DataFrame with eQTL results, must contains columns with gene and SNP IDs that can be mapped to column names in `dX` and `dG`, respectively
- `colG` - name or number of variant ID column in `dE`, default 1
- `colX` - name or number of gene ID column in `dE`, default 2

See also [`stackprobs`](@ref).
"""
function findr(dX::T, dG::T, dE::T; colG=1, colX=2, combination="IV") where T<:AbstractDataFrame
    if combination == "none"
        error("Returning posterior probabilities for individual tests is not supported with DataFrame inputs. Set combination argument to one of \"IV\", \"mediation\", or \"orig\", or use matrix inputs.")
    elseif combination in Set(["IV","mediation","orig"])
        # Create the array with SNP-Gene pairs
        pairGX = getpairs(dX, dG, dE; colG = colG, colX = colX)
        # Call Findr on numeric data
        PP = findr(Matrix(dX), Matrix(dG), pairGX; combination=combination)
        return stackprobs(PP, names(dX)[pairGX[:,2]], names(dX))   
    else
        error("Combination parameter must be one of \"IV\", \"mediation\", or \"orig\"")
    end
end



"""
    findr(X1::Matrix{T},X2::Array{T},G::Array{S},pairGX::Matrix{R}; combination="none")  where {T<:AbstractFloat, S<:Integer}

Compute posterior probabilities for nonzero causal relations from columns of input matrix `X2` to columns of input matrix `X1`. The probabilities are estimated for a subset of columns of `X2` that have a (discrete) instrumental variable in input matrix `G`. The matching between columns of `X2` and columns of `G` is given by `pairGX`, a two-column array where the first column corresponds to a column index in `G` and the second to a column index in `X2`.

Posterior probabilities are computed for the following tests

 - Test 2 (**Linkage test**) 
 - Test 3 (**Mediation test**)
 - Test 4 (**Relevance test**)
 - Test 5 (**Pleiotropy test**)

which can be combined into the mediation test (``P_2 P_3``; `combination="mediation"`), the instrumental variable or non-independence test (``P_2 P_5``; `combination="IV"`), or Findr's original combination (``\\frac{1}{2}(P_2 P_5 + P_4)``; `combination="orig"`). By default, individual probability matrices for all tests are returned (`combination="none"`).

All return matrices have size ncols(X1) x ncols(X2).

!!! note
    `G` is currently assumed to be an array (vector or matrix) of integers. I intend to use CategoricalArrays in the future.
"""
function findr(X1::Matrix{T}, X2::Array{T}, G::Array{S}, pairGX::Matrix{R}; combination="none")  where {T<:AbstractFloat, S<:Integer, R<:Integer}
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
        PP[:,:,col] = pprob_col(Y1, Y2[:,colX], G[:,colG])
    end
    return combineprobs(PP; combination = combination)
end


"""
    findr(dX1::T, dX2::T, dG::T, dE::T; colG=1, colX=2, combination="IV") where T<:AbstractDataFrame

Wrapper for `findr(Matrix(dX1), Matrix(dX2), Matrix(dG))` when the inputs `dX1`, `dX2`, and `dG` are in the form of a DataFrame. The output is then also wrapped in a DataFrame with `Source`, `Target` and `Posterior probability` columns. When DataFrames are used, only a combined posterior probabilities can be returned (`combination="IV"` (default), `"mediation"`, or `"orig"`).

See also [`stackprobs`](@ref).
"""
function findr(dX1::T, dX2::T, dG::T, dE::T; colG=1, colX=2, combination="IV") where T<:AbstractDataFrame
    if combination == "none"
        error("Returning posterior probabilities for individual tests is not supported with DataFrame inputs. Set combination argument to one of \"IV\", \"mediation\", or \"orig\", or use matrix inputs.")
    elseif combination in Set(["IV","mediation","orig"])
        # Create the array with SNP-Gene pairs
        pairGX = getpairs(dX2, dG, dE; colG = colG, colX = colX)
        # Call Findr on numeric data
        PP = findr(Matrix(dX1), Matrix(dX2), Matrix(dG), pairGX; combination=combination)
        return stackprobs(PP, names(dX2)[pairGX[:,2]], names(dX1))
    else
        error("Combination parameter must be one of \"IV\", \"mediation\", or \"orig\"")
    end
end

end
