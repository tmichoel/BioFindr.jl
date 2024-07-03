#####################################
#   Methods for causal inference    #
#####################################


"""
    findr(dX::T, dG::T, dE::T; colX=2, colG=1, method="moments", combination="IV", FDR=1.0, sorted=true) where T<:AbstractDataFrame

Wrapper for `findr(Matrix(dX), Matrix(dG), pairGX)` when the inputs are in the form of a DataFrame. The output is then also wrapped in a DataFrame with `Source`, `Target` (Posterior) `Probability`, and `qvalue` columns. When DataFrames are used, only combined posterior probabilities can be returned (`combination="IV"` (default), `"mediation"`, or `"orig"`).

The input dataframes are:

- `dX` - DataFrame with expression data, columns are genes
- `dG` - DataFrame with genotype data, columns are variants (SNPs)
- `dE` - DataFrame with eQTL results, must contain columns with gene and SNP IDs that can be mapped to column names in `dX` and `dG`, respectively

The numeric mapping between column indices in `Matrix(dG)` and `Matrix(dX)` is obtained from these inputs using the [`getpairs`](@ref) function and the optional parameters:

- `colG` - name or number of variant ID column in `dE`, default 1
- `colX` - name or number of gene ID column in `dE`, default 2
- `namesX` - names of a possible subset of columns in `dX` to be considered as potential causal regulators (default `names(dX)`)

The optional parameter `method` determines the LLR mixture distribution fitting method and can be either `moments` (default) for the method of moments, or `kde` for kernel-based density estimation.

The optional parameter `FDR` can be used to return only a subset of interactions with a desired expected FDR value (q-value threshold) (default 1.0, no filtering).

The optional parameter `sorted` determines if the output must be sorted by increasing q-value / decreasing posterior probability (`sorted=true`, the default) or by causal factor (column names of `dX`) (`sorted=false`)

See also [`findr(::Matrix,::Array,::Matrix)`](@ref), [`getpairs`](@ref), [`combineprobs`](@ref), [`stackprobs`](@ref), [`globalfdr!`](@ref).
"""
function findr(dX::T, dG::T, dE::T; colG=1, colX=2, namesX=[], method="moments", combination="IV", FDR=1.0, sorted=true) where T<:AbstractDataFrame
    if combination == "none"
        error("Returning posterior probabilities for individual tests is not supported with DataFrame inputs. Set combination argument to one of \"IV\", \"mediation\", or \"orig\", or use matrix inputs.")
    elseif combination in Set(["IV","mediation","orig"])
        # Create the array with SNP-Gene pairs
        pairGX = getpairs(dX, dG, dE; colG = colG, colX = colX, namesX = namesX)
        # Call BioFindr on numeric data
        PP = findr_matrix(Matrix(dX), Matrix(dG), pairGX; method = method, combination = combination)
        dP = stackprobs(PP, names(dX)[pairGX[:,2]], names(dX)) 
        globalfdr!(dP, FDR = FDR, sorted = sorted)
        return dP
    else
        error("Combination parameter must be one of \"IV\", \"mediation\", or \"orig\"")
    end
end

"""
    findr(dX1::T, dX2::T, dG::T, dE::T; colG=1, colX=2, method="moments", combination="IV", FDR=1.0, sorted=true) where T<:AbstractDataFrame

Wrapper for `findr(Matrix(dX1), Matrix(dX2), Matrix(dG), pairGX2)` when the inputs `dX1`, `dX2`, and `dG` are in the form of a DataFrame. The output is then also wrapped in a DataFrame with `Source`, `Target`, (Posterior) `Probability`, and `qvalue` columns. When DataFrames are used, only combined posterior probabilities can be returned (`combination="IV"` (default), `"mediation"`, or `"orig"`).

The numeric mapping between column indices in `Matrix(dG)` and `Matrix(dX2)` is obtained from these inputs using the [`getpairs`](@ref) function and the optional parameters:

- `colG` - name or number of variant ID column in `dE`, default 1
- `colX` - name or number of gene ID column in `dE`, default 2
- `namesX` - names of a possible subset of columns in `dX` to be considered as potential causal regulators (default `names(dX)`)

The optional parameter `method` determines the LLR mixture distribution fitting method and can be either `moments` (default) for the method of moments, or `kde` for kernel-based density estimation.

The optional parameter `FDR` can be used to return only a subset of interactions with a desired expected FDR value (q-value threshold) (default 1.0, no filtering).

The optional parameter `sorted` determines if the output must be sorted by increasing q-value / decreasing posterior probability (`sorted=true`, the default) or by causal factor (column names of `dX2`) (`sorted=false`)
    
See also [`findr(::Matrix,::Array,::Array,::Matrix)`](@ref), [`combineprobs`](@ref), [`stackprobs`](@ref), [`globalfdr!`](@ref).
"""
function findr(dX1::T, dX2::T, dG::T, dE::T; colG=1, colX=2, namesX=[], method="moments", combination="IV", FDR=1.0, sorted=true) where T<:AbstractDataFrame
    if combination == "none"
        error("Returning posterior probabilities for individual tests is not supported with DataFrame inputs. Set combination argument to one of \"IV\", \"mediation\", or \"orig\", or use matrix inputs.")
    elseif combination in Set(["IV","mediation","orig"])
        # Create the array with SNP-Gene pairs
        pairGX = getpairs(dX2, dG, dE; colG = colG, colX = colX, namesX = namesX)
        # Call BioFindr on numeric data
        PP = findr_matrix(Matrix(dX1), Matrix(dX2), Matrix(dG), pairGX; method = method, combination = combination)
        dP = stackprobs(PP, names(dX2)[pairGX[:,2]], names(dX1))
        globalfdr!(dP, FDR = FDR, sorted = sorted)
        return dP
    else
        error("Combination parameter must be one of \"IV\", \"mediation\", or \"orig\"")
    end
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


