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
    PP = findr_matrix(Matrix(dX); cols = cols, method = method, combination = combination)
    dP = stackprobs(PP, colnames, names(dX))
    globalfdr!(dP, FDR = FDR, sorted = sorted)
    return dP
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
    dP = stackprobs(findr_matrix(Matrix(dX), Matrix(dG); method = method), names(dG), names(dX))
    globalfdr!(dP, FDR = FDR, sorted = sorted)
    return dP
end

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
