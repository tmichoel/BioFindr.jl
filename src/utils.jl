"""
    coerce_scitypes!(df, scitype)

Coerce all columns of dataframe `df` to the [scientic type](https://github.com/JuliaAI/ScientificTypes.jl) `scitype`. If `df` contains gene expression data, `scitype` can be  `Continuous` or `Count`. If `df` contains genotype data (or categorical data more generally), `scitype` can be `Multiclass` or `OrderedFactor`. Note though that genotypes are always treated as *unordered* categorical variables in BioFindr, and the ordering of levels for `OrderedFactor` data is not used. Continuous genotypes (e.g. expected allele counts outputted by genotype imputation methods) are not supported and must be converted to integers before calling this function.
"""
function coerce_scitypes!(df, scitype)
    @assert scitype in [ScientificTypes.Continuous, ScientificTypes.Count, ScientificTypes.Multiclass, ScientificTypes.OrderedFactor]
    for col in names(df)
        df[!, col] = coerce(df[!, col], scitype)
    end
end


"""
    test_scitype(df, scitype)

Test if all columns of dataframe `df` have the scientific type `scitype`. If `scitype` is `Continuous` or `Count`, the function will return `true` if all columns are of that type. If `scitype` is `Multiclass` or `OrderedFactor`, the function will return `true` if all columns are of that type, although the columns may have different levels. If the columns have different types, the function will return `false`.
"""
function test_scitype(df, scitype)
    # get the scientific types of the columns of df
    sct = unique(schema(df).scitypes)
    if scitype == ScientificTypes.Continuous
        # all columns should be continuous
        return all(sct .== ScientificTypes.Continuous)
    elseif scitype == ScientificTypes.Count
        # all columns should be count
        return all(sct .== ScientificTypes.Count)
    elseif scitype == ScientificTypes.Multiclass
        # all columns should be multiclass but need not have the same levels
        return all(sct .<: ScientificTypes.Multiclass)
    elseif scitype == ScientificTypes.OrderedFactor
        # all columns should be ordered factor but need not have the same levels
        return all(sct .<: ScientificTypes.OrderedFactor)
    else
        error("scitype must be one of Continuous, Count, Multiclass, or OrderedFactor")
    end
end

"""
    getpairs(dX::T, dG::T, dE::T; colG=1, colX=2)

Get pairs of indices of matching columns from dataframes `dX` and `dG`, with column names that should be matched listed in dataframe `dE`. The optional parameters `colG` (default value 1) and `colX` (default value 2) indicate which columns of `dE` need to be used for matching, either as a column number (integer) or column name (string). The optional parameter `namesX` can be used to match rows in `dE` to only a subset of the column names of `dX`.
"""
function getpairs(dX::T, dG::T, dE::T; colG=1, colX=2, namesX=[]) where T<:AbstractDataFrame
    # if namesX is empty, use all column names of dX
    if isempty(namesX)
        namesX = names(dX)
    end
    
    # Extract dX ID column from dE
    idX = select(dE, colX)[:,1]
    # Extract dG ID column from dE
    idG = select(dE, colG)[:,1]

    # Keep only rows in idX and idG where idX is in namesX
    row_select = findall(x -> x in namesX, idX)
    idX = idX[row_select]
    idG = idG[row_select]

    # Create the array with idG-idX pairs
    pairsGX = zeros(Int64,length(idX),2);
    for rowE = axes(pairsGX,1)
        pairsGX[rowE,1] = findfirst(idG[rowE] .== names(dG))
        pairsGX[rowE,2] = findfirst(idX[rowE] .== names(dX))
    end
    # sort by dX ID to preserve order from gene expression matrix
    pairsGX = pairsGX[sortperm(pairsGX[:,2]),:]
end

"""
    symprobs(P; combination="prod")

Symmetrize a square matrix of posterior probabilities `P`. The optional parameter `combination` defines the symmetrization method:

- `none`: do nothing (default)
- `prod`: ``P'_{ij}=P_{ij}P_{ji}``
- `mean`: ``P'_{ij}=\\frac{1}{2}(P_{ij} + P_{ji})``
- `anti`: ``P'_{ij}=\\frac{1}{2}(P_{ij} + 1 - P_{ji})``

Note that the `anti` option defines "antisymmetric" probabilities, ``P'_{ij} +  P'_{ji} = 1``, where evidence *for* a causal interaction ``i\\to j`` is also considered evidence *against* the opposite interaction ``j\\to i``.
"""
function symprobs(P::Matrix{T}; combination="none") where T<:AbstractFloat
    if size(P,1) != size(P,2) && combination != "none"
        error("Input matrix must be square")
    end
    if combination == "none"
        return P
    elseif combination == "prod"
        return P .* P'
    elseif combination == "mean"
        return 0.5 .* (P .+ P')
    elseif combination == "anti"
        return 0.5 .* (P .+ 1 .- P')
    else
        error("Combination parameter must be one of \"prod\", \"mean\", or \"anti\"")
    end
end

"""
    combineprobs(P; combination="none")

Combine posterior probabilities `P` for multiple likelihood likelihood ratio tests in a single probability (local precision) value.

The optional parameter `combination` defines the combination test:

- `none`: do nothing, return the input `P` (default)
- `mediation`: the mediation test (``P_2 P_3``)
- `IV`: the instrumental variable or non-independence test (``P_2 P_5``)
- `orig`: BioFindr's original combination (``\\frac{1}{2}(P_2 P_5 + P_4)``

The input must be a three-dimensional array where the second dimension has size 4 and indexes the individual BioFindr tests (test 2-5). The output is a matrix of size `size(P,1) x size(P,3)`.
"""
function combineprobs(P; combination="none")
    if combination == "none"
        return P
    elseif combination == "IV"
        # P2 x P5
        return P[:,1,:] .* P[:,4,:]
    elseif combination == "mediation"
        # P2 x P3
        return P[:,1,:] .* P[:,2,:]
    elseif combination == "orig"
        # 0.5 (P2 x P5 + P4)
        return 0.5 .*( P[:,1,:] .* P[:,4,:] .+ P[:,3,:] )
    else
        error("Combination parameter must be one of \"none\", \"IV\", \"mediation\", or \"orig\"")
    end
end

"""
    stackprobs(P,colnames,rownames;nodiag=true)

Convert a matrix of pairwise posterior probabilities `P` with column and row names `colnames` and `rownames`, respectively, to a stacked dataframe with `Source`, `Target`, and `Probability` columns, corresponding respectively to a column name, a row name, and the value of `P` in the corresponding row and column pair.

The optional parameter `nodiag` determines if self-interactions (equal row and column name) are excluded (`nodiag=true`, default) or not (`nodiag=false`).
"""
function stackprobs(P,colnames,rownames;nodiag=true)
    # First put the matrix of probabilities in a dataframe
    dP = DataFrame(P, colnames)
    # Add column with row names
    insertcols!(dP, 1, "Target" => rownames)
    dP = stack(dP, Not(:Target), variable_name=:Source, value_name=:"Probability")
    if nodiag
        # remove rows where source and target are the same
        filter!(row -> row.Target != row.Source, dP)
    end
    return dP[!, [2,1,3]]
end

"""
    globalfdr(P::Array{T},FDR) where T<:AbstractFloat

For an array (matrix or vector) `P` of posterior probabilities (local precision values), compute their corresponding q-values `Q`, and return the indices of `P` with q-value less than a desired global false discovery rate `FDR`.

See also [`qvalue`](@ref)
"""
function globalfdr(P::Array{T},FDR) where T<:AbstractFloat
    Qvec = qvalue(vec(P))
    # return entries with q-value <= FDR
    if isa(P,Vector)  
        return findall(Qvec .<= FDR), Qvec
    else
        # reshape q-values in original shape
        Q = reshape(Qvec,size(P))
        return findall(Q .<= FDR), Q
    end
end

"""
    globalfdr!(dP::T; FDR=1.0, sorted=true) where T<:AbstractDataFrame

For a DataFrame `dP` of posterior probabilities (local precision values), compute their corresponding q-values and keep only the rows with q-value less than a desired global false discovery rate `FDR` (default value 1, no selection). `dP` is assumed to be the output of a `findr` run with columns `Source`, `Target`, and `Probability`. The output DataFrame mirrors the structure of `dP`, keeping only the selected rows, and with an additional column `qvalue`. The output is sorted by `qvalue` if the optional argument `sorted` is `true` (default). If `dP` already contains a column `qvalue`, only the filtering and optional sorting are performed. 
"""
function globalfdr!(dP::T; FDR=1.0, sorted=true) where T<:AbstractDataFrame
    # test if dP already has a q-value column, this allows repeated calling of the function for additional filtering or sorting
    if ∉("qvalue",names(dP))
        qval = qvalue(dP."Probability")
        insertcols!(dP,"qvalue" => qval)
    end
    subset!(dP, :"qvalue" => x -> x .<= FDR)
    if sorted
        sort!(dP, :"qvalue")
    end
end

"""
    qvalue(P::Vector{T}) where T<:AbstractFloat

Convert a vector `P` of posterior probabilities (local precisions) to a vector of q-values. For a threshold value `c` on the posterior probabilities `P`, the global FDR, ``FDR(c)`` is defined as one minus the average local precision:

``FDR(c) = 1 - \\frac{1}{N_c} \\sum_{i\\colon P_i\\leq c} P_i,``

where ``N_c=\\sharp\\{i\\colon P_i\\leq c\\}`` is the number of selected pairs. The q-value of a given index in `P` is then defined as the smallest FDR at which this pair is still called significant.
"""
function qvalue(P::Vector{T}) where T<:AbstractFloat
    # permutation to sort Pvec in descending order
    I = sortperm(P,rev=true)
    # the inverse permutation
    Iinv = invperm(I)
    # accumulate 1 - mean(Pvec[I])
    qval = 1 .- (cumsum(P[I])./(1:length(P)))
    # q-values must be ordered, if not, set qval[k] = minimum(qval[k:end]) using efficient operations
    if !issorted(qval)
        @info "Average local precisions needed sorting"
        reverse!(qval)
        accumulate!(min,qval,qval)
        reverse!(qval)
        # for k = eachindex(qval)
        #     qval[k] = minimum(qval[k:end])
        # end
    end
    @assert issorted(qval)
    # return to original order
    qval = qval[Iinv]
end