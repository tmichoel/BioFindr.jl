"""
    getpairs(dX::T, dG::T, dE::T; colG=1, colX=2)

Get pairs of indices of matching columns from dataframes `dX` and `dG`, with column names that should be matched listed in dataframe `dE`. The optional parameters `colG` (default value 1) and `colX` (default value 2) indicate which columns of `dE` need to be used for matching, either as a column number (integer) or column name (string).
"""
function getpairs(dX::T, dG::T, dE::T; colG=1, colX=2) where T<:AbstractDataFrame
    # Extract dX ID column from dE
    if typeof(colX) <: Int
        idX = dE[!,colX]
    else
        idX = dE.colX
    end
    # Extract dG ID column from dE
    if typeof(colG) <: Int
        idG = dE[!,colG]
    else
        idG = dE.colG
    end
    # Create the array with idG-idX pairs
    pairsGX = zeros(Int64,nrow(dE),2);
    for rowE = axes(pairsGX,1)
        pairsGX[rowE,1] = findfirst(idG[rowE] .== names(dG))
        pairsGX[rowE,2] = findfirst(idX[rowE] .== names(dX))
    end
    # sort by dX ID to preserve order from gene expression matrix
    pairsGX = pairsGX[sortperm(pairsGX[:,2]),:]
end

"""
    combineprobs(P; combination="none")

Combine posterior probabilities `P` for multiple likelihood likelihood ratio tests in a single probability (local precision) value.
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

Convert a matrix of pairwise posterior probabilities `P` with column and row names `colnames` and `rownames`, respectively, to a stacked dataframe with `Source`, `Target`, and `Posterior probability` columns, corresponding respectively to a column name, a row name, and the value of `P` in the corresponding row and column pair.
"""
function stackprobs(P,colnames,rownames;nodiag=true)
    # First put the matrix of probabilities in a dataframe
    dP = DataFrame(P, colnames)
    # Add column with row names
    insertcols!(dP, 1, "Target" => rownames)
    dP = stack(dP, Not(:Target), variable_name=:Source, value_name=:"Posterior probability")
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
    globalfdr!(dP::T;FDR=1.0, sorted=true) where T<:AbstractDataFrame

For a DataFrame `dP` of posterior probabilities (local precision values), compute their corresponding q-values `Q` and keep only the rows with `Q` less than a desired global false discovery rate `FDR`. `dP` is assumed to be the output of a `findr` run with columns `Source`, `Target`, and `Posterior probability`. The output DataFrame mirrors the structure of `dP`, keeping only the selected rows, and with an additional column `q-value`. The output is sorted by `q-value` if the optional argument `sorted` is `true` (default).
"""
function globalfdr!(dP::T; FDR=1.0, sorted=true) where T<:AbstractDataFrame
    # test if dP already has a q-value column, this allows repeated calling of the function for additional filtering or sorting
    if ∉("q-value",names(dP))
        qval = qvalue(dP."Posterior probability")
        insertcols!(dP,"q-value" => qval)
    end
    subset!(dP, :"q-value" => x -> x .<= FDR)
    if sorted
        sort!(dP, :"q-value")
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
    # q-values must be ordered
    if !issorted(qval)
        error("sorting")
        for k = eachindex(qval)
            qval[k] = minimum(qval[k:end])
        end
    end
    # return to original order
    qval[Iinv]
end