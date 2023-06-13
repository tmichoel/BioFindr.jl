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
    globalfdr(P::Matrix{T},FDR) where T<:AbstractFloat

For a matrix `P` of posterior probabilities (local precision values), find the threshold corresponding to global false discovery rate `FDR`. Return the selected pairs, their posterior probabilities, and their q-values. 

For a threshold value `c`, the global FDR, ``FDR(c)`` is defined as one minus the average local precision:

``FDR(c) = 1 - \frac{1}{N_c} \sum_{i\colon P_i\leq c} P_i,``

where ``N_c=\sharp\{i\colon P_i\leq c\}`` is the number of selected pairs. The q-value of a given pair is defined as the smallest FDR at which this pair is still selected.
"""

function globalfdr(P::Matrix{T},FDR) where T<:AbstractFloat
    
end