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
    combineprobs(PP; combination="none")
"""
function combineprobs(PP; combination="none")
    if combination == "none"
        return PP
    elseif combination == "IV"
        # P2 x P5
        return PP[:,1,:] .* PP[:,4,:]
    elseif combination == "mediation"
        # P2 x P3
        return PP[:,1,:] .* PP[:,2,:]
    elseif combination == "orig"
        # 0.5 (P2 x P5 + P4)
        return 0.5 .*( PP[:,1,:] .* PP[:,4,:] .+ PP[:,3,:] )
    else
        error("Combination parameter must be one of \"none\", \"IV\", \"mediation\", or \"orig\"")
    end
end