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
    globalfdr(P::Array{T},FDR) where T<:AbstractFloat

For an array (matrix or vector) `P` of posterior probabilities (local precision values), compute their corresponding q-values `Q`, and return the indices of `P` with q-value less than a desired global false discovery rate `FDR`. 

For a threshold value `c` on the posterior probabilities `P`, the global FDR, ``FDR(c)`` is defined as one minus the average local precision:

``FDR(c) = 1 - \\frac{1}{N_c} \\sum_{i\\colon P_i\\leq c} P_i,``

where ``N_c=\\sharp\\{i\\colon P_i\\leq c\\}`` is the number of selected pairs. The q-value of a given index in `P` is defined as the smallest FDR at which this pair is still called significant.
"""
function globalfdr(P::Array{T},FDR) where T<:AbstractFloat
    # operate on vector of probabilities
    Pvec = vec(P)
    # permutation to sort Pvec in descending order
    I = sortperm(Pvec,rev=true)
    # the inverse permutation
    Iinv = invperm(I)
    # accumulate 1 - mean(Pvec[I])
    Qvec = 1 .- (cumsum(Pvec[I])./(1:length(Pvec)))
    # q-values must be sorted
    if !issorted(Qvec)
        error("sorting")
        for k = eachindex(Qvec)
            Qvec[k] = minimum(Qvec[k:end])
        end
    end
    # return to original order
    Qvec = Qvec[Iinv]
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
    globalfdr(dP::T,FDR) where T<:AbstractDataFrame

Wrapper for `globalfdr(Matrix(dP),FDR)` when the input `dP` is a DataFrame. `dP` is assumed to be the output of a `findr` where column names correspond to source variables, rows to target variables, and the first column contains the target variable names.

With a DataFrame input, only the selected pairs and their posterior probabilities and q-values are returned, in the form of a DataFrame.
"""
function globalfdr(dP::T,FDR) where T<:AbstractDataFrame
    # Extract matrix of posterior probability values
    P = Matrix(select(dP,Not(1)))
    # call globalfdr on array input
    I, Q = globalfdr(P,FDR)
    # create new dataframe with selected pairs
    df = DataFrame(P = P[I], q = Q[I])
end