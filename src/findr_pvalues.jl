#####################################################################
#   Methods to extract p-values instead of posterior probabilities  #
#####################################################################


function findrpval_matrix(X::Matrix{T}) where T<:AbstractFloat
    # Inverse-normal transformation and standardization for each columns of X
    Y = supernormalize(X)
    ns = size(Y,1)
    # Matrix to store posterior probabilities
    ncols = size(Y,2)
    pv = ones(ncols,ncols) # this sets the diagonal elements to one
    # Compute LLRs and p-values for each column separately
    Threads.@threads for col = axes(Y,2)
        llr = realLLR_col(Y[:,Not(col)],Y[:,col])
        pv[Not(col),col] = nulllog10pval(llr,ns)
    end
    return pv
end


