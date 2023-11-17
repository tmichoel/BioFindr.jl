"""
    supernormalize(X[, c])

Convert each column of matrix or dataframe `X` of reals into standard normally distributed values using a rank-based inverse normal transformation. Then scale each column to have variance one.

Note that after the inverse normal transformation, each column has mean zero and identical variance (if we use ordinal ranking). Hence rescaling can be done once on the whole matrix.

The formula and default value for the paramater `c` come from [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2921808/)
"""
function supernormalize(X, c=0.375)
    nd = Normal()
    n = size(X,1)
    Y = zeros(size(X))
    Threads.@threads for i = axes(X,2)
        Y[:,i] = map(x -> quantile(nd,(x-c)/(n-2c+1)), ordinalrank(X[:,i]))
    end
    σ = std(Y[:,1];corrected=false)
    Y = Y/σ
    return Y
end