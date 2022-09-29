"""
    invnorm(X[, c])

Convert each row of matrix `X` of reals into standard normally distributed values using a rank-based inverse normal transformation.

The formula and default value for the paramater `c` come from this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2921808/
"""
function invnorm(X::Matrix{AbstractFloat}, c=0.375)
    nd = Normal()
    n = size(X,2)
    Y = zeros(size(X))
    for i = axes(X)[1]
        Y[i,:] = map(x -> quantile(nd,(x-c)/(n-2c+1)), tiedrank(X[i,:]))     
    end
    return Y
end