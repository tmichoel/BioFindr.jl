# number of samples used for testing
ns = 100;

# supernormalization is rank-based, create a vector of ranks
X₁ = vec(1:ns);

# supernormalize X\_1
Y₁ = Findr.supernormalize(X₁);

# supernormalized data must have mean 0, std 1
@test mean(Y₁) ≈ 0 atol=1e-10
@test std(Y₁,corrected=false) ≈ 1 atol=1e-10 

# supernormalization must preserve ranking
@test issorted(Y₁)

# any random vector of the same size as X\_1 must have the same supernormalization after sorting
@test sort(Findr.supernormalize(rand(ns))) ≈ Y₁ atol=1e-10

# supernormalization must preserve the shape and type of the original data
X₂ = rand(ns,10);
Y₂ = Findr.supernormalize(X₂);
@test size(Y₂) == size(X₂)
@test typeof(Y₂) == typeof(X₂)