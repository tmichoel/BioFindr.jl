@testset "pi0est tests" begin
    # Test if pi0est returns correct estimate for known input
    pval = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    @test Findr.pi0est(pval) ≈ 1.0
end

@testset "fit_kde tests" begin
    # Test if fit_kde returns correct pdf for known input
    llr = log.(1:10)
    pd = Findr.fit_kde(llr)
    z = log.(exp.(2 .* llr) .- 1)
    dfit = kde(z[llr .> 0])
    @test all(pd .≈ 2 * pdf(dfit, z) .* (1 .+ exp.(-z)))
end