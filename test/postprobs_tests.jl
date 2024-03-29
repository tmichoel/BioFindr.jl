# Generate supernormalized test data

# We're setting a random seed because the tests should pass for most random data, but occasionally data are generated that fail the LBeta test_moments condition, while this condition must be satisfied to be able to test the fit_mixdist_mom functions. See the [reproducibility section of the Random documentation](https://docs.julialang.org/en/v1/stdlib/Random/#Reproducibility) for justification.
 
Random.seed!(1234)
nA = 1
nB = 1000
fB = 0.2
ns = 100
ng = 3
maf = 0.3
bGA = 1.
bAB = 1.
ϵ = 0.5
δ = 0.5
ρ = 0.1
supernormalize = true
E, Ycol, Y, istarget = BioFindr.generate_test_data(nA, nB, fB, ns, ng, maf, bGA, bAB, supernormalize)

# compute likelihood ratios
llr = BioFindr.realLLR_col(Y,Ycol);
llr1 = BioFindr.realLLR_col(Y,E);
llr2, llr3, llr4, llr5 = BioFindr.realLLR_col(Y,Ycol,E);

@testset "pprob_col test 0" begin
    # test method of moments
    pp_mom = BioFindr.pprob_col(Y,Ycol,method="moments");
    @test all(pp_mom .≈ BioFindr.fit_mixdist_mom(llr,ns)[1])
    # test kde method
    pp_kde = BioFindr.pprob_col(Y,Ycol,method="kde");
    @test all(pp_kde .≈ BioFindr.fit_mixdist_KDE(llr,ns))
end

@testset "pprob_col test 2" begin
    # test method of moments
    pp_mom = BioFindr.pprob_col(Y,E,method="moments");
    @test all(pp_mom .≈ BioFindr.fit_mixdist_mom(llr1,ns,ng,:link)[1])
    # test kde method
    pp_kde = BioFindr.pprob_col(Y,E,method="kde");
    @test all(pp_kde .≈ BioFindr.fit_mixdist_KDE(llr1,ns,ng,:link))
end

@testset "pprob_col test 2-5" begin 
    # test method of moments
    pp_mom = BioFindr.pprob_col(Y,Ycol,E,method="moments");
    @test all(pp_mom[:,1] .≈ BioFindr.fit_mixdist_mom(llr2,ns,ng,:link)[1])
    @test all(pp_mom[:,2] .≈ 1 .- BioFindr.fit_mixdist_mom(llr3,ns,ng,:med)[1])
    @test all(pp_mom[:,3] .≈ BioFindr.fit_mixdist_mom(llr4,ns,ng,:relev)[1])
    @test all(pp_mom[:,4] .≈ BioFindr.fit_mixdist_mom(llr5,ns,ng,:pleio)[1])
    # test kde method
    pp_kde = BioFindr.pprob_col(Y,Ycol,E,method="kde");
    @test all(pp_kde[:,1] .≈ BioFindr.fit_mixdist_KDE(llr2,ns,ng,:link))
    @test all(pp_kde[:,2] .≈ 1 .- BioFindr.fit_mixdist_KDE(llr3,ns,ng,:med))
    @test all(pp_kde[:,3] .≈ BioFindr.fit_mixdist_KDE(llr4,ns,ng,:relev))
    @test all(pp_kde[:,4] .≈ BioFindr.fit_mixdist_KDE(llr5,ns,ng,:pleio))
end

@testset "fit_mixdist_mom test" begin
    # fit mixture distributions to llr values
    pp0, d0 = BioFindr.fit_mixdist_mom(llr,ns);
    pp2, d2 = BioFindr.fit_mixdist_mom(llr2,ns,ng,:link);
    pp3, d3 = BioFindr.fit_mixdist_mom(llr3,ns,ng,:med);
    pp4, d4 = BioFindr.fit_mixdist_mom(llr4,ns,ng,:relev);
    pp5, d5 = BioFindr.fit_mixdist_mom(llr5,ns,ng,:pleio);

    # test if fit_mixdist_mom returns correct null distribution
    @test d0.components[1] == BioFindr.LBeta(1,ns-2)
    @test d2.components[1] == BioFindr.LBeta(ng-1,ns-ng)
    @test d3.components[1] == BioFindr.LBeta(ng-1,ns-ng-1)
    @test d4.components[1] == BioFindr.LBeta(ng,ns-ng-1)
    @test d5.components[1] == BioFindr.LBeta(1,ns-ng-1)

    # test if fit_mixdist_mom returns correct mixture proportions
    @test d0.prior.p[1] ≈ BioFindr.pi0est(BioFindr.nullpval(llr,ns))
    @test d2.prior.p[1] ≈ BioFindr.pi0est(BioFindr.nullpval(llr2,ns,ng,:link))
    @test d3.prior.p[1] ≈ BioFindr.pi0est(BioFindr.nullpval(llr3,ns,ng,:med))
    @test d4.prior.p[1] ≈ BioFindr.pi0est(BioFindr.nullpval(llr4,ns,ng,:relev))
    @test d5.prior.p[1] ≈ BioFindr.pi0est(BioFindr.nullpval(llr5,ns,ng,:pleio))

    # test if fit_mixdist_mom returns alternative distribution with correct limits
    @test d0.components[2].α ≥ d0.components[1].α
    @test d0.components[2].β ≤ d0.components[1].β
    @test d2.components[2].α ≥ d2.components[1].α
    @test d2.components[2].β ≤ d2.components[1].β
    @test d3.components[2].α ≥ d3.components[1].α
    @test d3.components[2].β ≤ d3.components[1].β
    @test d4.components[2].α ≥ d4.components[1].α
    @test d4.components[2].β ≤ d4.components[1].β
    @test d5.components[2].α ≥ d5.components[1].α
    @test d5.components[2].β ≤ d5.components[1].β

    # test if posterior probabilities are between 0 and 1
    @test all(pp0 .≥ 0)
    @test all(pp0 .≤ 1)
    @test all(pp2 .≥ 0)
    @test all(pp2 .≤ 1)
    @test all(pp3 .≥ 0)
    @test all(pp3 .≤ 1)
    @test all(pp4 .≥ 0)
    @test all(pp4 .≤ 1)
    @test all(pp5 .≥ 0)
    @test all(pp5 .≤ 1) 

    # test if posterior probabilities increase (test 0,2,4,5) / decrease (test 3) monotonically with llr values
    @test issorted(round.(pp0,digits=3)[sortperm(llr)])
    @test issorted(round.(pp2,digits=3)[sortperm(llr2)])
    @test issorted(round.(pp3,digits=3)[sortperm(llr3)])
    @test issorted(round.(pp4,digits=3)[sortperm(llr4)])
    @test issorted(round.(pp5,digits=3)[sortperm(llr5)])

    # @test sortperm(llr) == sortperm(pp0)
    # @test sortperm(llr2) == sortperm(pp2)
    # @test sortperm(llr3) == sortperm(pp3)
    # @test sortperm(llr4) == sortperm(pp4)
    # @test sortperm(llr5) == sortperm(pp5)
end 

@testset "fit_mixdist_KDE" begin
    # fit mixture distributions to llr values
    pp0 = BioFindr.fit_mixdist_KDE(llr,ns);
    pp2 = BioFindr.fit_mixdist_KDE(llr2,ns,ng,:link);
    pp3 = BioFindr.fit_mixdist_KDE(llr3,ns,ng,:med);
    pp4 = BioFindr.fit_mixdist_KDE(llr4,ns,ng,:relev);
    pp5 = BioFindr.fit_mixdist_KDE(llr5,ns,ng,:pleio);
    
    # test if posterior probabilities are between 0 and 1
    @test all(pp0 .≥ 0)
    @test all(pp0 .≤ 1)
    @test all(pp2 .≥ 0)
    @test all(pp2 .≤ 1)
    @test all(pp3 .≥ 0)
    @test all(pp3 .≤ 1)
    @test all(pp4 .≥ 0)
    @test all(pp4 .≤ 1)
    @test all(pp5 .≥ 0)
    @test all(pp5 .≤ 1)

    # test if posterior probabilities increase (test 0,2,4,5) / decrease (test 3) monotonically with llr values
    @test issorted(round.(pp0,digits=3)[sortperm(llr)])
    @test issorted(round.(pp2,digits=3)[sortperm(llr2)])
    @test issorted(round.(pp3,digits=3)[sortperm(llr3)])
    @test issorted(round.(pp4,digits=3)[sortperm(llr4)])
    @test issorted(round.(pp5,digits=3)[sortperm(llr5)])
end

@testset "fit_kde tests" begin
    pd = BioFindr.fit_kde(llr)
    z = log.(exp.(2 .* llr) .- 1)
    dfit = kde(z[llr .> 0])
    @test all(pd .≈ 2 * pdf(dfit, z) .* (1 .+ exp.(-z)))
end

@testset "pi0est tests" begin
    # Test if pi0est returns correct estimate for known input
    pval = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    @test BioFindr.pi0est(pval) ≈ 1.0
end

