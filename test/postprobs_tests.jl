# Generate supernormalized test data
nA = 1
nB = 1000
fB = 0.2
ns = 100
ng = 3
maf = 0.3
bGA = 2.
bAB = 2.
supernormalize = true
E, Ycol, Y, istarget = Findr.generate_test_data(nA, nB, fB, ns, ng, maf, bGA, bAB, supernormalize)

@testset "pprob_col test 0" begin
    # compute likelihood ratios
    llr = Findr.realLLR_col(Y,Ycol);
    # test method of moments
    pp_mom = Findr.pprob_col(Y,Ycol,method="moments");
    @test all(pp_mom .≈ Findr.fit_mixdist_mom(llr,ns)[1])
    # test kde method
    pp_kde = Findr.pprob_col(Y,Ycol,method="kde");
    @test all(pp_kde .≈ Findr.fit_mixdist_KDE(llr,ns))
end

@testset "pprob_col test 2" begin
    # compute likelihood ratios
    llr2 = Findr.realLLR_col(Y,E);
    # test method of moments
    pp_mom = Findr.pprob_col(Y,E,method="moments");
    @test all(pp_mom .≈ Findr.fit_mixdist_mom(llr2,ns,ng,:link)[1])
    # test kde method
    pp_kde = Findr.pprob_col(Y,E,method="kde");
    @test all(pp_kde .≈ Findr.fit_mixdist_KDE(llr2,ns,ng,:link))
end

@testset "pprob_col test 2-5" begin
    # compute likelihood ratios
    llr2, llr3, llr4, llr5 = Findr.realLLR_col(Y,Ycol,E);
    # test method of moments
    pp_mom = Findr.pprob_col(Y,Ycol,E,method="moments");
    @test all(pp_mom[:,1] .≈ Findr.fit_mixdist_mom(llr2,ns,ng,:link)[1])
    @test all(pp_mom[:,2] .≈ Findr.fit_mixdist_mom(llr3,ns,ng,:med)[1])
    @test all(pp_mom[:,3] .≈ Findr.fit_mixdist_mom(llr4,ns,ng,:relev)[1])
    @test all(pp_mom[:,4] .≈ Findr.fit_mixdist_mom(llr5,ns,ng,:pleio)[1])
    # test kde method
    pp_kde = Findr.pprob_col(Y,Ycol,E,method="kde");
    @test all(pp_kde[:,1] .≈ Findr.fit_mixdist_KDE(llr2,ns,ng,:link))
    @test all(pp_kde[:,2] .≈ Findr.fit_mixdist_KDE(llr3,ns,ng,:med))
    @test all(pp_kde[:,3] .≈ Findr.fit_mixdist_KDE(llr4,ns,ng,:relev))
    @test all(pp_kde[:,4] .≈ Findr.fit_mixdist_KDE(llr5,ns,ng,:pleio))
end



@testset "pi0est tests" begin
    # Test if pi0est returns correct estimate for known input
    pval = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    @test Findr.pi0est(pval) ≈ 1.0
end

# @testset "fit_kde tests" begin
#     # Test if fit_kde returns correct pdf for known input
#     llr = log.(1:10)
#     pd = Findr.fit_kde(llr)
#     z = log.(exp.(2 .* llr) .- 1)
#     dfit = kde(z[llr .> 0])
#     @test all(pd .≈ 2 * pdf(dfit, z) .* (1 .+ exp.(-z)))
# end