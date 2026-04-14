@testset "nulldist tests" begin
    ns = 100
    ng = 10

    # Test if nulldist returns correct distributions
    @test BioFindr.nulldist(ns, ng, :corr) == BioFindr.LBeta(1, ns-2)
    @test BioFindr.nulldist(ns, ng, :link) == BioFindr.LBeta(ng-1, ns-ng)
    @test BioFindr.nulldist(ns, ng, :med) == BioFindr.LBeta(ng-1, ns-ng-1)
    @test BioFindr.nulldist(ns, ng, :relev) == BioFindr.LBeta(ng, ns-ng-1)
    @test BioFindr.nulldist(ns, ng, :pleio) == BioFindr.LBeta(1, ns-ng-1)

    # Test if nulldist throws an error for invalid test argument
    @test_throws ErrorException BioFindr.nulldist(ns, ng, :invalid)
end

# Test nullpval function
@testset "nullpval tests" begin
    ns = 100
    ng = 10
    llr = exp.(rand(ns))

    # Test if nullpval returns correct p-values
    @test all(BioFindr.nullpval(llr, ns, ng, :corr) .== ccdf.(BioFindr.LBeta(1, ns-2), llr))
    @test all(BioFindr.nullpval(llr, ns, ng, :link) .== ccdf.(BioFindr.LBeta(ng-1, ns-ng), llr))
    @test all(BioFindr.nullpval(llr, ns, ng, :med) .== ccdf.(BioFindr.LBeta(ng-1, ns-ng-1), llr))
    @test all(BioFindr.nullpval(llr, ns, ng, :relev) .== ccdf.(BioFindr.LBeta(ng, ns-ng-1), llr))
    @test all(BioFindr.nullpval(llr, ns, ng, :pleio) .== ccdf.(BioFindr.LBeta(1, ns-ng-1), llr))

    # Test if nullpval throws an error for invalid test argument
    @test_throws ErrorException BioFindr.nullpval(llr, ns, ng, :invalid)
end

# Test nulllog10pval function
@testset "nulllog10pval tests" begin
    ns = 100
    ng = 10
    llr = exp.(rand(ns))

    # Test if nullpval returns correct p-values
    @test all(BioFindr.nulllog10pval(llr, ns, ng, :corr) .== -logccdf.(BioFindr.LBeta(1, ns-2), llr)/log(10)) 
    @test all(BioFindr.nulllog10pval(llr, ns, ng, :link) .== -logccdf.(BioFindr.LBeta(ng-1, ns-ng), llr)/log(10)) 
    @test all(BioFindr.nulllog10pval(llr, ns, ng, :med) .== -logccdf.(BioFindr.LBeta(ng-1, ns-ng-1), llr)/log(10))
    @test all(BioFindr.nulllog10pval(llr, ns, ng, :relev) .== -logccdf.(BioFindr.LBeta(ng, ns-ng-1), llr)/log(10)) 
    @test all(BioFindr.nulllog10pval(llr, ns, ng, :pleio) .== -logccdf.(BioFindr.LBeta(1, ns-ng-1), llr)/log(10))

    # Test if nullpval throws an error for invalid test argument
    @test_throws ErrorException BioFindr.nulllog10pval(llr, ns, ng, :invalid)
end

# Test nullpdf function
@testset "nullpdf tests" begin
    ns = 100
    ng = 10
    llr = exp.(rand(ns))

    # Test if nullpval returns correct p-values
    @test all(BioFindr.nullpdf(llr, ns, ng, :corr) .== BioFindr.pdf.(BioFindr.LBeta(1, ns-2), llr))
    @test all(BioFindr.nullpdf(llr, ns, ng, :link) .== BioFindr.pdf.(BioFindr.LBeta(ng-1, ns-ng), llr))
    @test all(BioFindr.nullpdf(llr, ns, ng, :med) .== BioFindr.pdf.(BioFindr.LBeta(ng-1, ns-ng-1), llr))
    @test all(BioFindr.nullpdf(llr, ns, ng, :relev) .== BioFindr.pdf.(BioFindr.LBeta(ng, ns-ng-1), llr))
    @test all(BioFindr.nullpdf(llr, ns, ng, :pleio) .== BioFindr.pdf.(BioFindr.LBeta(1, ns-ng-1), llr))

    # Test if nullpval throws an error for invalid test argument
    @test_throws ErrorException BioFindr.nullpdf(llr, ns, ng, :invalid)
end