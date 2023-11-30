@testset "nulldist tests" begin
    ns = 100
    ng = 10

    # Test if nulldist returns correct distributions
    @test Findr.nulldist(ns, ng, :corr) == Findr.LBeta(1, ns-2)
    @test Findr.nulldist(ns, ng, :link) == Findr.LBeta(ng-1, ns-ng)
    @test Findr.nulldist(ns, ng, :med) == Findr.LBeta(ng-1, ns-ng-1)
    @test Findr.nulldist(ns, ng, :relev) == Findr.LBeta(ng, ns-ng-1)
    @test Findr.nulldist(ns, ng, :pleio) == Findr.LBeta(1, ns-ng-1)

    # Test if nulldist throws an error for invalid test argument
    @test_throws ErrorException Findr.nulldist(ns, ng, :invalid)
end

# Test nullpval function
@testset "nullpval tests" begin
    ns = 100
    ng = 10
    llr = exp.(rand(ns))

    # Test if nullpval returns correct p-values
    @test all(Findr.nullpval(llr, ns, ng, :corr) .== ccdf.(Findr.LBeta(1, ns-2), llr))
    @test all(Findr.nullpval(llr, ns, ng, :link) .== ccdf.(Findr.LBeta(ng-1, ns-ng), llr))
    @test all(Findr.nullpval(llr, ns, ng, :med) .== ccdf.(Findr.LBeta(ng-1, ns-ng-1), llr))
    @test all(Findr.nullpval(llr, ns, ng, :relev) .== ccdf.(Findr.LBeta(ng, ns-ng-1), llr))
    @test all(Findr.nullpval(llr, ns, ng, :pleio) .== ccdf.(Findr.LBeta(1, ns-ng-1), llr))

    # Test if nullpval throws an error for invalid test argument
    @test_throws ErrorException Findr.nullpval(llr, ns, ng, :invalid)
end