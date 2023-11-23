
# Test groupmeans function
@testset "groupmeans tests" begin
    Y = [1 2 3; 4 5 6; 7 8 9; 10 11 12]
    E = [1, 2, 1, 2]

    # Test groupmeans with 2 arguments
    gs, μ = Findr.groupmeans(Y, E)

    # Test if groupmeans returns correct group sizes
    @test gs == [2, 2]

    # Test if groupmeans returns correct group means
    @test μ == [4 7; 5 8; 6 9]

    # Test groupmeans with 3 arguments
    Ycol = [1, 2, 3, 4]
    gs, μ, μcol = Findr.groupmeans(Y, Ycol, E)

    # Test if groupmeans returns correct group sizes
    @test gs == [2, 2]

    # Test if groupmeans returns correct group means for Y
    @test μ == [4 7; 5 8; 6 9]

    # Test if groupmeans returns correct group means for Ycol
    @test μcol == [2, 3]
end

Y = Findr.supernormalize([1 2 3 4; 4 3 2 1; 1 3 2 4]')
Ycol = Findr.supernormalize([4, 2, 1, 3])
E = [1, 2, 1, 2]

# Test llrstats_col function
@testset "llrstats_col tests" begin
    
    # Test llrstats_col with 2 arguments
    stats = Findr.llrstats_col(Y, E)

    # Test if llrstats_col returns a vector of the correct length
    @test length(stats) == size(Y, 2)

    # Test if llrstats_col returns correct statistics
    @test all(stats .== 1 .- sum((Findr.groupmeans(Y, E)[2]).^2, pweights([2, 2]/4), dims=2))

    # Test llrstats_col with 3 arguments
    ρ, σ, σcol = Findr.llrstats_col(Y, Ycol, E)

    # Test if llrstats_col returns a vector of the correct length
    @test length(ρ) == size(Y, 2)

    # Test if llrstats_col returns correct ρ
    @test all(ρ .== vec(cov(Y, Ycol, corrected=false)))

    # Test if llrstats_col returns correct σ
    gs, μ, μcol = Findr.groupmeans(Y, Ycol, E)
    w = pweights(gs/length(E))
    @test all(σ .== [vec(1 .- sum(μ.^2,w,dims=2)) vec(1 .- sum(μ.*μcol',w,dims=2))])

    # Test if llrstats_col returns correct σcol
    @test σcol == 1 - sum(μcol.^2,w)
end

# Test realLLR_col function
@testset "realLLR_col tests" begin
    # Test realLLR_col with 2 continuous arguments
    LLR = Findr.realLLR_col(Y, Ycol)

    # Test if realLLR_col returns a vector of the correct length
    @test length(LLR) == size(Y, 2)

    # Test if realLLR_col returns correct LLR
    ρ = vec(cov(Y, Ycol, corrected=false))
    @test all(LLR .== -0.5*log.(abs.(1 .- ρ.^2)))

    
    # Test realLLR_col with 1 continuous, 1 discrete argument
    LLR = Findr.realLLR_col(Y, E)

    # Test if realLLR_col returns a vector of the correct length
    @test length(LLR) == size(Y, 2)

    # Test if realLLR_col returns correct LLR
    σ = Findr.llrstats_col(Y, E)
    @test all(LLR .== -0.5*log.(σ))

    # Test realLLR_col with 3 arguments
    llr2, llr3, llr4, llr5 = Findr.realLLR_col(Y, Ycol, E)

    # Test if realLLR_col returns vectors of the correct length
    @test length(llr2) == size(Y, 2)
    @test length(llr3) == size(Y, 2)
    @test length(llr4) == size(Y, 2)
    @test length(llr5) == size(Y, 2)

    # Test if realLLR_col returns correct LLR values
    ρ, σ, σcol = Findr.llrstats_col(Y, Ycol, E)
    @test all(llr4 .== -0.5*log.(abs.(σcol.*σ[:,1] .- (ρ .+ σ[:,2] .- 1).^2)) .+ 0.5*log.(σcol))
    @test all(llr3 .== llr4 .+ 0.5*log.(abs.(1 .- ρ.^2)))
    @test all(llr5 .== llr4 .+ 0.5*log.(σ[:,1]))
end

