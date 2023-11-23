# parameters for LBeta distributions, integer and float
α₀  = 1
β₀  = 100
α₁  = 2.5
β₁  = 101.5

# create LBeta distributions with parameters of different types
d₀  = Findr.LBeta(α₀,β₀)
d₁  = Findr.LBeta(α₁,β₁)
d₂  = Findr.LBeta(α₀,β₁)
d₃  = Findr.LBeta(α₁,β₀)

# test that the parameters are correct
@testset "Parameter values" begin
    @test Findr.params(d₀) == (α₀,β₀)
    @test Findr.params(d₁) == (α₁,β₁)
    @test Findr.params(d₂) == (α₀,β₁)
    @test Findr.params(d₃) == (α₁,β₀)
end


# LBeta distribution parameters should always be float regardless of input type
@testset "Parameter type" begin
    @test typeof(Findr.params(d₀)) == Tuple{Float64, Float64}
    @test typeof(Findr.params(d₁)) == Tuple{Float64, Float64}
    @test typeof(Findr.params(d₂)) == Tuple{Float64, Float64}
    @test typeof(Findr.params(d₃)) == Tuple{Float64, Float64}
end

# Test that the pdf is correct
@testset "Pdf range and values" begin
    # LBeta distribution has support on x >= 0
    @test Findr.pdf(d₀, -1) == 0
    # With parameter α<2, the pdf diverges at x=0
    @test Findr.pdf(d₀, 0) == Inf
    # With parameter α>2, the pdf is zero at x=0
    @test Findr.pdf(d₁, 0) == 0
    # With parameter α=2, the pdf is 2/Beta(1,0.5β)
    @test Findr.pdf(Findr.LBeta(2,β₁), 0) == 2/beta(1,0.5β₁)
end