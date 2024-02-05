# parameters for LBeta distributions, integer and float
α₀  = 1
β₀  = 100
α₁  = 2.5
β₁  = 101.5

# create LBeta distributions with parameters of different types
d₀  = BioFindr.LBeta(α₀,β₀)
d₁  = BioFindr.LBeta(α₁,β₁)
d₂  = BioFindr.LBeta(α₀,β₁)
d₃  = BioFindr.LBeta(α₁,β₀)

# test that the parameters are correct
@testset "Parameter values" begin
    @test BioFindr.params(d₀) == (α₀,β₀)
    @test BioFindr.params(d₁) == (α₁,β₁)
    @test BioFindr.params(d₂) == (α₀,β₁)
    @test BioFindr.params(d₃) == (α₁,β₀)
end


# LBeta distribution parameters should always be float regardless of input type
@testset "Parameter type" begin
    @test typeof(BioFindr.params(d₀)) == Tuple{Float64, Float64}
    @test typeof(BioFindr.params(d₁)) == Tuple{Float64, Float64}
    @test typeof(BioFindr.params(d₂)) == Tuple{Float64, Float64}
    @test typeof(BioFindr.params(d₃)) == Tuple{Float64, Float64}
end

# Test that the pdf is correct
@testset "Pdf range and values" begin
    # LBeta distribution has support on x >= 0
    @test BioFindr.pdf(d₀, -1) == 0
    # With parameter α<2, the pdf diverges at x=0
    @test BioFindr.pdf(d₀, 0) == Inf
    # With parameter α>2, the pdf is zero at x=0
    @test BioFindr.pdf(d₁, 0) == 0
    # With parameter α=2, the pdf is 2/Beta(1,0.5β)
    @test BioFindr.pdf(BioFindr.LBeta(2,β₁), 0) == 2/beta(1,0.5β₁)
end