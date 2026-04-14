using BioFindr
using Test
using Distributions
using Statistics
using StatsBase
using SpecialFunctions
using KernelDensity
using Printf
using Random
using ScientificTypes
using DataFrames
using Graphs


@testset "BioFindr.jl" begin
    # Write your tests here.
   
    @testset "Supernormalization" begin
        include("supernormalization_tests.jl");
    end

    @testset "Real LLRs" begin
        include("real_llr_tests.jl");
    end

    @testset "LBeta" begin
        include("lbeta_tests.jl");
    end

    @testset "Random LLRs" begin
        include("random_llr_tests.jl");
    end

    @testset "Posterior probabilities" begin
        include("postprobs_tests.jl");
    end

    @testset "Utils" begin
        include("utils_tests.jl");
    end

    @testset "findr and findr_matrix" begin
        include("findr_tests.jl");
    end

    @testset "dagfindr!" begin
        include("dagfindr_tests.jl");
    end
end
