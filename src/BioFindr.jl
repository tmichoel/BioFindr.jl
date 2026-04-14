module BioFindr

# External packages

using DataFrames
using CategoricalArrays
using Statistics
using StatsBase
using Distributions
using Random
using SpecialFunctions
using KernelDensity
using LinearAlgebra
using InvertedIndices
using Graphs
using MetaGraphsNext
using ScientificTypes

# Export statements
export findr, dagfindr!, LBeta, generate_test_data, coerce_scitypes!

# Load the different code files
include("findr.jl");

include("findr_matrix.jl");

include("findr_pvalues.jl");

include("generate_test_data.jl");

include("supernormalization.jl")

include("lbeta.jl")

include("random_llr.jl")

include("real_llr.jl")

include("posteriorprobs.jl")

include("bayesiannets.jl")

include("utils.jl")

end