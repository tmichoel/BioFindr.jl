module BioFindr

# External packages

using DataFrames
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
export findr, dagfindr!

# Define test names as variables

corr = "corr"
link = "link"
med = "med"
relev = "relev"
pleio = "pleio"

# Load the different code files
include("findr_coexpression.jl");

include("findr_association.jl");

include("findr_causal_inference.jl");

include("findr_pvalues.jl");

include("generate_test_data.jl");

include("supernormalization.jl")

include("lbeta.jl")

include("randomLLR.jl")

include("realLLR.jl")

include("posteriorprobs.jl")

include("bayesiannets.jl")

include("utils.jl")

end