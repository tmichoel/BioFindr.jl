module Findr

# Write your package code here.

using Statistics
using StatsBase
using Distributions
using SpecialFunctions
using KernelDensity
using LinearAlgebra
using InvertedIndices

# Define test names as variables
corr = "corr"
link = "link"
med = "med"
relev = "relev"
pleio = "pleio"

# Load the different code files

include("supernormalization.jl")

include("randomLLR.jl")

include("realLLR.jl")

include("posteriorprobs.jl")

include("bayesiannets.jl")

end
