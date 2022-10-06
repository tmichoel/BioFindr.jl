module Findr

# Write your package code here.

using Statistics
using StatsBase
using Distributions
using SpecialFunctions

include("supernormalization.jl")

include("randomLLR.jl")

include("realLLR.jl")

include("posteriorprobs.jl")

include("bayesiannets.jl")

end
