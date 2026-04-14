push!(LOAD_PATH,"../src/")

using BioFindr
using Documenter

DocMeta.setdocmeta!(BioFindr, :DocTestSetup, :(using BioFindr); recursive=true)

makedocs(;
    modules=[BioFindr],
    authors="tmichoel <11647967+tmichoel@users.noreply.github.com> and contributors",
    repo="https://github.com/tmichoel/BioFindr.jl/blob/{commit}{path}#{line}",
    sitename="BioFindr",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tmichoel.github.io/BioFindr.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "General inference algorithm" => "inference.md",
        "Likelihood ratio tests" => "realLLR.md",
        "Null distributions of the LLRs" => "randomLLR.md",
        "Bayesian inference of posterior probabilities" => "posteriorprobs.md",
        "Tests to evaluate" => "testLLR.md",
        "Bayesian networks" => "bayesiannets.md",
        "Utilities" => "utils.md",
        "List of functions" => "listfunctions.md"       
    ],
)

deploydocs(;
    repo="github.com/tmichoel/BioFindr.jl",
    devbranch="main",
)
