using Findr
using Documenter

DocMeta.setdocmeta!(Findr, :DocTestSetup, :(using Findr); recursive=true)

makedocs(;
    modules=[Findr],
    authors="tmichoel <11647967+tmichoel@users.noreply.github.com> and contributors",
    repo="https://github.com/tmichoel/Findr.jl/blob/{commit}{path}#{line}",
    sitename="Findr.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tmichoel.github.io/Findr.jl",
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
    repo="github.com/tmichoel/Findr.jl",
    devbranch="main",
)
