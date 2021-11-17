using SigmaShift
using Documenter

DocMeta.setdocmeta!(SigmaShift, :DocTestSetup, :(using SigmaShift); recursive=true)

makedocs(;
    modules=[SigmaShift],
    authors="G Jake Gebbie",
    repo="https://github.com/ggebbie/SigmaShift.jl/blob/{commit}{path}#{line}",
    sitename="SigmaShift.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ggebbie.github.io/SigmaShift.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ggebbie/SigmaShift.jl",
    devbranch="main",
)
