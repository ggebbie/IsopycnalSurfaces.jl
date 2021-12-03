using IsopycnalSurfaces
using Documenter

DocMeta.setdocmeta!(IsopycnalSurfaces, :DocTestSetup, :(using IsopycnalSurfaces); recursive=true)

makedocs(;
    modules=[IsopycnalSurfaces],
    authors="G Jake Gebbie",
    repo="https://github.com/ggebbie/IsopycnalSurfaces.jl/blob/{commit}{path}#{line}",
    sitename="IsopycnalSurfaces.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ggebbie.github.io/IsopycnalSurfaces.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ggebbie/IsopycnalSurfaces.jl",
    devbranch="main",
)
