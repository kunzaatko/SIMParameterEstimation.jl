using SIMParameterEstimation
using Documenter

DocMeta.setdocmeta!(SIMParameterEstimation, :DocTestSetup, :(
        include(joinpath(@__DIR__, "../test/doctestsetup.jl"))
    ); recursive=true)

makedocs(;
    modules=[SIMParameterEstimation],
    authors="Martin Kunz <martinkunz@email.cz> and contributors",
    repo=Remotes.GitHub("kunzaatko","SIMParameterEstimation.jl"),
    sitename="SIMParameterEstimation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kunzaatko.github.io/SIMParameterEstimation.jl",
        edit_link="trunk",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Utilities" => [
            "Cross Correlation" => "cross-correlation.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/kunzaatko/SIMParameterEstimation.jl",
    devbranch="trunk",
)
