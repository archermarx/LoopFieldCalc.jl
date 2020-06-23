using LoopFieldCalc
using Documenter

makedocs(;
    modules=[LoopFieldCalc],
    authors="Thomas Marks <marksta@umich.edu> and contributors",
    repo="https://github.com/archermarx/LoopFieldCalc.jl/blob/{commit}{path}#L{line}",
    sitename="LoopFieldCalc.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://archermarx.github.io/LoopFieldCalc.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/archermarx/LoopFieldCalc.jl",
)
