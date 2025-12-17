using CFAnelastics
using Documenter

DocMeta.setdocmeta!(CFAnelastics, :DocTestSetup, :(using CFAnelastics); recursive=true)

makedocs(;
    modules=[CFAnelastics],
    authors="Thomas Dubos <thomas.dubos@polytechnique.edu> and contributors",
    sitename="CFAnelastics.jl",
    format=Documenter.HTML(;
        canonical="https://dubosipsl.github.io/CFAnelastics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dubosipsl/CFAnelastics.jl",
    devbranch="main",
)
