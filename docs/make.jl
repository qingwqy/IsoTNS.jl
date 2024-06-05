using IsoTNS
using Documenter

DocMeta.setdocmeta!(IsoTNS, :DocTestSetup, :(using IsoTNS); recursive=true)

makedocs(;
    modules=[IsoTNS],
    authors="Qingyuan <qingyuan@QingyuandeMacBook-Air.local> and contributors",
    sitename="IsoTNS.jl",
    format=Documenter.HTML(;
        canonical="https://qingwqy.github.io/IsoTNS.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/qingwqy/IsoTNS.jl",
    devbranch="main",
)
