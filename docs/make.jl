using MOLGA
using Documenter

DocMeta.setdocmeta!(MOLGA, :DocTestSetup, :(using MOLGA); recursive=true)

makedocs(;
    modules=[MOLGA],
    authors="Michael Gatt, Gabriel Schöpfer et al.",
    sitename="MOLGA.jl",
    format=Documenter.HTML(;
        canonical="https://photophys.github.io/MOLGA.jl",
        edit_link="main",
        assets=["assets/favicon.ico"],
    ),
    pages=[
        "Home" => "index.md",
        "Features" => [
            "Overview" => "features/overview.md",
            "features/genetic-algorithm.md",
            "features/structure-clustering.md",
            "Interfaces" => "features/interfaces.md",
        ],
        "Getting Started" => ["getting-started/installation.md", "getting-started/tutorial.md"],
        "Parameters" => ["parameters/input-file.md", "parameters/struct.md"],
        "reference.md",
    ],
)

deploydocs(; repo="github.com/photophys/MOLGA.jl", devbranch="main")
