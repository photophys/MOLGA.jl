using MOLGA
using Documenter

DEMO_STRUCTURE_Ag = MOLGA.Types.Structure(;
    atoms=[
        MOLGA.Types.Atom(8, [1.834176000, -1.086997000, -0.145691000], 2, 2),
        MOLGA.Types.Atom(1, [2.642219000, -0.527623000, 0.117024000], 2, 2),
        MOLGA.Types.Atom(1, [2.082181000, -1.613762000, -0.930716000], 2, 2),
        MOLGA.Types.Atom(8, [3.732413000, 0.511988000, 0.663892000], 3, 2),
        MOLGA.Types.Atom(1, [4.173112000, 0.327927000, 1.518267000], 3, 2),
        MOLGA.Types.Atom(1, [4.423184000, 0.920848000, 0.104580000], 3, 2),
        MOLGA.Types.Atom(47, [0.000012000, -0.000259000, -0.210701000], 1, 1),
        MOLGA.Types.Atom(8, [-1.833937000, 1.086716000, -0.146848000], 4, 2),
        MOLGA.Types.Atom(1, [-2.081857000, 1.612214000, -0.932792000], 4, 2),
        MOLGA.Types.Atom(1, [-2.642220000, 0.527829000, 0.116511000], 4, 2),
        MOLGA.Types.Atom(8, [-3.732675000, -0.510572000, 0.664366000], 5, 2),
        MOLGA.Types.Atom(1, [-4.423943000, -0.919580000, 0.105806000], 5, 2),
        MOLGA.Types.Atom(1, [-4.173075000, -0.324762000, 1.518501000], 5, 2),
    ],
)

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
