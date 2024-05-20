using MOLGA, MOLGA.Types
using Test, Documenter

@testset verbose = true "MOLGA.jl" begin
    @testset "IO" begin
        @test MOLGA.IO.load_xyz("samples/co2.xyz") == [
            BaseAtom(6, [-1.1619373, 0.0, -1.704638554]),
            BaseAtom(8, [-0.610097929, -0.126832008, -0.681728275]),
            BaseAtom(8, [-1.713776672, 0.126832008, -2.727548833]),
        ]

        atoms = [
            Atom(8, [-1.674872668, 0.000000000, -0.984966492], 1, 1),
            Atom(1, [-1.674872668, 0.759337000, -0.388923492], 2, 2),
            Atom(1, [-1.674872668, -0.759337000, -0.388923492], 3, 2),
        ]
        MOLGA.IO.export_xyz(
            [
                Structure(; atoms, energy=56.13537, age=10),
                Structure(; atoms, nuclear_repulsion_energy=1.006),
            ],
            "temp-h2o.xyz",
        )
        @test success(`diff temp-h2o.xyz 'samples/h2o-structure-example.xyz'`)
    end

    doctest(MOLGA)
end
