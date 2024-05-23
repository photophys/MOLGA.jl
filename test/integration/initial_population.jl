using Random

@testset verbose = true "Initial Population" begin
    function test_population(name, atom_config)
        config = InitialPopulationConfiguration(10, Vec(8, 8, 8), atom_config)
        population = MOLGA.GeneticAlgorithm.InitialPopulation.create(
            config, DistanceThresholds(0.6, 5); rng=Xoshiro(123)
        )

        filename = "temp-testfile-population-$(name).xyz"

        MOLGA.IOm.export_xyz(population, filename)

        # check if files are equal
        return success(`diff $filename integration/xyz/population-$name.xyz`)
    end

    water_atoms = [
        BaseAtom(8, Vec(-1.674872668, 0.000000000, -0.984966492))
        BaseAtom(1, Vec(-1.674872668, 0.759337000, -0.388923492))
        BaseAtom(1, Vec(-1.674872668, -0.759337000, -0.388923492))
    ]

    @test test_population("C6H6", [ConfigAtom(6, 6), ConfigAtom(6, 1)])
    @test test_population("O2H2O", [ConfigAtom(2, 8), ConfigMolecule(4, water_atoms)])
    @test test_population("AgH2O", [ConfigAtom(1, 47), ConfigMolecule(5, water_atoms)])
end
