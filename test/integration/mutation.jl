using Random

@testset verbose = true "Mutation" begin
    function test_mutation(name, atom_config, method; preserve=false)
        rng = Xoshiro(123)

        # create sample population
        box_size = Vec(8, 8, 8)
        population_config = InitialPopulationConfiguration(10, box_size, atom_config)
        distance_thresholds = DistanceThresholds(0.6, 5)
        population = MOLGA.GeneticAlgorithm.InitialPopulation.create(
            population_config, distance_thresholds; rng
        )

        # mutate all structures (propability=1) and keep original ones
        mutation_config = MutationConfiguration(1, true)
        MOLGA.GeneticAlgorithm.Mutation.mutate!(
            population,
            mutation_config,
            preserve,
            box_size,
            distance_thresholds;
            mutation_method_distribution=method == "atom_switch" ? [1, 0] : [0, 1],
            rng,
        )

        # reorder population so that you have pairs of the original and mutated structures
        function reorder(v)
            n = length(v)
            half_n = div(n, 2)
            reordered = Vector{eltype(v)}(undef, n)

            for i in 1:half_n
                reordered[2 * i - 1] = v[half_n + i]
                reordered[2 * i] = v[i]
            end

            return reordered
        end
        population = reorder(population)

        # export to multiple xyz file
        filename = "$(name)-$(method)$(preserve ? "-preserve" : "")"
        testfile_name = "temp-testfile-mutation-$(filename).xyz"
        MOLGA.IOm.export_xyz(population, testfile_name)

        # check if equal to manually revised file
        return success(`diff $testfile_name integration/xyz/mutation-$filename.xyz`)
    end

    water_atoms = [
        BaseAtom(8, Vec(-1.674872668, 0.000000000, -0.984966492))
        BaseAtom(1, Vec(-1.674872668, 0.759337000, -0.388923492))
        BaseAtom(1, Vec(-1.674872668, -0.759337000, -0.388923492))
    ]

    c6h6_config = [ConfigAtom(6, 6), ConfigAtom(6, 1)]
    @test test_mutation("C6H6", c6h6_config, "atom_switch")
    @test test_mutation("C6H6", c6h6_config, "atom_switch"; preserve=true)
    @test test_mutation("C6H6", c6h6_config, "atom_reposition")
    @test test_mutation("C6H6", c6h6_config, "atom_reposition"; preserve=true)

    o2h2o_config = [ConfigAtom(2, 8), ConfigMolecule(4, water_atoms)]
    @test test_mutation("O2H2O", o2h2o_config, "atom_switch")
    @test test_mutation("O2H2O", o2h2o_config, "atom_switch"; preserve=true)
    @test test_mutation("O2H2O", o2h2o_config, "atom_reposition")
    @test test_mutation("O2H2O", o2h2o_config, "atom_reposition"; preserve=true)

    agh2o_config = [ConfigAtom(1, 47), ConfigMolecule(5, water_atoms)]
    @test test_mutation("AgH2O", agh2o_config, "atom_switch")
    @test test_mutation("AgH2O", agh2o_config, "atom_switch"; preserve=true)
    @test test_mutation("AgH2O", agh2o_config, "atom_reposition")
    @test test_mutation("AgH2O", agh2o_config, "atom_reposition"; preserve=true)
end
