using MOLGA
using Test

@testset verbose = true "MOLGA.jl" begin
    @testset "IO" begin
        @test MOLGA.IO.load_xyz("samples/co2.xyz") == [
            MOLGA.Types.BaseAtom(6, [-1.1619373, 0.0, -1.704638554]),
            MOLGA.Types.BaseAtom(8, [-0.610097929, -0.126832008, -0.681728275]),
            MOLGA.Types.BaseAtom(8, [-1.713776672, 0.126832008, -2.727548833]),
        ]
    end

    @testset "GA-Utils" begin
        pos = MOLGA.Types.Vec([1.282131000, 0.737863000, -0.041933000])
        atoms = [MOLGA.Types.BaseAtom(6, [1.619507883, -0.147537000, -0.210363518])]

        # distance between both atoms is 0.96235 (Chemcraft)
        test_check_distance(min::Float64, max::Float64) =
            MOLGA.GeneticAlgorithm.Utils.check_distance(
                pos, atoms, MOLGA.Configuration.DistanceThresholds(min, max)
            )

        @test test_check_distance(0.5, 5.0) == true
        @test test_check_distance(1.0, 5.0) == false
        @test test_check_distance(0.5, 0.75) == false
    end
end
