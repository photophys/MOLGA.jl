using MOLGA, MOLGA.Types, MOLGA.Configuration
using Test, Documenter

@testset verbose = true "MOLGA.jl" begin
    @testset "Unit Tests" begin
        include("./unit/helpers.jl")
    end

    @testset "Integration Tests" begin
        include("./integration/main.jl")
        include("./integration/initial_population.jl")
        include("./integration/mutation.jl")
    end

    doctest(MOLGA)
end
