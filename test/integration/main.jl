using Random

@testset verbose = true "Main" begin
    run_genetic_algorithm("integration/config/example.toml")
    @test success(`diff results.xyz integration/xyz/main.xyz`)
end
