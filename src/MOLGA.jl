module MOLGA

using Random

include("./logging.jl")

include("./helpers.jl")
include("./types.jl")

include("./configuration.jl")
include("./interfaces/.module.jl")
include("./genetic_algorithm/.module.jl")

export run_genetic_algorithm

"""
    run_genetic_algorithm(configuration_object [; debug_log = false, rng])

Run the molecular genetic algorithm using the provided
[`ConfigurationObject`](@ref Configuration.ConfigurationObject).

# Parameters

  - `configuration_object::`[`ConfigurationObject`](@ref Configuration.ConfigurationObject)
  - `debug_log::Bool`: Defines if you want to enable the logging of debug-level events.
  - `rng::AbstractRNG`: Random number generator. If you need consistent results for testing
    purposes, pass a seeded pseudorandom number generator here, eg.
    [`Xoshiro(seed)`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.Xoshiro).
    Defaults to
    [`Random.default_rng()`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.default_rng).
"""
function run_genetic_algorithm(
    config::Configuration.ConfigurationObject;
    debug_log::Bool=false,
    rng::AbstractRNG=Random.default_rng(),
)
    Logging.setup_logger(debug_log)

    population = GeneticAlgorithm.InitialPopulation.create(
        config.initial_population, config.distance_thresholds; rng
    )

    for big_cycle_idx in 1:1
        GeneticAlgorithm.Mutation.mutate!(
            population,
            config.mutation,
            config.preserve_fragments,
            config.initial_population.box_size,
            config.distance_thresholds;
            rng,
        )
    end

    @info length(population)
    Helpers.export_xyz(population, "results.xyz")

    return nothing
end

"""
run_genetic_algorithm(configuration_file [, debug_log = false, rng])

Start the molecular genetic algorithm and load the necessary configuration parameters from the
specified TOML file.

# Parameters

  - `configuration_file::String`: The filename and path to the TOML configuration file (see
    specification [here](parameters/input-file.md)).
  - `debug_log::Bool`: Defines if you want to enable the logging of debug-level events.
  - `rng::AbstractRNG`: Random number generator. If you need consistent results for testing
    purposes, pass a seeded pseudorandom number generator here, eg.
    [`Xoshiro(seed)`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.Xoshiro).
    Defaults to
    [`Random.default_rng()`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.default_rng).
"""
function run_genetic_algorithm(
    configuration_file::String; debug_log::Bool=false, rng::AbstractRNG=Random.default_rng()
)
    Logging.setup_logger(debug_log)

    config = Configuration.load_config(configuration_file)
    run_genetic_algorithm(config; debug_log, rng)

    return nothing
end

end
