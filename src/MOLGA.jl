module MOLGA

include("./logging.jl")

include("./helpers.jl")
include("./types.jl")
include("./samples.jl")
include("./io.jl")

include("./configuration.jl")
include("./interfaces/.module.jl")
include("./genetic_algorithm/.module.jl")

export run_genetic_algorithm

"""
    run_genetic_algorithm(configuration_object [; debug_log = false])

Run the molecular genetic algorithm using the provided
[`ConfigurationObject`](@ref Configuration.ConfigurationObject).

# Parameters

  - `configuration_object::`[`ConfigurationObject`](@ref Configuration.ConfigurationObject)
  - `debug_log::Bool`: Defines if you want to enable the logging of debug-level events.
"""
function run_genetic_algorithm(config::Configuration.ConfigurationObject; debug_log::Bool=false)
    Logging.setup_logger(debug_log)

    population = GeneticAlgorithm.InitialPopulation.create(config)

    for big_cycle_idx in 1:1
        #  genetic algorithm here
    end

    @info length(population)

    return nothing
end

"""
run_genetic_algorithm(configuration_file [, debug_log = false])

Start the molecular genetic algorithm and load the necessary configuration parameters from the
specified TOML file.

# Parameters

  - `configuration_file::String`: The filename and path to the TOML configuration file (see
    specification [here](parameters/input-file.md)).
  - `debug_log::Bool`: Defines if you want to enable the logging of debug-level events.
"""
function run_genetic_algorithm(configuration_file::String; debug_log::Bool=false)
    Logging.setup_logger(debug_log)

    config = Configuration.load_config(configuration_file)
    run_genetic_algorithm(config; debug_log)

    return nothing
end

end
