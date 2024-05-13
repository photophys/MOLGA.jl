module MOLGA

include("./logging.jl")

include("./helpers.jl")
include("./types.jl")
include("./io.jl")

include("./configuration.jl")
include("./interfaces/_index.jl")
include("./genetic_algorithm/_index.jl")

export run

"""
    run(configuration_file [, debug_log::Bool = false])

Start the Molecular Genetic Algorithm.
"""
function run(configuration_file::String, debug_log::Bool=false)
    Logging.setup_logger(debug_log)

    CONFIG = Configuration.load_config(configuration_file)
    population = GeneticAlgorithm.InitialPopulation.create(CONFIG)

    @info population

    return nothing
end

end

using ArgParse

if abspath(PROGRAM_FILE) == @__FILE__
    arg_parser = ArgParseSettings()
    @add_arg_table! arg_parser begin
        "--config", "-c"
        help = "path to the configuration file"
        arg_type = String
        default = "config.yaml"
        metavar = "path/to/config.yaml"

        "--debug", "-d"
        help = "enable logging of debug-level events"
        action = :store_true
    end

    PARSED_ARGS = parse_args(ARGS, arg_parser)

    MOLGA.run(PARSED_ARGS["config"], PARSED_ARGS["debug"])
end
