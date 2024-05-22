using MOLGA
using ArgParse

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

PARSED_ARGS = @time parse_args(ARGS, arg_parser)

run_genetic_algorithm(PARSED_ARGS["config"]; debug_log=PARSED_ARGS["debug"])
