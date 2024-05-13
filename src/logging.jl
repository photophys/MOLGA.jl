module Logging

using LoggingExtras
using Dates

export setup_logger

const DATE_FORMAT = "yyyy-mm-dd HH:MM:SS"

function timestamp_logger(logger)
    TransformerLogger(logger) do log
        merge(log, (; message="$(Dates.format(now(), DATE_FORMAT)) $(log.message)"))
    end
end

function setup_logger(enable_debug::Bool)
    if enable_debug
        global_logger(timestamp_logger(ConsoleLogger(stdout, Base.CoreLogging.Debug)))
        @info "Enabled debug log"
    else
        global_logger(timestamp_logger(ConsoleLogger(stdout)))
    end
end

end
