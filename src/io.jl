module IO

using PeriodicTable

using ..Types
using ..Helpers

export load_xyz

"""
    load_xyz(xyz_file::String)

Get all atoms from a structure in XYZ-format (with element symbols).
"""
function load_xyz(xyz_file::String)
    @debug "Read atoms from xyz-file $(xyz_file)"

    atoms = BaseAtom[]

    file = open(xyz_file, "r")

    try
        for line in eachline(file)
            parts = split(line)
            push!(atoms, BaseAtom(atomic_number(parts[1]), parse.(Float64, parts[2:4])))
        end
    finally
        close(file)
    end

    return atoms
end

end
