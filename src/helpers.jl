module Helpers

using PeriodicTable

using ..Types

export element_symbol, atomic_number, load_xyz, export_xyz

"""
    element_symbol(atomic_number)

Get the symbol of an element of the periodic table as `String`.

# Arguments

  - `atomic_number::Int`

# Example

```jldoctest; setup=:(using MOLGA.Helpers)
julia> element_symbol(1)
"H"
```
"""
element_symbol(atomic_number::Int) = String(elements[atomic_number].symbol)

"""
    atomic_number(element_symbol)

Get the atomic number of an element of the periodic table as `Int`.

# Arguments

  - `element_symbol::AbstractString`: The case-sensitive element symbol.

# Example

```jldoctest; setup=:(using MOLGA.Helpers)
julia> atomic_number("Ar")
18
```
"""
atomic_number(element_symbol::AbstractString) = elements[Symbol(element_symbol)].number

"""
    load_xyz(xyz_file::String)

Get all atoms from a structure in XYZ-format (with element symbols).
"""
function load_xyz(xyz_file::String)
    @debug "Read atoms from xyz-file $(xyz_file)"

    atoms = BaseAtom[]

    file = open(xyz_file, "r")

    try
        for line in Iterators.drop(eachline(file), 2)
            parts = split(line)
            push!(atoms, BaseAtom(atomic_number(parts[1]), parse.(Float64, parts[2:4])))
        end
    finally
        close(file)
    end

    return atoms
end

"""
    export_xyz(structures, xyz_file)

Export a set of structures to a xyz files (with element symbols).

# Arguments

  - `structures::AbstractVector{`[`Structure`](@ref)`}`
  - `xyz_file::String`: The expected name/path of the file.
"""
function export_xyz(structures::AbstractVector{Structure}, xyz_file::String)
    @debug "Save structures to $(xyz_file)"

    file = open(xyz_file, "w")
    write(xs::Vararg{Any}) = println(file, xs...)

    function format(x::AbstractFloat)
        raw = split(string(x), ".")
        return lpad(raw[1], 3) * "." * rpad(raw[2], 17)
    end

    for structure in structures
        write(length(structure.atoms))
        write(
            "energy ",
            structure.energy,
            ", nre ",
            structure.nuclear_repulsion_energy,
            ", age ",
            structure.age,
        )

        for atom in structure.atoms
            write(rpad(element_symbol(atom.element), 2), " ", join(format.(atom.position), " "))
        end
        write()
    end

    close(file)

    return nothing
end

end
