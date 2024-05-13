module Types

using StaticArrays
using DocStringExtensions
using PeriodicTable

using ..Helpers

export AbstractAtom, Atom, BaseAtom, Structure, Vec

"""
Three-dimensional vector.
"""
const Vec = SVector{3,Float64}

function Base.show(io::IO, vec::Vec)
    print(io, "[", join(rpad.(vec, 20), " "), "]")
    return nothing
end

abstract type AbstractAtom end

"""
    BaseAtom(element, position)

Stores a basic atom.
$(TYPEDFIELDS)
"""
struct BaseAtom <: AbstractAtom
    """Atomic number of the element, eg. 6 for Carbon."""
    element::Int

    """3D position [vector](@ref Vec) of the atom in angstroms."""
    position::Vec
end

function Base.show(io::IO, atom::BaseAtom)
    return println(io, element_symbol(atom.element), " ", atom.position)
end

"""
    Atom(element, position, permanent_index, fragment_kind)

Stores properties of a single atom.
$(TYPEDFIELDS)
"""
struct Atom <: AbstractAtom
    """Atomic number of the element, eg. 6 for Carbon."""
    element::Int

    """3D position [vector](@ref Vec) of the atom in angstroms."""
    position::Vec

    """Unique index attached to an atom."""
    permanent_index::Int

    """
    Kind of fragment this atom belongs to, e.g. all atoms that belong to a H2O molecule have the
    same value here, also all single O atoms have the same value here, respectively.
    """
    fragment_kind::Int
end

function Base.show(io::IO, atom::Atom; table::Bool=false)
    if table
        println(
            io,
            rpad(element_symbol(atom.element), 2),
            " ",
            lpad(atom.fragment_kind, 4),
            " ",
            lpad(atom.permanent_index, 4),
            " ",
            atom.position,
        )
    else
        println(
            io,
            elements[atom.element].name,
            " fragment_kind=",
            atom.fragment_kind,
            " permanent_index=",
            atom.permanent_index,
            " position=",
            atom.position,
        )
    end

    return nothing
end

function Base.show(io::IO, atom_list::AbstractVector{Atom})
    println(io, rpad("", 3), "frag", " ", "perm", " ", "position")
    for atom in atom_list
        show(io, atom; table=true)
    end
    return nothing
end

"""
    Structure(atoms)
    Structure(atoms [, energy [, nuclear_repulsion_energy [, age]]])

Stores properties of a molecular structure.
$(TYPEDFIELDS)
"""
@kwdef mutable struct Structure
    """Array of [atoms](@ref Atom) that the structure consists of."""
    atoms::Vector{Atom}

    """Energy of the structure in eV."""
    energy::Float64 = 0.0

    """Nuclear repulsion energy of the structure in eV."""
    nuclear_repulsion_energy::Float64 = 0.0

    """Number of cycles that this structure is already present."""
    age::Int = 0
end

function Base.show(io::IO, structure::Structure)
    println(io, "Structure with $(length(structure.atoms)) atoms")
    println(io, "Energy: ", structure.energy)
    println(io, "Nuclear repulsion energy: ", structure.nuclear_repulsion_energy)
    println(io, "Age: ", structure.age)
    println(io, structure.atoms)

    return nothing
end

function Base.show(io::IO, structure_list::AbstractVector{Structure})
    println(io, "List with ", length(structure_list), " structures")

    for structure in structure_list
        println(io, structure)
    end

    return nothing
end

end
