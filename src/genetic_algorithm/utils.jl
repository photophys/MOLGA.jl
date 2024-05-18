module Utils

using LinearAlgebra
using Random

using ...Types
using ...Configuration

export check_distance

"""
    check_distance(new_position, atoms, thresholds)

Checks if the new position is sufficiently far away from all existing atoms in the atom list.

# Arguments

  - `new_position::`[`Vec`](@ref): The new position to be checked.
  - `atoms::Vector{<:`[`AbstractAtom`](@ref)`}`: The current list of atoms.
  - `thresholds::`[`DistanceThresholds`](@ref)

# Example

The atom list contains the atoms of a water molecule. The task is to check if the chosen position is
suitable. The distances between the new position and the hydrogen atoms are 1.80392 and 1.30274; the
distance between `pos` and the carbon atom is 1.40762 (all measured with Chemcraft).

```jldoctest; setup=:(using MOLGA.Types, MOLGA.Configuration, MOLGA.GeneticAlgorithm.Utils)
julia> pos = Vec([0.959329065, 0.512611795, -0.106788996]);

julia> atoms = [
           BaseAtom(8, [-0.299181461, 0.000000000, -0.473896792]),
           BaseAtom(1, [-0.299181461, 0.759337000, 0.122146208]),
           BaseAtom(1, [-0.299181461, -0.759337000, 0.122146208]),
       ];

julia> check_distance(pos, atoms, DistanceThresholds(0.5, 5.0))
true

julia> check_distance(pos, atoms, DistanceThresholds(1.4, 5.0))
false

julia> check_distance(pos, atoms, DistanceThresholds(0.5, 1))
false
```
"""
function check_distance(
    new_position::Vec, atoms::Vector{<:AbstractAtom}, thresholds::DistanceThresholds
)
    if length(atoms) < 1
        @debug "No distance check because atom list has less than one atom"
        return true
    end

    for atom in atoms
        distance = norm(new_position - atom.position)
        if distance < thresholds.min || distance > thresholds.max
            return false
        end
    end

    return true
end

end
