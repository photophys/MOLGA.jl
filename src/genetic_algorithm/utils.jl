module Utils

using LinearAlgebra
using Random

using ...Types
using ...Configuration

export check_distance

"""
    check_distance(new_position, atoms, thresholds)

Checks if the new position is sufficiently far away from all existing atoms in the atom list.
"""
function check_distance(
    new_position::AbstractArray, atoms::Vector{<:AbstractAtom}, thresholds::DistanceThresholds
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
