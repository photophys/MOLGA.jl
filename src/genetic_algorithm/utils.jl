module Utils

using LinearAlgebra
using Distributions
using Random

using ...Types
using ...Configuration

export check_distance, shuffle_rows!, rotate_position, random_position, random_angles

"""
    check_distance(new_position, atoms, thresholds, check_only_lower_threshold=false)

Checks if the new position is sufficiently far away from all existing atoms in the atom list.

# Arguments

  - `new_position::`[`Vec`](@ref): The new position to be checked.
  - `atoms::Vector{<:`[`AbstractAtom`](@ref)`}`: The current list of atoms.
  - `thresholds::`[`DistanceThresholds`](@ref)
  - `check_only_lower_threshold::Bool`: True if only the lower distance threshold should be checked.

# Example

The atom list contains the atoms of a water molecule. The task is to check if the chosen position is
suitable. The distances between the new position and the hydrogen atoms are 1.80392 and 1.30274; the
distance between `pos` and the carbon atom is 1.40762 (all measured with Chemcraft).

```jldoctest; setup=:(using MOLGA.Types, MOLGA.Configuration, MOLGA.GeneticAlgorithm.Utils)
julia> pos = Vec([0.959329065, 0.512611795, -0.106788996]);

julia> check_distance(pos, Atom[], DistanceThresholds(0.5, 5.0))
true

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

julia> check_distance(pos, atoms, DistanceThresholds(0.5, 1), true)
true
```
"""
function check_distance(
    new_position::Vec,
    atoms::Vector{<:AbstractAtom},
    thresholds::DistanceThresholds,
    check_only_lower_threshold::Bool=false,
)
    if length(atoms) < 1
        @debug "No distance check because atom list has less than one atom"
        return true
    end

    for atom in atoms
        distance = norm(new_position - atom.position)
        if distance < thresholds.min || (!check_only_lower_threshold && distance > thresholds.max)
            return false
        end
    end

    return true
end

"""
    shuffle_rows!(m::AbstractMatrix [; rng])

Randomly reorder the rows of the given matrix in-place. You can pass a random number generator to
`rng` like described [here](@ref create).

# Example

```jldoctest; setup=:(import MOLGA.GeneticAlgorithm.InitialPopulation.shuffle_rows!; using Random)
julia> mat = [1 10; 2 20; 3 30; 4 40; 5 50];

julia> shuffle_rows!(mat; rng=Xoshiro(123));

julia> mat
5×2 Matrix{Int64}:
 5  50
 4  40
 2  20
 3  30
 1  10
```
"""
function shuffle_rows!(m::AbstractMatrix; rng::AbstractRNG=Random.default_rng())
    row_order = randperm(rng, size(m, 1))
    m .= m[row_order, :]

    return m
end

"""
    rotate_position(position, angle)

Rotate a 3D position vector by a 3D angle vector around the origin using
[rotation matrices](https://en.wikipedia.org/wiki/Rotation_matrix#Basic_3D_rotations).

# Arguments

  - `position::Vector{Float64}`: The 3D position vector to be rotated.
  - `angle::Vector{Float64}`: The 3D angle vector representing rotation angles around the ``x``,
    ``y``, and ``z`` axes in radians.

# Example

```jldoctest; setup=:(import MOLGA.GeneticAlgorithm.InitialPopulation.rotate_position; using MOLGA.Types)
julia> rotate_position(Vec([1, 0, 0.5]), Vec([0, 0, π]))
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -1.0
  1.2246467991473532e-16
  0.5
```
"""
function rotate_position(position::Vec, angle::Vec)
    a, b, c = angle

    R = [
        cos(b)cos(c) cos(c)sin(a)sin(b)-cos(a)sin(c) cos(a)cos(c)sin(b)+sin(a)sin(c)
        cos(b)sin(c) cos(a)cos(c)+sin(a)sin(b)sin(c) -cos(c)sin(a)+cos(a)sin(b)sin(c)
        -sin(b) cos(b)sin(a) cos(a)cos(b)
    ]

    return Vec(R * position)
end

"""
    random_position(box_size [; rng])

Generate a random position within the specified box. The origin is the box's center point, so that
a box size of ``\\mathbf{s}=\\left(s_x,s_y,s_z\\right)`` leads to a random position ``\\mathbf{x}``

```math
-\\frac{1}{2}\\,\\mathbf{s}\\leq\\mathbf{x}\\leq\\frac{1}{2}\\,\\mathbf{s} \\text{.}
```

# Example

```jldoctest; setup=:(import MOLGA.GeneticAlgorithm.InitialPopulation.random_position; using MOLGA.Types; using Random)
julia> random_position(Vec([8, 5, 5]); rng=Xoshiro(1))
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -3.413069164245657
 -0.7537925522140694
  0.9941334184573423
```
"""
function random_position(box_size::Vec; rng::AbstractRNG=Random.default_rng())
    return Vec(rand(rng, Uniform(-s / 2, s / 2)) for s in box_size)
end

"""
    random_angles([; rng])

Generate a vector of three random angles (from ``0`` to ``2\\pi``, in radians).

# Example

```jldoctest; setup=:(import MOLGA.GeneticAlgorithm.InitialPopulation.random_angles; using Random)
julia> random_angles(; rng=Xoshiro(123))
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 3.2748828620072237
 3.687015596584574
 5.597555946335841
```
"""
random_angles(; rng::AbstractRNG=Random.default_rng()) = Vec(rand(rng, 3) * 2π)

end
