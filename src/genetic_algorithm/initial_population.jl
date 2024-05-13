module InitialPopulation

using Random
using Distributions

using ...Configuration
using ...Types

using ..Utils

"""
    create(config)

Initialize the population by creating the desired amount of random structures using the [`InitialPopulationConfiguration`](@ref).

# Arguments

  - `config::`[`ConfigurationObject`](@ref)

# Returns

A vector of [structures](@ref Structure).
"""
function create(config::ConfigurationObject)
    population = Structure[]

    for _ in 1:(config.initial_population.num_structures)
        atoms = generate_atoms(config.initial_population, config.distance_thresholds, 10_000)
        push!(population, Structure(; atoms=atoms))
    end

    return population
end

"""
    generate_atoms(population_config, distance_thresholds, max_iterations)

Distribute the atoms and fragments randomly within a cuboid periodic boundary condition box with
respect to the distance thresholds.

# Arguments

  - `population_config::`[`InitialPopulationConfiguration`](@ref)
  - `distance_thresholds::`[`DistanceThresholds`](@ref)
  - `max_iterations::Int`: The maximum number of iterations to find proper positions while
    distributing the atoms and fragments.

# Returns

A vector of [atoms](@ref Atom).
"""
function generate_atoms(
    population_config::InitialPopulationConfiguration,
    distance_thresholds::DistanceThresholds,
    max_iterations::Int,
)
    atoms = Vector{Atom}()
    indices = generate_indices(population_config.atom_config)

    for it in 1:max_iterations
        failed = false
        shuffle!(indices)

        for (config_idx, permanent_index) in indices
            result = find_position!(
                population_config.atom_config[config_idx],
                atoms,
                population_config.box_size,
                distance_thresholds,
                max_iterations;
                permanent_index=permanent_index,
                fragment_kind=config_idx,
            )
            if !result
                empty!(atoms)
                failed = true
                break
            end
        end
        if !failed
            @debug "Needed $(it) iterations for generate_atoms"
            return atoms
        end
    end

    @error "Atom list generation failed."
    return atoms
end

"""
    generate_indices(atom_config)

Create a list of the indices using the [`AtomConfig`](@ref).

The tuples in the list contain the `fragment_kind` first and the `permanent_index` second. Both
indices are integers starting from ``1``. The entries are inserted as often as specified in the
[`AtomConfig`](@ref) `quantity`-property.

Each fragment or _single_ atom (one or more) have their own `fragment_kind` which is equivalent to
the index of the [`AtomConfig`](@ref) list. To provide an option to preserve fragments, all atoms of
an instance of a fragment have a common `permanent_index`.

# Example

Small example distributing two hydrogen atoms, one oxygen atom, as well as two carbon dioxide
fragments for demonstration purposes.

```jldoctest; setup=:(import MOLGA.GeneticAlgorithm.InitialPopulation.generate_indices;using MOLGA.Configuration,MOLGA.Types)
atom_config::AtomConfig = [
    ConfigAtom(2, 1), # 2x H
    ConfigAtom(1, 8), # 1x O
    ConfigMolecule(
        2, # 2x CO2
        [
            BaseAtom(6, [-1.161937300 0.000000000 -1.704638554]),
            BaseAtom(8, [-0.610097929 -0.126832008 -0.681728275]),
            BaseAtom(8, [-1.713776672 0.126832008 -2.727548833]),
        ],
    ),
]
generate_indices(atom_config)

# output

5-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (1, 2)
 (2, 3)
 (3, 4)
 (3, 5)

```
"""
function generate_indices(atom_config::AtomConfig)
    indices = Vector{Tuple{Int,Int}}() # (fragment_index=atom_config_index, permanent_index)
    permanent_index = 1
    for (i, entry) in enumerate(atom_config)
        for _ in 1:(entry.quantity)
            push!(indices, (i, permanent_index))
            permanent_index += 1
        end
    end

    return indices
end

"""
    find_position!(atom, atom_list, box_size, distance_thresholds, max_attempts;
        permanent_index, fragment_kind)

Try to find a position for an atom within the periodic boundary condition box by respecting distance
thresholds. If a proper position is found, the list of atoms is updated in-place by adding the
specified atom with its `permanent_index` and `fragment_kind` to the `atom_list`.

# Arguments

  - `atom::`[`ConfigAtom`](@ref): The atom to be placed into the structure.
  - `atom_list::Vector{`[`Atom`](@ref)`}`: The list of atoms of the structure.
  - `box_size::AbstractVector`
  - `distance_thresholds::`[`DistanceThresholds`](@ref)
  - `max_attempts::Int`: Maximum number of attempts to find a proper random position.
  - `permanent_index::Int`
  - `fragment_kind::Int`

# Returns

  - `true` if the atom could be placed,
  - `false` if no suitable position could be found within the maximum number of attempts.
"""
function find_position!(
    atom::ConfigAtom,
    atom_list::Vector{Atom},
    box_size::AbstractVector,
    distance_thresholds::DistanceThresholds,
    max_attempts::Int;
    permanent_index::Int,
    fragment_kind::Int,
)
    for it in 1:max_attempts
        new_pos = random_position(box_size) # TODO: only update new_pos value, but don't allocate a new one in each iteration
        proper_space = check_distance(new_pos, atom_list, distance_thresholds)

        if proper_space
            push!(atom_list, Atom(atom.element, new_pos, permanent_index, fragment_kind))
            @debug "Needed $(it) iterations for find_position!(new_atom, ...)"
            return true
        end
    end

    return false
end

"""
    find_position!(fragment, atom_list, box_size, distance_thresholds, max_attempts;
        permanent_index, fragment_kind)

Place a fragment ([`ConfigMolecule`](@ref)) into the periodic boundary condition box. Perform
distance checks for all atoms of the fragment and assign the specified permanent index and fragment
kind to all atoms. This method upates the `atom_list` in-place.

See [`find_position!(atom, ...)`](@ref find_position!) for more details.
"""
function find_position!(
    fragment::ConfigMolecule,
    atom_list::Vector{Atom},
    box_size::AbstractVector,
    distance_thresholds::DistanceThresholds,
    max_attempts::Int;
    permanent_index::Int,
    fragment_kind::Int,
)
    fragment_atom_list = Vector{Atom}()

    for it in 1:max_attempts
        center_pos = random_position(box_size)
        angles = random_angles()

        proper_space = false

        for atom in fragment.atoms
            new_pos = rotate_position(atom.position, angles) + center_pos
            if check_distance(new_pos, atom_list, distance_thresholds)
                push!(
                    fragment_atom_list, Atom(atom.element, new_pos, permanent_index, fragment_kind)
                )
                proper_space = true
            else
                empty!(fragment_atom_list)
                proper_space = false
                break
            end
        end

        if proper_space
            append!(atom_list, fragment_atom_list)
            @debug "Needed $(it) iterations for find_position!(fragment)"
            return true
        end
    end

    return false
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

```jldoctest; setup=:(import MOLGA.GeneticAlgorithm.InitialPopulation.rotate_position)
julia> rotate_position([1, 0, 0.5], [0, 0, π])
3-element Vector{Float64}:
 -1.0
  1.2246467991473532e-16
  0.5
```
"""
function rotate_position(position::AbstractVector, angle::AbstractVector)
    # TODO: use pre-computed form to avoid matrix multiplications ?
    Rx = [1 0 0; 0 cos(angle[1]) -sin(angle[1]); 0 sin(angle[1]) cos(angle[1])]
    Ry = [cos(angle[2]) 0 sin(angle[2]); 0 1 0; -sin(angle[2]) 0 cos(angle[2])]
    Rz = [cos(angle[3]) -sin(angle[3]) 0; sin(angle[3]) cos(angle[3]) 0; 0 0 1]

    R = Rz * Ry * Rx

    return R * position
end

"""
    random_position(box_size)

Generate a random position within the specified box. The origin is the box's center point, so that
a box size of ``\\mathbf{s}=\\left(s_x,s_y,s_z\\right)`` leads to a random position ``\\mathbf{x}``

```math
-\\frac{1}{2}\\,\\mathbf{s}\\leq\\mathbf{x}\\leq\\frac{1}{2}\\,\\mathbf{s} \\text{.}
```

# Example

```@repl
julia> random_position([8, 5, 5])
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
  3.6197648788491597
  1.4205599794740316
 -1.2442172005294312
```
"""
random_position(box_size::AbstractVector) = Vec(rand(Uniform(-s / 2, s / 2)) for s in box_size)

"""
    random_angles()

Generate a vector of three random angles (from ``0`` to ``2\\pi``, in radians).

# Example

```@repl
julia> random_angles()
3-element Vector{Float64}:
 5.892339839566983
 3.042066787260773
 0.8485697032585366
```
"""
random_angles() = rand(3) * 2π

end
