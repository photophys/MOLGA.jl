module InitialPopulation

using Random

using ...Configuration
using ...Types

using ..Utils

"""
    create(config [; rng])

Initialize the population by creating the desired amount of random structures using the
[`InitialPopulationConfiguration`](@ref).

# Arguments

  - `config::`[`ConfigurationObject`](@ref)
  - `rng::AbstractRNG`: Random number generator. If you need consistent results for testing
    purposes, pass a seeded pseudorandom number generator here, eg.
    [`Xoshiro(seed)`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.Xoshiro).
    Defaults to
    [`Random.default_rng()`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.default_rng).

# Returns

A vector of [structures](@ref Structure).
"""
function create(config::ConfigurationObject; rng::AbstractRNG=Random.default_rng())
    population = Structure[]

    indices = generate_indices(config.initial_population.atom_config)

    for _ in 1:(config.initial_population.num_structures)
        atoms = generate_atoms(
            indices, config.initial_population, config.distance_thresholds, 10_000; rng
        )
        push!(population, Structure(; atoms=atoms))
    end

    return population
end

"""
    generate_atoms(indices, population_config, distance_thresholds, max_iterations [; rng])

Distribute the atoms and fragments randomly within a cuboid periodic boundary condition box with
respect to the distance thresholds.

# Arguments

  - `indices::Matrix{Int}`: A two-dimensional matrix generated by [`generate_indices`](@ref)
    representing the atom configuration.
  - `population_config::`[`InitialPopulationConfiguration`](@ref)
  - `distance_thresholds::`[`DistanceThresholds`](@ref)
  - `max_iterations::Int`: The maximum number of iterations to find proper positions while
    distributing the atoms and fragments.
  - `rng::AbstractRNG`: Random number generator. If you need consistent results for testing
    purposes, pass a seeded pseudorandom number generator here, eg.
    [`Xoshiro(seed)`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.Xoshiro).
    Defaults to
    [`Random.default_rng()`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.default_rng).

# Returns

A vector of [atoms](@ref Atom).
"""
function generate_atoms(
    indices::Matrix{Int},
    population_config::InitialPopulationConfiguration,
    distance_thresholds::DistanceThresholds,
    max_iterations::Int;
    rng::AbstractRNG=Random.default_rng(),
)
    atoms = Vector{Atom}()

    for it in 1:max_iterations
        failed = false
        shuffle_rows!(indices; rng)

        for (config_idx, permanent_index) in eachrow(indices)
            result = find_position!(
                population_config.atom_config[config_idx],
                atoms,
                permanent_index,
                config_idx,
                population_config.box_size,
                distance_thresholds,
                max_iterations;
                rng,
            )
            if !result
                empty!(atoms)
                failed = true
                break
            end
        end
        if !failed
            @debug "Iterations needed for generate_atoms", it
            return atoms
        end
    end

    @error "Atom list generation failed."
    return atoms
end

"""
    generate_indices(atom_config)

Create a list of the indices using the [`AtomConfig`](@ref).

The resulting two-dimensional matrix contains the `fragment_kind` in the first column and the
`permanent_index` in the second column. Both indices are integers starting from ``1``. The entries
are inserted as often as specified in the [`AtomConfig`](@ref) `quantity`-property.

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

5×2 Matrix{Int64}:
 1  1
 1  2
 2  3
 3  4
 3  5
```
"""
function generate_indices(atom_config::AtomConfig)
    total_quantity = sum(entry.quantity for entry in atom_config)

    # (fragment_index=atom_config_index, permanent_index)
    indices = Matrix{Int}(undef, total_quantity, 2)

    permanent_index = 1
    for (i, entry) in enumerate(atom_config)
        for _ in 1:(entry.quantity)
            indices[permanent_index, 1] = i
            indices[permanent_index, 2] = permanent_index
            permanent_index += 1
        end
    end

    return indices
end

"""
    find_position!(atom, atom_list, permanent_index, fragment_kind box_size, distance_thresholds,
        max_attempts, [; rng])

Try to find a position for an atom within the periodic boundary condition box by respecting distance
thresholds. If a proper position is found, the list of atoms is updated in-place by adding the
specified atom with its `permanent_index` and `fragment_kind` to the `atom_list`.

# Arguments

  - `atom::`[`ConfigAtom`](@ref): The atom to be placed into the structure.
  - `atom_list::Vector{`[`Atom`](@ref)`}`: The list of atoms of the structure.
  - `permanent_index::Int`
  - `fragment_kind::Int`
  - `box_size::Vec`
  - `distance_thresholds::`[`DistanceThresholds`](@ref)
  - `max_attempts::Int`: Maximum number of attempts to find a proper random position.
  - `rng::AbstractRNG`: Random number generator. If you need consistent results for testing
    purposes, pass a seeded pseudorandom number generator here, eg.
    [`Xoshiro(seed)`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.Xoshiro).
    Defaults to
    [`Random.default_rng()`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.default_rng).

# Returns

  - `true` if the atom could be placed,
  - `false` if no suitable position could be found within the maximum number of attempts.

# Example

```jldoctest; setup=:(using MOLGA.Types, MOLGA.Configuration, Random; import MOLGA.GeneticAlgorithm.InitialPopulation.find_position!)
julia> atom_list = [
           Atom(8, [-0.299181461, 0.000000000, -0.473896792], 1, 1),
           Atom(1, [-0.299181461, 0.759337000, 0.122146208], 1, 1),
           Atom(1, [-0.299181461, -0.759337000, 0.122146208], 1, 1),
       ];

julia> find_position!(
           ConfigAtom(1, 1),
           atom_list,
           2,
           2,
           Vec([8, 8, 8]),
           DistanceThresholds(1, 5),
           5_000;
           rng=Xoshiro(123),
       )
true

julia> atom_list
4-element Vector{Atom}:
 O     1    1 [-0.299181461         0.0                  -0.473896792        ]

 H     1    1 [-0.299181461         0.759337             0.122146208         ]

 H     1    1 [-0.299181461         -0.759337            0.122146208         ]

 H     2    2 [0.1697103642830644   0.6944540596267874   3.1270295847422487  ]
```
"""
function find_position!(
    atom::ConfigAtom,
    atom_list::Vector{Atom},
    permanent_index::Int,
    fragment_kind::Int,
    box_size::Vec,
    distance_thresholds::DistanceThresholds,
    max_attempts::Int;
    rng::AbstractRNG=Random.default_rng(),
)
    for it in 1:max_attempts
        new_pos = random_position(box_size; rng)
        proper_space = check_distance(new_pos, atom_list, distance_thresholds)

        if proper_space
            push!(atom_list, Atom(atom.element, new_pos, permanent_index, fragment_kind))
            @debug "Iterations needed for find_position!(atom, ...)", it
            return true
        end
    end

    return false
end

"""
    find_position!(fragment, atom_list, permanent_index, fragment_kind, box_size,
        distance_thresholds, max_attempts [; rng])

Place a fragment ([`ConfigMolecule`](@ref)) into the periodic boundary condition box. Perform
distance checks for all atoms of the fragment and assign the specified permanent index and fragment
kind to all atoms. This method upates the `atom_list` in-place.

See [`find_position!(atom, ...)`](@ref find_position!) for more details.
"""
function find_position!(
    fragment::ConfigMolecule,
    atom_list::Vector{Atom},
    permanent_index::Int,
    fragment_kind::Int,
    box_size::Vec,
    distance_thresholds::DistanceThresholds,
    max_attempts::Int;
    rng::AbstractRNG=Random.default_rng(),
)
    fragment_atom_list = Vector{Atom}()

    for it in 1:max_attempts
        center_pos = random_position(box_size; rng)
        angles = random_angles(; rng)

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
            @debug "Iterations needed for find_position!(fragment, ...)", it
            return true
        end
    end

    return false
end

end
