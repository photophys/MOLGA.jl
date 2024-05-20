module Mutation

using Random

using ...Configuration
using ...Types
using ..Utils

"""
    mutate!(population, config [; rng])

Perform mutation on all structures of the population in-place.

Currently, there are two mutations (atom switch, atom reposition). Each mutation function is applied
individually to each structure with the specified mutation probability, i.e. it is possible for a
structure to be mutated twice (first atom switch, then atom reposition).

If the `keep_non_mutated` switch is activated, the original structure is added to the population
before performing the mutation(s).

Depending on the choice of the `preserve_fragments` option, the appropriate function is selected for
both mutations.

# Arguments

  - `population::Vector{`[`Structure`](@ref)`}`
  - `config::`[`ConfigurationObject`](@ref)
  - `rng::AbstractRNG`: Random number generator. If you need consistent results for testing
    purposes, pass a seeded pseudorandom number generator here, eg.
    [`Xoshiro(seed)`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.Xoshiro).
    Defaults to
    [`Random.default_rng()`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.default_rng).
"""
function mutate!(
    population::Vector{Structure},
    config::ConfigurationObject;
    rng::AbstractRNG=Random.default_rng(),
)
    original_structures = Vector{Structure}()

    for structure in population
        # TODO: do we want to control each mutation probability separately?
        should_mutate = rand(2) .< config.mutation.probability

        if config.mutation.keep_non_mutated && any(should_mutate)
            push!(original_structures, deepcopy(structure))
        end

        if should_mutate[1]
            if config.preserve_fragments
                switch_atom_preserving!(structure, config.distance_thresholds; rng)
            else
                switch_atom!(structure; rng)
            end
        end

        if should_mutate[2]
            if config.preserve_fragments
                reposition_atom_preserving!(
                    structure, config.initial_population.box_size, config.distance_thresholds; rng
                )
            else
                reposition_atom!(structure, config.initial_population.box_size; rng)
            end
        end
    end

    append!(population, original_structures)

    return nothing
end

"""
    switch_atom!(structure [; rng])

Exchange the position of two randomly selected atoms (of different chemical elements) from the
provided structure in-place without preserving any fragments.

# Arguments

  - `structure::`[`Structure`](@ref)
  - `rng::AbstractRNG`: Random number generator. If you need consistent results for testing
    purposes, pass a seeded pseudorandom number generator here, eg.
    [`Xoshiro(seed)`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.Xoshiro).
    Defaults to
    [`Random.default_rng()`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.default_rng).

# Example

```jldoctest; setup=:(using Random, MOLGA.Types; import MOLGA.GeneticAlgorithm.Mutation.switch_atom!)
julia> atoms = [
           Atom(8, [-0.702728547, 0.000000000, -3.588916226], 1, 1),
           Atom(6, [-0.702735805, -0.000007025, -2.382267226], 2, 2),
           Atom(1, [-0.702724918, 0.937663512, -1.787406767], 3, 3),
           Atom(1, [-0.702724918, -0.937656488, -1.787385685], 4, 3),
       ];

julia> structure = Structure(; atoms)
Structure with 4 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 0
   frag perm position
O     1    1 [-0.702728547         0.0                  -3.588916226        ]
C     2    2 [-0.702735805         -7.025e-6            -2.382267226        ]
H     3    3 [-0.702724918         0.937663512          -1.787406767        ]
H     3    4 [-0.702724918         -0.937656488         -1.787385685        ]

julia> switch_atom!(structure; rng=Xoshiro(123))

julia> structure
Structure with 4 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 1
   frag perm position
O     1    1 [-0.702724918         0.937663512          -1.787406767        ]
C     2    2 [-0.702735805         -7.025e-6            -2.382267226        ]
H     3    3 [-0.702728547         0.0                  -3.588916226        ]
H     3    4 [-0.702724918         -0.937656488         -1.787385685        ]
```
"""
function switch_atom!(structure::Structure; rng::AbstractRNG=Random.default_rng())
    num_atoms = length(structure.atoms)
    idx_1 = rand(rng, 1:num_atoms)

    for idx_2 in randperm(rng, num_atoms)
        # TODO: do we need a distance check here?
        if idx_1 !== idx_2 && structure.atoms[idx_1].element != structure.atoms[idx_2].element
            switch_atom_by_indices!(structure, idx_1, idx_2)
            return nothing
        end
    end

    @warn "Skipping switch_atom! because two different atoms cannot be found."
    return nothing
end

"""
    switch_atom_by_indices!(structure, index_1, index_2)

Exchange the position of the index-specified atoms from the provided structure in-place without
preserving any fragments.

# Arguments

  - `structure::`[`Structure`](@ref)
  - `index_1::Int`: Index of the first atom.
  - `index_2::Int`: Index of the second atom.

# Example

```jldoctest; setup=:(using MOLGA.Types; import MOLGA.GeneticAlgorithm.Mutation.switch_atom_by_indices!)
julia> atoms = [
           Atom(8, [-0.702728547, 0.000000000, -3.588916226], 1, 1),
           Atom(6, [-0.702735805, -0.000007025, -2.382267226], 2, 2),
           Atom(1, [-0.702724918, 0.937663512, -1.787406767], 3, 3),
           Atom(1, [-0.702724918, -0.937656488, -1.787385685], 4, 3),
       ];

julia> structure = Structure(; atoms)
Structure with 4 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 0
   frag perm position
O     1    1 [-0.702728547         0.0                  -3.588916226        ]
C     2    2 [-0.702735805         -7.025e-6            -2.382267226        ]
H     3    3 [-0.702724918         0.937663512          -1.787406767        ]
H     3    4 [-0.702724918         -0.937656488         -1.787385685        ]

julia> switch_atom_by_indices!(structure, 1, 4)

julia> structure
Structure with 4 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 1
   frag perm position
O     1    1 [-0.702724918         -0.937656488         -1.787385685        ]
C     2    2 [-0.702735805         -7.025e-6            -2.382267226        ]
H     3    3 [-0.702724918         0.937663512          -1.787406767        ]
H     3    4 [-0.702728547         0.0                  -3.588916226        ]
```
"""
function switch_atom_by_indices!(structure::Structure, index_1::Int, index_2::Int)
    temp_position = structure.atoms[index_1].position
    structure.atoms[index_1].position = structure.atoms[index_2].position
    structure.atoms[index_2].position = temp_position

    structure.age += 1

    return nothing
end

"""
    switch_atom_preserving!(structure, distance_thresholds [; rng])

Exchange the position of two randomly selected atoms and/or fragments from the provided structure
in-place while preserving fragments.

We loop over all possible combinations in a random order and try the switch. If it does't work, the
structure is not modified.

# Arguments

  - `structure::`[`Structure`](@ref)
  - `distance_thresholds::`[`DistanceThresholds`](@ref)
  - `rng::AbstractRNG`: Random number generator. If you need consistent results for testing
    purposes, pass a seeded pseudorandom number generator here, eg.
    [`Xoshiro(seed)`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.Xoshiro).
    Defaults to
    [`Random.default_rng()`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.default_rng).

# Example

```jldoctest; setup=:(using Random, MOLGA, MOLGA.Configuration; import MOLGA.GeneticAlgorithm.Mutation.switch_atom_preserving!)
julia> structure = deepcopy(MOLGA.Samples.DEMO_STRUCTURE_Ag)
Structure with 13 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 0
   frag perm position
O     2    2 [1.834176             -1.086997            -0.145691           ]
H     2    2 [2.642219             -0.527623            0.117024            ]
H     2    2 [2.082181             -1.613762            -0.930716           ]
O     2    3 [3.732413             0.511988             0.663892            ]
H     2    3 [4.173112             0.327927             1.518267            ]
H     2    3 [4.423184             0.920848             0.10458             ]
Ag    1    1 [1.2e-5               -0.000259            -0.210701           ]
O     2    4 [-1.833937            1.086716             -0.146848           ]
H     2    4 [-2.081857            1.612214             -0.932792           ]
H     2    4 [-2.64222             0.527829             0.116511            ]
O     2    5 [-3.732675            -0.510572            0.664366            ]
H     2    5 [-4.423943            -0.91958             0.105806            ]
H     2    5 [-4.173075            -0.324762            1.518501            ]

julia> switch_atom_preserving!(structure, DistanceThresholds(0.6, 5); rng=Xoshiro(123))

julia> structure
Structure with 13 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 0
   frag perm position
O     2    2 [-0.2479929999999999  0.5265059999999999   0.574324            ]
H     2    2 [0.5600499999999999   1.08588              0.837039            ]
H     2    2 [1.1999999999900979e-5 -0.0002590000000000092 -0.21070100000000003]
O     2    3 [3.732413             0.511988             0.663892            ]
H     2    3 [4.173112             0.327927             1.518267            ]
H     2    3 [4.423184             0.920848             0.10458             ]
Ag    1    1 [2.082181             -1.613762            -0.930716           ]
O     2    4 [-1.833937            1.086716             -0.146848           ]
H     2    4 [-2.081857            1.612214             -0.932792           ]
H     2    4 [-2.64222             0.527829             0.116511            ]
O     2    5 [-3.732675            -0.510572            0.664366            ]
H     2    5 [-4.423943            -0.91958             0.105806            ]
H     2    5 [-4.173075            -0.324762            1.518501            ]
```
"""
function switch_atom_preserving!(
    structure::Structure,
    distance_thresholds::DistanceThresholds;
    rng::AbstractRNG=Random.default_rng(),
)
    num_atoms = length(structure.atoms)

    # max number of combinations is Binomial(num_atoms, 2),
    # preferring looping through combinations over creating an array with all combinations
    for idx_1 in shuffle(rng, 1:num_atoms)
        for idx_2 in shuffle(rng, (idx_1 + 1):num_atoms)
            # we are not trying to exchange two atoms from the same fragment kind,
            # including an atom with itself
            if (
                idx_1 !== idx_2 &&
                structure.atoms[idx_1].fragment_kind != structure.atoms[idx_2].fragment_kind
            )
                if try_switch_atom_preserving!(structure, idx_1, idx_2, distance_thresholds)
                    return nothing
                end
            end
        end
    end

    @warn "Skipping switch_atom_preserving! because two different atoms cannot be found."
    return nothing
end

"""
    try_switch_atom_preserving!(structure, index_1, index_2, distance_thresholds)

Try to exchange the position of the atoms and/or fragments, specified by their indices, in-place
with preservation of fragments.

While building up the new atom list, distance checks are subsequently performed. Here, we only check
the lower distance threshold.

# Arguments

  - `structure::`[`Structure`](@ref)
  - `index_1::Int`: Index of the first atom.
  - `index_2::Int`: Index of the second atom.
  - `distance_thresholds::`[`DistanceThresholds`](@ref)

# Returns

  - `true` if the switch is successful,
  - `false` if a distance check fails.

# Example

```jldoctest; setup=:(using MOLGA,MOLGA.Configuration; import MOLGA.GeneticAlgorithm.Mutation.try_switch_atom_preserving!)
julia> structure = deepcopy(MOLGA.Samples.DEMO_STRUCTURE_Ag)
Structure with 13 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 0
   frag perm position
O     2    2 [1.834176             -1.086997            -0.145691           ]
H     2    2 [2.642219             -0.527623            0.117024            ]
H     2    2 [2.082181             -1.613762            -0.930716           ]
O     2    3 [3.732413             0.511988             0.663892            ]
H     2    3 [4.173112             0.327927             1.518267            ]
H     2    3 [4.423184             0.920848             0.10458             ]
Ag    1    1 [1.2e-5               -0.000259            -0.210701           ]
O     2    4 [-1.833937            1.086716             -0.146848           ]
H     2    4 [-2.081857            1.612214             -0.932792           ]
H     2    4 [-2.64222             0.527829             0.116511            ]
O     2    5 [-3.732675            -0.510572            0.664366            ]
H     2    5 [-4.423943            -0.91958             0.105806            ]
H     2    5 [-4.173075            -0.324762            1.518501            ]

julia> try_switch_atom_preserving!(structure, 2, 4, DistanceThresholds(0.6, 5))
false

julia> try_switch_atom_preserving!(structure, 10, 11, DistanceThresholds(0.6, 5))
false

julia> try_switch_atom_preserving!(structure, 10, 13, DistanceThresholds(1, 5))
false

julia> try_switch_atom_preserving!(structure, 2, 7, DistanceThresholds(1.6, 5))
false

julia> try_switch_atom_preserving!(structure, 2, 7, DistanceThresholds(0.6, 5))
true

julia> structure
Structure with 13 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 0
   frag perm position
O     2    2 [-0.8080309999999999  -0.559633            -0.47341599999999995]
H     2    2 [1.1999999999900979e-5 -0.0002590000000000092 -0.21070099999999997]
H     2    2 [-0.5600260000000001  -1.086398            -1.258441           ]
O     2    3 [3.732413             0.511988             0.663892            ]
H     2    3 [4.173112             0.327927             1.518267            ]
H     2    3 [4.423184             0.920848             0.10458             ]
Ag    1    1 [2.642219             -0.527623            0.11702399999999999 ]
O     2    4 [-1.833937            1.086716             -0.146848           ]
H     2    4 [-2.081857            1.612214             -0.932792           ]
H     2    4 [-2.64222             0.527829             0.116511            ]
O     2    5 [-3.732675            -0.510572            0.664366            ]
H     2    5 [-4.423943            -0.91958             0.105806            ]
H     2    5 [-4.173075            -0.324762            1.518501            ]
```
"""
function try_switch_atom_preserving!(
    structure::Structure, index_1::Int, index_2::Int, distance_thresholds::DistanceThresholds
)
    displacement = structure.atoms[index_1].position - structure.atoms[index_2].position

    permanent_idx_1 = structure.atoms[index_1].permanent_index
    permanent_idx_2 = structure.atoms[index_2].permanent_index

    new_atoms = Vector{Atom}()

    for atom in structure.atoms
        new_pos = atom.position

        if atom.permanent_index == permanent_idx_1
            new_pos -= displacement
        elseif atom.permanent_index == permanent_idx_2
            new_pos += displacement
        end

        # we only check the lower distance threshold here
        if check_distance(new_pos, new_atoms, distance_thresholds, true)
            push!(new_atoms, Atom(atom.element, new_pos, atom.permanent_index, atom.fragment_kind))
        else
            @debug(
                "Distance check within switch_atom_by_indices_preserving! failed.",
                structure,
                index_1,
                index_2,
            )
            return false
        end
    end

    structure.atoms = new_atoms

    return true
end

"""
    reposition_atom!(structure, box_size [; rng])

Move a randomly chosen atom from the provided structure in-place to another random position within
the cuboid box. This doesn't preserve any fragments.

# Arguments

  - `structure::`[`Structure`](@ref)
  - `box_size::`[`Vec`](@ref)
  - `rng::AbstractRNG`: Random number generator. If you need consistent results for testing
    purposes, pass a seeded pseudorandom number generator here, eg.
    [`Xoshiro(seed)`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.Xoshiro).
    Defaults to
    [`Random.default_rng()`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.default_rng).

# Example

```jldoctest; setup=:(using Random, MOLGA.Types; import MOLGA.GeneticAlgorithm.Mutation.reposition_atom!)
julia> atoms = [
           Atom(8, [-0.702728547, 0.000000000, -3.588916226], 1, 1),
           Atom(6, [-0.702735805, -0.000007025, -2.382267226], 2, 2),
           Atom(1, [-0.702724918, 0.937663512, -1.787406767], 3, 3),
           Atom(1, [-0.702724918, -0.937656488, -1.787385685], 4, 3),
       ];

julia> structure = Structure(; atoms)
Structure with 4 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 0
   frag perm position
O     1    1 [-0.702728547         0.0                  -3.588916226        ]
C     2    2 [-0.702735805         -7.025e-6            -2.382267226        ]
H     3    3 [-0.702724918         0.937663512          -1.787406767        ]
H     3    4 [-0.702724918         -0.937656488         -1.787385685        ]

julia> reposition_atom!(structure, Vec(8, 8, 8); rng=Xoshiro(123))

julia> structure
Structure with 4 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 0
   frag perm position
O     1    1 [-0.702728547         0.0                  -3.588916226        ]
C     2    2 [-0.702735805         -7.025e-6            -2.382267226        ]
H     3    3 [0.6944540596267874   3.1270295847422487   -2.472746407793897  ]
H     3    4 [-0.702724918         -0.937656488         -1.787385685        ]
```
"""
function reposition_atom!(
    structure::Structure, box_size::Vec; rng::AbstractRNG=Random.default_rng()
)
    num_atoms = length(structure.atoms)
    atom_index = rand(rng, 1:num_atoms)

    # TODO: do we need a distance check here as well?
    structure.atoms[atom_index].position = random_position(box_size; rng)

    return nothing
end

"""
    reposition_atom_preserving!(structure, box_size, distance_thresholds [; rng])

Move a randomly chosen atom or fragment from the provided structure in-place to another random
position and preserve any fragments.

All atoms are tried out in a randomized order. If there is no successful reposition, the structure
is not modified.

# Arguments

  - `structure::`[`Structure`](@ref)
  - `box_size::`[`Vec`](@ref)
  - `distance_thresholds::`[`DistanceThresholds`](@ref)
  - `rng::AbstractRNG`: Random number generator. If you need consistent results for testing
    purposes, pass a seeded pseudorandom number generator here, eg.
    [`Xoshiro(seed)`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.Xoshiro).
    Defaults to
    [`Random.default_rng()`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.default_rng).

# Example

```jldoctest; setup=:(using Random, MOLGA, MOLGA.Types, MOLGA.Configuration; import MOLGA.GeneticAlgorithm.Mutation.reposition_atom_preserving!)
julia> structure = deepcopy(MOLGA.Samples.DEMO_STRUCTURE_Ag)
Structure with 13 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 0
   frag perm position
O     2    2 [1.834176             -1.086997            -0.145691           ]
H     2    2 [2.642219             -0.527623            0.117024            ]
H     2    2 [2.082181             -1.613762            -0.930716           ]
O     2    3 [3.732413             0.511988             0.663892            ]
H     2    3 [4.173112             0.327927             1.518267            ]
H     2    3 [4.423184             0.920848             0.10458             ]
Ag    1    1 [1.2e-5               -0.000259            -0.210701           ]
O     2    4 [-1.833937            1.086716             -0.146848           ]
H     2    4 [-2.081857            1.612214             -0.932792           ]
H     2    4 [-2.64222             0.527829             0.116511            ]
O     2    5 [-3.732675            -0.510572            0.664366            ]
H     2    5 [-4.423943            -0.91958             0.105806            ]
H     2    5 [-4.173075            -0.324762            1.518501            ]

julia> thresholds = DistanceThresholds(0.6, 5);
       reposition_atom_preserving!(structure, Vec(8, 8, 8), thresholds; rng=Xoshiro(123))

julia> structure
Structure with 13 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 0
   frag perm position
O     2    2 [1.834176             -1.086997            -0.145691           ]
H     2    2 [2.642219             -0.527623            0.117024            ]
H     2    2 [2.082181             -1.613762            -0.930716           ]
O     2    3 [0.14825493974911375  1.4998088298154295   -2.4914429619616025 ]
H     2    3 [0.5889539397491133   1.3157478298154297   -1.6370679619616026 ]
H     2    3 [0.8390259397491135   1.9086688298154297   -3.050754961961603  ]
Ag    1    1 [1.2e-5               -0.000259            -0.210701           ]
O     2    4 [-1.833937            1.086716             -0.146848           ]
H     2    4 [-2.081857            1.612214             -0.932792           ]
H     2    4 [-2.64222             0.527829             0.116511            ]
O     2    5 [-3.732675            -0.510572            0.664366            ]
H     2    5 [-4.423943            -0.91958             0.105806            ]
H     2    5 [-4.173075            -0.324762            1.518501            ]
```
"""
function reposition_atom_preserving!(
    structure::Structure,
    box_size::Vec,
    distance_thresholds::DistanceThresholds;
    rng::AbstractRNG=Random.default_rng(),
)
    num_atoms = length(structure.atoms)

    for atom_idx in randperm(rng, num_atoms)
        rand_pos = random_position(box_size; rng)

        if try_reposition_atom_preserving!(structure, atom_idx, rand_pos, distance_thresholds)
            return nothing
        end
    end

    @warn "Skipping reposition_atom_preserving! because a suitable atom cannot be found."
    return nothing
end

"""
    try_reposition_atom_preserving!(structure, atom_index, new_position, distance_thresholds)

Try to move the specified atom or fragment from the provided structure to the provided position,
preserving fragments.

While building up the new atom list, distance checks are subsequently performed. Here, we only check
the lower distance threshold.

# Arguments

  - `structure::`[`Structure`](@ref)
  - `atom_index::Int`: Index of the atom to be moved.
  - `new_position::`[`Vec`](@ref)
  - `distance_thresholds::`[`DistanceThresholds`](@ref)

# Returns

  - `true` if the reposition is successful,
  - `false` if a distance check fails.

# Example

```jldoctest; setup=:(using MOLGA, MOLGA.Types, MOLGA.Configuration; import MOLGA.GeneticAlgorithm.Mutation.try_reposition_atom_preserving!)
julia> structure = deepcopy(MOLGA.Samples.DEMO_STRUCTURE_Ag)
Structure with 13 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 0
   frag perm position
O     2    2 [1.834176             -1.086997            -0.145691           ]
H     2    2 [2.642219             -0.527623            0.117024            ]
H     2    2 [2.082181             -1.613762            -0.930716           ]
O     2    3 [3.732413             0.511988             0.663892            ]
H     2    3 [4.173112             0.327927             1.518267            ]
H     2    3 [4.423184             0.920848             0.10458             ]
Ag    1    1 [1.2e-5               -0.000259            -0.210701           ]
O     2    4 [-1.833937            1.086716             -0.146848           ]
H     2    4 [-2.081857            1.612214             -0.932792           ]
H     2    4 [-2.64222             0.527829             0.116511            ]
O     2    5 [-3.732675            -0.510572            0.664366            ]
H     2    5 [-4.423943            -0.91958             0.105806            ]
H     2    5 [-4.173075            -0.324762            1.518501            ]

julia> thresholds = DistanceThresholds(0.6, 5);

julia> try_reposition_atom_preserving!(structure, 3, Vec(2.7924, -0.8167, -0.1357), thresholds)
false

julia> try_reposition_atom_preserving!(structure, 7, Vec(-1.3685, 0.8138, -0.2107), thresholds)
false

julia> try_reposition_atom_preserving!(structure, 11, Vec(-1.5886, -0.9054, 0.6644), thresholds)
true

julia> structure
Structure with 13 atoms
Energy: 0.0
Nuclear repulsion energy: 0.0
Age: 0
   frag perm position
O     2    2 [1.834176             -1.086997            -0.145691           ]
H     2    2 [2.642219             -0.527623            0.117024            ]
H     2    2 [2.082181             -1.613762            -0.930716           ]
O     2    3 [3.732413             0.511988             0.663892            ]
H     2    3 [4.173112             0.327927             1.518267            ]
H     2    3 [4.423184             0.920848             0.10458             ]
Ag    1    1 [1.2e-5               -0.000259            -0.210701           ]
O     2    4 [-1.833937            1.086716             -0.146848           ]
H     2    4 [-2.081857            1.612214             -0.932792           ]
H     2    4 [-2.64222             0.527829             0.116511            ]
O     2    5 [-1.5886              -0.9054              0.6644              ]
H     2    5 [-2.2798680000000004  -1.3144079999999998  0.10583999999999998 ]
H     2    5 [-2.029               -0.71959             1.518535            ]
```
"""
function try_reposition_atom_preserving!(
    structure::Structure,
    atom_index::Int,
    new_position::Vec,
    distance_thresholds::DistanceThresholds,
)
    displacement = structure.atoms[atom_index].position - new_position
    permanent_idx = structure.atoms[atom_index].permanent_index

    new_atoms = Vector{Atom}()

    for atom in structure.atoms
        pos = atom.position

        if atom.permanent_index == permanent_idx
            pos -= displacement
        end

        # we only check the lower distance threshold here
        if check_distance(pos, new_atoms, distance_thresholds, true)
            push!(new_atoms, Atom(atom.element, pos, atom.permanent_index, atom.fragment_kind))
        else
            @debug "Distance check within reposition_atom_preserving! failed" structure atom_index
            return false
        end
    end

    structure.atoms = new_atoms

    return true
end

end
