module Configuration

using Base.Filesystem
using TOML
using DocStringExtensions
using PeriodicTable

using ..Types
using ..IOm

export ConfigurationObject,
    DistanceThresholds,
    MutationConfiguration,
    InitialPopulationConfiguration,
    ConfigAtom,
    ConfigMolecule,
    AtomConfig

"""
    DistanceThresholds(min, max)

Defines the thresholds for distances between atoms.
$(TYPEDFIELDS)
The values are taken into account at several points, including in the generation of the initial 
population and in mutation and recombination. They are given in angstroms.
"""
struct DistanceThresholds
    """Minimum distance that two atoms need."""
    min::Float64

    """Maximum distance that two atoms need."""
    max::Float64
end

"""
    ConfigAtom(num, element)

Atom to define initial atom configuration.
$(TYPEDFIELDS)
"""
struct ConfigAtom
    """Number of atoms of this kind."""
    quantity::Int

    """Chemical element denoted as atomic number."""
    element::Int
end

"""
    ConfigMolecule(num, atoms)

Molecule to define initial atom configuration.
$(TYPEDFIELDS)
"""
struct ConfigMolecule
    """Number of molecules of this kind."""
    quantity::Int

    """Array of [atoms](@ref Types.BaseAtom) that the molecule consists of."""
    atoms::Vector{BaseAtom}
end

"""
The atom configuration consists of single [atoms](@ref ConfigAtom) and
[molecules](@ref ConfigMolecule).
"""
const AtomConfig = Array{Union{ConfigAtom,ConfigMolecule}}

"""
Defines the configuration of the initial population.
$(TYPEDFIELDS)
"""
struct InitialPopulationConfiguration
    """The number of initial structures to be created."""
    num_structures::Int

    """The three dimensions of the periodic boundary condition box (in angstroms)."""
    box_size::Types.Vec

    """See [`AtomConfig`](@ref)."""
    atom_config::AtomConfig
end

"""
Configure the mutation parameters.
$(TYPEDFIELDS)
"""
struct MutationConfiguration
    """
    Probability in the interval [0,1] with which each individual structure mutates. Each mutation
    function is applied to the structure with this probability.
    """
    probability::Float64

    """
    If True, when a mutation happens, the structure before and after the mutation is kept; otherwise
    only the mutated structure is kept.
    """
    keep_non_mutated::Bool
end

"""
Stores the entire configuration parameters.
$(FIELDS)
"""
struct ConfigurationObject
    """See [`InitialPopulationConfiguration`](@ref)."""
    initial_population::InitialPopulationConfiguration

    """Charge of the structure."""
    charge::Int

    """Spin multiplicity of the structure."""
    spin_multiplicity::Int

    """See [`DistanceThresholds`](@ref)."""
    distance_thresholds::DistanceThresholds

    """True if fragments should be preserved during recombination and mutation."""
    preserve_fragments::Bool

    """See [`MutationConfiguration`](@ref)."""
    mutation::MutationConfiguration
end

"""
    create_atom_config(atoms::Dict{String,Int}, molecules::Dict{String,Int}, config_path::String)

Creates [`AtomConfig`](@ref) from the TOML raw values (atoms, molecules). Elements are internally
represented via their atomic number. The molecules are loaded and parsed from their xyz-files.

The `config_path` is the folder where the config file is stored. The molecule file paths are
expected to be relative to this folder.
"""
function create_atom_config(
    atoms::Dict{String,Int}, molecules::Dict{String,Int}, config_path::String
)
    return [
        [ConfigAtom(num, elements[Symbol(element_symb)].number) for (element_symb, num) in atoms]
        [
            ConfigMolecule(num, load_xyz(Filesystem.joinpath(config_path, molecule_path))) for
            (molecule_path, num) in molecules
        ]
    ]
end

"""
    load_config(configuration_file)

Loads the configuration parameters from the specified TOML file.
"""
function load_config(configuration_file::String)
    @debug "Read config from $(configuration_file)"

    parsed_config = TOML.tryparsefile(configuration_file)

    if isa(parsed_config, TOML.ParserError)
        error("Error parsing configuration file $(configuration_file)")
    end

    initial_population = InitialPopulationConfiguration(
        parsed_config["initial_population"]["num_structures"],
        parsed_config["initial_population"]["box_size"],
        create_atom_config(
            Dict{String,Int}(parsed_config["initial_population"]["atoms"]),
            Dict{String,Int}(parsed_config["initial_population"]["molecules"]),
            Filesystem.dirname(configuration_file),
        ),
    )

    return ConfigurationObject(
        initial_population,
        parsed_config["charge"],
        parsed_config["spin_multiplicity"],
        DistanceThresholds(
            parsed_config["distance_thresholds"]["min"], parsed_config["distance_thresholds"]["max"]
        ),
        parsed_config["preserve_fragments"],
        MutationConfiguration(
            parsed_config["mutation"]["probability"], parsed_config["mutation"]["keep_non_mutated"]
        ),
    )
end

end
