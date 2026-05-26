"""
    RDImmunology

Comprehensive computational immunology ecosystem base types.
Provides abstract types and standard representations for immune components
(cells, cytokines, receptors, pathogens, epitopes) to facilitate mapping to
reaction networks (`Catalyst.jl`) and agent-based models (`Agents.jl`).
"""
module RDImmunology

using RareDiseaseCore

export
    AbstractImmuneComponent,
    AbstractImmuneCell,
    AbstractCytokine,
    AbstractReceptor,
    AbstractPathogen,
    Epitope,
    Macrophage,
    TCell,
    BCell

"""
    AbstractImmuneComponent

Base type for all modeled entities in the immunology ecosystem.
"""
abstract type AbstractImmuneComponent end

"""
    AbstractImmuneCell <: AbstractImmuneComponent

Base type for all immune cells (e.g., Macrophage, TCell, BCell).
"""
abstract type AbstractImmuneCell <: AbstractImmuneComponent end

"""
    AbstractCytokine <: AbstractImmuneComponent

Base type for signaling molecules (e.g., IL-1, TNF-alpha).
"""
abstract type AbstractCytokine <: AbstractImmuneComponent end

"""
    AbstractReceptor <: AbstractImmuneComponent

Base type for cell-surface receptors (e.g., TCR, BCR, TLR).
"""
abstract type AbstractReceptor <: AbstractImmuneComponent end

"""
    AbstractPathogen <: AbstractImmuneComponent

Base type for invading entities (viruses, bacteria, etc.).
"""
abstract type AbstractPathogen <: AbstractImmuneComponent end

"""
    Epitope

Represents a specific antigenic determinant recognized by the immune system.
"""
struct Epitope
    sequence::String
    source_organism::String
end

# Example cell subtypes
struct Macrophage <: AbstractImmuneCell end
struct TCell <: AbstractImmuneCell end
struct BCell <: AbstractImmuneCell end

end # module
