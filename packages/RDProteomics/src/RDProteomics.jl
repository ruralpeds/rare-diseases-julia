"""
    RDProteomics

UniProt, AlphaFold, InterPro, Pfam. Maps gene-coordinate variants onto
UniProt residues and 3-D structures; exposes per-residue features
(pLDDT, secondary structure, domain membership, distance to active site).

Built on `BioStructures.jl` for PDB/mmCIF parsing and `FASTX.jl` for
sequence I/O.

Phase 5 of the build plan.

# source: UniProt
# source: AlphaFoldDB
# source: PDB
# source: InterPro
# source: Pfam
"""
module RDProteomics

using BioStructures
using FASTX
using RareDiseaseCore

export
    ResidueFeature,
    load_uniprot, fetch_alphafold, read_structure,
    map_variant_to_residue, residue_features

"""
    ResidueFeature

Per-residue feature row keyed by `(accession, position)`.
"""
struct ResidueFeature
    accession::UniProtAcc
    position::Int
    plddt::Union{Nothing,Float64}
    secondary::Union{Nothing,Symbol}    # :helix | :sheet | :loop
    domains::Vector{String}
    distance_to_active_site::Union{Nothing,Float64}
end

"""
    read_structure(path; format=BioStructures.PDBFormat) -> MolecularStructure

Thin pass-through to `BioStructures.read`. Use this entry point so callers
don't need to import BioStructures explicitly.
"""
read_structure(path::AbstractString; format=PDBFormat) =
    BioStructures.read(path, format)

load_uniprot(::AbstractString)             = error("load_uniprot not yet implemented (Phase 5)")
fetch_alphafold(::UniProtAcc; kw...)       = error("fetch_alphafold not yet implemented (Phase 5)")
map_variant_to_residue(::Variant)          = error("map_variant_to_residue not yet implemented (Phase 5)")
residue_features(::UniProtAcc, ::Int)      = error("residue_features not yet implemented (Phase 5)")

end # module
