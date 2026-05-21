# BioModels loader.
#
# Curated disease-relevant SBML models live in EBI BioModels. We use
# SBMLToolkit.jl to convert SBML into a ModelingToolkit ReactionSystem
# that plays nicely with the Catalyst-based stack.
#
# source: BioModels

"""
    load_biomodel(path) -> ReactionSystem

Parse a BioModels SBML file into a `ModelingToolkit.ReactionSystem`
ready to feed into `ODEProblem`. Equivalent to
`SBMLToolkit.ReactionSystem(SBML.readSBML(path))`.
"""
function load_biomodel(path::AbstractString)
    isfile(path) || throw(ArgumentError("no such file: $path"))
    model = SBML.readSBML(path)
    return SBMLToolkit.ReactionSystem(model)
end
