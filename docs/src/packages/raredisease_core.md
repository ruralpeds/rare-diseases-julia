# RareDiseaseCore

```@meta
CurrentModule = RareDiseaseCore
```

Core identifier types, base records, and provenance utilities shared
across every package in the monorepo.

## Identifier types

```@docs
HPOId
MondoId
OrphaId
OmimId
DoId
MeshId
HgncId
EnsemblGeneId
EnsemblTranscriptId
RefSeqId
UniProtAcc
ClinVarId
DbSnpId
PubMedId
PmcId
RxCui
ChemblId
PubChemCid
ReactomeId
PfamId
InterProId
```

## Provenance

```@docs
DataProvenance
cite
sha256_file
```

## Base records

```@docs
Phenotype
Disease
Gene
Variant
Protein
Drug
Trial
Pathway
```
