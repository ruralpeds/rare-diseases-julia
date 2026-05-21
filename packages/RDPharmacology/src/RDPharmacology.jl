"""
    RDPharmacology

ChEMBL, DrugBank-open, PubChem, BindingDB, RxNorm, FDA Orange Book, FDA OOPD,
DailyMed. Drug-target bioactivity, approval/orphan status, and structured
product labels.

Phase 7 of the build plan.

# source: ChEMBL
# source: DrugBank-open
# source: PubChem
# source: BindingDB
# source: RxNorm
# source: FDA-OrangeBook
# source: FDA-OOPD
# source: DailyMed
"""
module RDPharmacology

using RareDiseaseCore

export
    Bioactivity,
    drugs_for_target, orphan_designations_for, indications_for

"""
    Bioactivity

Drug–target affinity record. `value_um` is concentration in μM; `kind`
distinguishes IC50/Ki/Kd/EC50.
"""
struct Bioactivity
    drug::ChemblId
    target::UniProtAcc
    kind::Symbol           # :IC50 | :Ki | :Kd | :EC50
    value_um::Float64
    source::String
end

drugs_for_target(::UniProtAcc; max_um::Float64=1.0) =
    error("drugs_for_target not yet implemented (Phase 7)")
orphan_designations_for(::MondoId) =
    error("orphan_designations_for not yet implemented (Phase 7)")
indications_for(::Union{RxCui,ChemblId}) =
    error("indications_for not yet implemented (Phase 7)")

end # module
