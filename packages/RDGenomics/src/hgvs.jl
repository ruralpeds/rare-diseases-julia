# Minimal HGVS parser.
#
# We support the small subset of HGVS coding (.c) and protein (.p) syntax
# that appears in ClinVar's `Name` field for SNVs and short indels — enough
# to round-trip and key into a normalized representation. We do NOT
# implement the full HGVS grammar (Mutalyzer covers that), and we deliberately
# refuse exotic forms rather than guess.
#
# Examples we handle:
#   NM_000277.3:c.1241A>G
#   NM_000277.3:c.842+1G>A
#   NM_000277.3:c.143_144del
#   NM_000277.3:c.143_144insGT
#   NM_000277.3:c.143_144delinsAT
#   NP_000268.1:p.Arg408Trp
#   NP_000268.1:p.R408W
#   NP_000268.1:p.Arg408*
#   NP_000268.1:p.Arg408Ter
#   NP_000268.1:p.Arg408del
#   NP_000268.1:p.Arg408fs

"""
    HgvsCoding

Parsed `c.` (coding-sequence) HGVS expression. `kind` is one of
`:substitution`, `:deletion`, `:insertion`, `:delins`, `:duplication`.
"""
struct HgvsCoding
    reference::String          # e.g. "NM_000277.3"
    kind::Symbol
    start::Int
    stop::Int                  # == start for SNV/dup-of-1
    intronic_offset::Int       # 0 for fully exonic
    ref::String                # may be "" for ins
    alt::String                # may be "" for del
end

"""
    HgvsProtein

Parsed `p.` (protein) HGVS expression. `kind` is one of `:missense`,
`:nonsense`, `:synonymous`, `:deletion`, `:frameshift`, `:insertion`,
`:delins`, `:initiator_loss`.
"""
struct HgvsProtein
    reference::String
    kind::Symbol
    position::Int
    ref_aa::String       # one-letter, e.g. "R"
    alt_aa::String       # one-letter, e.g. "W"; "*" for stop, "" for del
end

const _AA3_TO_1 = Dict(
    "Ala"=>"A","Arg"=>"R","Asn"=>"N","Asp"=>"D","Cys"=>"C",
    "Gln"=>"Q","Glu"=>"E","Gly"=>"G","His"=>"H","Ile"=>"I",
    "Leu"=>"L","Lys"=>"K","Met"=>"M","Phe"=>"F","Pro"=>"P",
    "Ser"=>"S","Thr"=>"T","Trp"=>"W","Tyr"=>"Y","Val"=>"V",
    "Sec"=>"U","Pyl"=>"O","Ter"=>"*","*"=>"*",
)
const _AA1_VALID = Set("ACDEFGHIKLMNPQRSTVWY*U")

"""
    parse_hgvs_c(s) -> HgvsCoding

Parse a coding HGVS string. Throws `ArgumentError` if the string does not
match the supported subset.
"""
function parse_hgvs_c(s::AbstractString)
    parts = split(s, ":"; limit=2)
    length(parts) == 2 || throw(ArgumentError("HGVS: missing ':' in '$s'"))
    ref = String(parts[1])
    body = String(parts[2])
    startswith(body, "c.") ||
        throw(ArgumentError("HGVS: '$s' is not a c. expression"))
    expr = body[3:end]

    # Substitution: 1241A>G   or   842+1G>A   or   842-2C>T
    m = match(r"^(\d+)([+-]\d+)?([ACGT])>([ACGT])$", expr)
    if m !== nothing
        pos = parse(Int, m.captures[1])
        off = m.captures[2] === nothing ? 0 : parse(Int, m.captures[2])
        return HgvsCoding(ref, :substitution, pos, pos, off,
                          String(m.captures[3]), String(m.captures[4]))
    end

    # delins (must come before del/ins because it contains both substrings)
    m = match(r"^(\d+)(?:_(\d+))?delins([ACGT]+)$", expr)
    if m !== nothing
        a = parse(Int, m.captures[1])
        b = m.captures[2] === nothing ? a : parse(Int, m.captures[2])
        return HgvsCoding(ref, :delins, a, b, 0, "", String(m.captures[3]))
    end

    # Deletion: 143del  or  143_144del  or  143_144delAG
    m = match(r"^(\d+)(?:_(\d+))?del([ACGT]*)$", expr)
    if m !== nothing
        a = parse(Int, m.captures[1])
        b = m.captures[2] === nothing ? a : parse(Int, m.captures[2])
        return HgvsCoding(ref, :deletion, a, b, 0, String(m.captures[3]), "")
    end

    # Insertion: 143_144insGT — by HGVS, ins requires a flanking pair
    m = match(r"^(\d+)_(\d+)ins([ACGT]+)$", expr)
    if m !== nothing
        a = parse(Int, m.captures[1])
        b = parse(Int, m.captures[2])
        b == a + 1 ||
            throw(ArgumentError("HGVS ins: flanks must be adjacent in '$s'"))
        return HgvsCoding(ref, :insertion, a, b, 0, "", String(m.captures[3]))
    end

    # Duplication: 143dup or 143_144dup
    m = match(r"^(\d+)(?:_(\d+))?dup([ACGT]*)$", expr)
    if m !== nothing
        a = parse(Int, m.captures[1])
        b = m.captures[2] === nothing ? a : parse(Int, m.captures[2])
        return HgvsCoding(ref, :duplication, a, b, 0, String(m.captures[3]), "")
    end

    throw(ArgumentError("HGVS: unsupported c. expression '$s'"))
end

"""
    parse_hgvs_p(s) -> HgvsProtein

Parse a protein HGVS string. Accepts both 1-letter and 3-letter amino-acid
codes. Throws `ArgumentError` on unsupported forms.
"""
function parse_hgvs_p(s::AbstractString)
    parts = split(s, ":"; limit=2)
    length(parts) == 2 || throw(ArgumentError("HGVS: missing ':' in '$s'"))
    ref = String(parts[1])
    body = String(parts[2])
    startswith(body, "p.") ||
        throw(ArgumentError("HGVS: '$s' is not a p. expression"))
    expr = strip(body[3:end])
    # Strip optional surrounding parentheses (predicted)
    if startswith(expr, "(") && endswith(expr, ")")
        expr = expr[2:end-1]
    end

    # 3-letter substitution: Arg408Trp / Arg408Ter / Arg408*
    m = match(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)$", expr)
    if m !== nothing
        ref_aa = _aa1(String(m.captures[1]))
        pos = parse(Int, m.captures[2])
        raw = String(m.captures[3])
        alt_aa = raw == "*" ? "*" : _aa1(raw)
        kind = alt_aa == "*" ? :nonsense :
               ref_aa == alt_aa ? :synonymous : :missense
        return HgvsProtein(ref, kind, pos, ref_aa, alt_aa)
    end

    # 1-letter substitution: R408W
    m = match(r"^([ACDEFGHIKLMNPQRSTVWY\*U])(\d+)([ACDEFGHIKLMNPQRSTVWY\*U])$", expr)
    if m !== nothing
        ref_aa = String(m.captures[1])
        pos = parse(Int, m.captures[2])
        alt_aa = String(m.captures[3])
        kind = alt_aa == "*" ? :nonsense :
               ref_aa == alt_aa ? :synonymous : :missense
        return HgvsProtein(ref, kind, pos, ref_aa, alt_aa)
    end

    # del: Arg408del or R408del
    m = match(r"^([A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])(\d+)del$", expr)
    if m !== nothing
        ref_aa = _aa1(String(m.captures[1]))
        pos = parse(Int, m.captures[2])
        return HgvsProtein(ref, :deletion, pos, ref_aa, "")
    end

    # frameshift: Arg408fs or Arg408AlafsTer5
    m = match(r"^([A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])(\d+).*fs.*$", expr)
    if m !== nothing
        ref_aa = _aa1(String(m.captures[1]))
        pos = parse(Int, m.captures[2])
        return HgvsProtein(ref, :frameshift, pos, ref_aa, "")
    end

    # initiator loss: Met1?
    if occursin(r"^(Met1|M1)\?$", expr)
        return HgvsProtein(ref, :initiator_loss, 1, "M", "")
    end

    throw(ArgumentError("HGVS: unsupported p. expression '$s'"))
end

function _aa1(code::AbstractString)
    if length(code) == 1
        c = uppercase(code)
        first(c) in _AA1_VALID && return String(c)
        throw(ArgumentError("invalid 1-letter amino acid '$code'"))
    end
    haskey(_AA3_TO_1, code) || throw(ArgumentError("invalid 3-letter amino acid '$code'"))
    return _AA3_TO_1[code]
end
