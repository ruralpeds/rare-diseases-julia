# OBO 1.2 / 1.4 stanza parser.
#
# We parse the subset of OBO actually used by HPO / MONDO / Orphanet / DO:
#   * `[Typedef]` and `[Term]` stanzas
#   * repeated tags collected into a vector
#   * comments (` ! ...`) stripped from values
#   * escaped colons in IDs not supported (not used by these ontologies)
#
# Streaming line-based parse; allocates one `OBOStanza` per term.

struct OBOStanza
    kind::String                          # "Term" | "Typedef" | "Header"
    tags::Dict{String,Vector{String}}
end
OBOStanza(kind) = OBOStanza(kind, Dict{String,Vector{String}}())

"""
    parse_obo(path) -> Vector{OBOStanza}

Parse an OBO file into stanzas. The header is returned as the first
stanza with `kind == "Header"`; subsequent stanzas follow the order in
the file.
"""
function parse_obo(path::AbstractString)
    out = OBOStanza[]
    header = OBOStanza("Header")
    current = header
    push!(out, header)

    open(path, "r") do io
        for raw_line in eachline(io)
            line = strip(raw_line)
            isempty(line) && continue
            startswith(line, "!") && continue   # comment line

            if startswith(line, "[") && endswith(line, "]")
                kind = String(line[2:end-1])
                current = OBOStanza(kind)
                push!(out, current)
                continue
            end

            # tag: value [! comment]
            sep = findfirst(':', line)
            sep === nothing && continue
            tag = String(strip(SubString(line, 1, sep - 1)))
            value = SubString(line, sep + 1)

            # Strip an unescaped trailing ` ! comment`. OBO escapes ! with \!.
            value = _strip_trailing_comment(value)
            value = strip(value)
            isempty(tag) && continue
            push!(get!(current.tags, tag, String[]), String(value))
        end
    end
    return out
end

function _strip_trailing_comment(s::AbstractString)
    # Find the first unescaped '!'
    i = 1
    chars = collect(s)
    n = length(chars)
    while i <= n
        c = chars[i]
        if c == '\\' && i < n
            i += 2
            continue
        end
        if c == '!'
            return String(strip(SubString(s, 1, i - 1)))
        end
        i += 1
    end
    return String(s)
end
