function fg_abelian_group_string(value::Vector{Int64})
    isempty(value) && return "0"

    summands = String[]
    for n in value
        push!(summands, n == 0 ? "Z" : "Z/$n")
    end

    return join(summands, " + ")
end

function fg_abelian_group_tex(value::Vector{Int64})
    isempty(value) && return "0"

    summands = String[]
    for n in value
        push!(summands, n == 0 ? "\\mathbb{Z}" : "\\mathbb{Z}/$n")
    end

    return join(summands, " \\oplus ")
end

function subgroup_label_tex(label::AbstractString)
    label == "1" && return "1"
    return replace(label, r"([A-Za-z]+)([0-9]+)" => s"\1_{\2}")
end

matrix_rows(M::AbstractMatrix{<:Integer}) =
    [[Int(M[i, j]) for j in axes(M, 2)] for i in axes(M, 1)]

function gap_string_or_missing(f, fallback::String)
    try
        return string(f())
    catch
        return fallback
    end
end

function subgroup_structure_label(H::GapObj, i::Int)
    structure = gap_string_or_missing(() -> GAP.Globals.StructureDescription(H), "H$i")
    return structure == "fail" ? "H$i" : structure
end

function unique_subgroup_labels(subgroups::Vector{GapObj})
    base_labels = [subgroup_structure_label(H, i) for (i, H) in enumerate(subgroups)]
    counts = Dict{String,Int}()

    for label in base_labels
        counts[label] = get(counts, label, 0) + 1
    end

    seen = Dict{String,Int}()
    labels = String[]

    for label in base_labels
        if counts[label] == 1
            push!(labels, label)
        else
            seen[label] = get(seen, label, 0) + 1
            push!(labels, "$(label)_$(seen[label])")
        end
    end

    return labels
end

function subgroup_index(M::MackeyFunctor, H::GapObj)
    for (i, K) in enumerate(M.subgroups)
        K === H && return i
    end

    throw(ArgumentError("Subgroup is not stored in this Mackey functor"))
end

function matrix_data_entry(M::MackeyFunctor, H::GapObj, K::GapObj, matrix)
    return Dict{String,Any}(
        "source" => subgroup_index(M, H),
        "target" => subgroup_index(M, K),
        "matrix" => matrix_rows(matrix),
    )
end

function visualizer_data(M::MackeyFunctor)
    labels = unique_subgroup_labels(M.subgroups)

    nodes = Dict{String,Any}[]
    for (i, H) in enumerate(M.subgroups)
        id_group = gap_string_or_missing(() -> GAP.Globals.IdGroup(H), "unavailable")
        generators = gap_string_or_missing(() -> GAP.Globals.GeneratorsOfGroup(H), "[]")
        gap_name = string(H)

        push!(nodes, Dict{String,Any}(
            "id" => i,
            "label" => labels[i],
            "label_tex" => subgroup_label_tex(labels[i]),
            "value" => fg_abelian_group_string(M.values[H]),
            "value_tex" => fg_abelian_group_tex(M.values[H]),
            "raw_value" => copy(M.values[H]),
            "gap_name" => gap_name,
            "orbit" => "G/$(labels[i])",
            "orbit_tex" => "G/$(subgroup_label_tex(labels[i]))",
            "size" => Int(GAP.Globals.Size(H)),
            "id_group" => id_group,
            "generators" => generators,
        ))
    end

    covers = [
        Dict{String,Any}(
            "source" => subgroup_index(M, H),
            "target" => subgroup_index(M, K),
        )
        for (H, K) in M.covers
    ]

    restrictions = [
        matrix_data_entry(M, K, H, matrix)
        for ((H, K), matrix) in M.restrictions
    ]

    transfers = [
        matrix_data_entry(M, H, K, matrix)
        for ((H, K), matrix) in M.transfers
    ]

    conjugations = Dict{String,Any}[]
    generators = collect(GAP.Globals.GeneratorsOfGroup(M.G))

    for ((H, g), matrix) in M.conjugations
        any(gen -> gen == g, generators) || continue

        K = canonical_subgroup(M.subgroups, conjugate_subgroup(H, g))
        push!(conjugations, Dict{String,Any}(
            "source" => subgroup_index(M, H),
            "target" => subgroup_index(M, K),
            "element" => string(g),
            "matrix" => matrix_rows(matrix),
        ))
    end

    return Dict{String,Any}(
        "coefficient_ring" => string(M.coefficient_ring),
        "nodes" => nodes,
        "covers" => covers,
        "restrictions" => restrictions,
        "transfers" => transfers,
        "conjugations" => conjugations,
    )
end

function json_escape(s::AbstractString)
    result = IOBuffer()
    for c in s
        if c == '\\'
            print(result, "\\\\")
        elseif c == '"'
            print(result, "\\\"")
        elseif c == '\n'
            print(result, "\\n")
        elseif c == '\r'
            print(result, "\\r")
        elseif c == '\t'
            print(result, "\\t")
        else
            print(result, c)
        end
    end
    return String(take!(result))
end

function visualizer_json(x)
    if x isa AbstractString
        return "\"$(json_escape(x))\""
    elseif x isa Integer
        return string(x)
    elseif x isa Bool
        return x ? "true" : "false"
    elseif x isa AbstractVector
        return "[" * join(visualizer_json.(x), ",") * "]"
    elseif x isa AbstractDict
        pairs = [
            visualizer_json(string(k)) * ":" * visualizer_json(v)
            for (k, v) in sort(collect(x); by=first)
        ]
        return "{" * join(pairs, ",") * "}"
    elseif x === nothing
        return "null"
    else
        return visualizer_json(string(x))
    end
end

visualizer_json(M::MackeyFunctor) = visualizer_json(visualizer_data(M))
