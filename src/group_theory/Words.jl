const Group = GapObj
const GroupElement = GapObj
const GeneratorIndex::DataType = Int
const GeneratorWord = Vector{Tuple{GeneratorIndex,Int}}

function _generator_word_from_extrep(extrep)::GeneratorWord
    word = Vector{Int}(extrep)

    iseven(length(word)) ||
        throw(ArgumentError("GAP returned an invalid free-group word representation"))

    return [
        (
            word[i],
            word[i+1],
        )
        for i in 1:2:length(word)
    ]
end


function generator_word(group::Group, element::GroupElement)::GeneratorWord
    return _generator_word_from_extrep(
        GAP.Globals.ExtRepOfObj(GAP.Globals.Factorization(group, element)),
    )
end

"""
    generator_relations(group)
    generator_relations(group, generators)

Return GAP's relators for a finite presentation of `group` on `generators`.

Each relator is returned as a `GeneratorWord`, with the same `(generator_index,
exponent)` convention.
"""
function generator_relations(group::Group, generators::Vector{GroupElement})::Vector{GeneratorWord}
    fp_isomorphism =
        GAP.Globals.IsomorphismFpGroupByGenerators(group, GapObj(generators; recursive=true))
    fp_group = GAP.Globals.Image(fp_isomorphism)

    return [
        _generator_word_from_extrep(GAP.Globals.ExtRepOfObj(relator))
        for relator in GAP.Globals.RelatorsOfFpGroup(fp_group)
    ]
end

function generator_relations(group::Group)::Vector{GeneratorWord}
    return generator_relations(
        group,
        Vector{GroupElement}(GAP.Globals.GeneratorsOfGroup(group)),
    )
end
