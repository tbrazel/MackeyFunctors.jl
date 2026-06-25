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

function generator_word_from_isomorphism(
    fp_isomorphism::GapObj,
    element::GroupElement,
)::GeneratorWord
    return _generator_word_from_extrep(
        GAP.Globals.ExtRepOfObj(
            GAP.Globals.ImagesRepresentative(fp_isomorphism, element),
        ),
    )
end

function generator_relations_from_isomorphism(fp_isomorphism::GapObj)::Vector{GeneratorWord}
    fp_group = GAP.Globals.Image(fp_isomorphism)
    return [
        _generator_word_from_extrep(GAP.Globals.ExtRepOfObj(relator))
        for relator in GAP.Globals.RelatorsOfFpGroup(fp_group)
    ]
end