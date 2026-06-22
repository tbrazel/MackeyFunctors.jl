# Return a Julia Vector of all subgroups of a GAP group
function subgroups(G::GapObj)

    GAP.Globals.IsGroup(G) ||

        throw(ArgumentError("subgroups(G) can only be called when G is a group"))

    return Vector{GapObj}(GAP.Globals.AllSubgroups(G))

end