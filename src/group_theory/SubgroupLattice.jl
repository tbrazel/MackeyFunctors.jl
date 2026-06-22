# Given H<G, return all K<G for which H<K
function all_supergroups_from_lattice(G::GapObj, H::GapObj)
    GAP.Globals.IsSubgroup(G, H) || throw(ArgumentError("Second argument must be a subgroup of the first"))
    L = GAP.Globals.LatticeSubgroups(G)

    classes = GAP.Globals.ConjugacyClassesSubgroups(L)

    supergroups = GapObj[]

    for cls in classes
        for K in GAP.Globals.Elements(cls)
            if Bool(GAP.Globals.IsSubgroup(K, H)) && K != H
                push!(supergroups, K)
            end
        end
    end

    return supergroups
end