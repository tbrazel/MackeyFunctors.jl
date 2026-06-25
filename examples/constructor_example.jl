using GAP
using AbstractAlgebra
using MackeyFunctors







"""

Construct the free Mackey functor on the underlying subgroup H = context.subgroups[H_idx].

The value at K is the free Z-module on the double cosets.

"""
function cohomological__mackey_functor(context::MackeyContext, H_idx::Int)
    G          = context.group
    subgroups  = context.subgroups
    generators = context.generators
    H          = subgroups[H_idx]


    double_cosets = map(subgroups) do K
        gap_list = GAP.Globals.DoubleCosetRepsAndSizes(G, K, H)
    # Convert to a Julia Vector of (element, size) tuples
        [( gap_list[i][1], gap_list[i][2] ) for i in 1:length(gap_list)] 
    end

    values = map(double_cosets) do dc
        free_module(ZZ, length(dc))
    end


cover_restrictions = map(context.covers) do (J_idx, K_idx)
    dc_J = double_cosets[J_idx]
    dc_K = double_cosets[K_idx]
    K    = subgroups[K_idx]

    mat = zero_matrix(ZZ, length(dc_J), length(dc_K))
    for (k_col, (x, _)) in enumerate(dc_K)
        KxH = GAP.Globals.DoubleCoset(subgroups[K_idx], x, H)
        # Every J-double coset rep y that lies inside KxH contributes a 1
        for (j_row, (y, _)) in enumerate(dc_J)
            if y in KxH
                mat[j_row, k_col] = 1
            end
        end
    end

    ModuleHomomorphism(values[K_idx], values[J_idx], mat)
end


    # Each J-double coset should map into exactly one K-double coset

cover_transfers = map(context.covers) do (J_idx, K_idx)
    K    = subgroups[K_idx]
    dc_J = double_cosets[J_idx]
    dc_K = double_cosets[K_idx]

    mat = zero_matrix(ZZ, length(dc_K), length(dc_J))
    for (j_col, (x, _)) in enumerate(dc_J)
        k_row = findfirst(dc_K) do (y, _)
            x in GAP.Globals.DoubleCoset(K, y, H)
        end
        mat[k_row, j_col] = 1
    end

    ModuleHomomorphism(values[J_idx], values[K_idx], mat)
end
    # g acts on coset reps by left multiplication: KxH |-> gKg^{-1} · gx · H
    # generator_left_conjugation_matrix[n, K_idx] gives the index of gKg^{-1}

    generator_conjugations = Matrix{Generic.ModuleIsomorphism}(
        undef, length(generators), length(subgroups)
    )

    for n in eachindex(generators), K_idx in eachindex(subgroups)
        g         = generators[n]
        gKginv_idx = context.generator_left_conjugation_matrix[n, K_idx]
        gKginv    = subgroups[gKginv_idx]

        dc_K   = double_cosets[K_idx]
        dc_gKg = double_cosets[gKginv_idx]

        # g · (KxH) = gKg^{-1} · gx · H, so we map rep x to rep gx
        mat = zero_matrix(ZZ, length(dc_gKg), length(dc_K))
        for (k_col, (x, _)) in enumerate(dc_K)
            gx    = g * x
            i_row = findfirst(dc_gKg) do (y, _)
                gx in GAP.Globals.DoubleCoset(gKginv, y, H)
            end
            mat[i_row, k_col] = 1
        end

        generator_conjugations[n, K_idx] = ModuleIsomorphism(
            values[K_idx], values[gKginv_idx], mat
        )
    end

    return MackeyFunctor(context, values, cover_restrictions, cover_transfers, generator_conjugations)
end


""" 
Constructs the zero Mackey functor for the group, 
which assigns the zero module to every finite G-set. 
"""


function zero_module_isomorphism(m::AbstractAlgebra.FPModule)
    # Explicitly construct a ModuleIsomorphism for a zero (rank-0) module
    return Generic.ModuleIsomorphism(m, m, ZZ[;])
end

function zero_mackey_functor(context::MackeyContext)

        zero_underlying = free_module(ZZ, 0)
        zero_fixed = free_module(ZZ, 0)

        res = hom(zero_underlying, zero_fixed, matrix(0))
        tr = hom(zero_fixed, zero_underlying, matrix(0))
end 
