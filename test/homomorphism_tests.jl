using MackeyFunctors.AbstractAlgebraLocal: Hom,
    as_hom_module_element,
    as_homomorphism,
    underlying_module,
    zero_homomorphism

function _zero_mackey_functor_homomorphism(
    domain_mf::MackeyFunctor,
    codomain_mf::MackeyFunctor,
)
    return MackeyFunctorHomomorphism(
        domain_mf,
        codomain_mf,
        Generic.ModuleHomomorphism[
            zero_homomorphism(domain_mf.values[i], codomain_mf.values[i])
            for i in eachindex(domain_mf.context.subgroups)
        ],
    )
end

@testset "Module map composition order" begin
    M = free_module(ZZ, 2)
    A_hom = ModuleHomomorphism(M, M, matrix(ZZ, [1 1; 0 1]))
    B_hom = ModuleHomomorphism(M, M, matrix(ZZ, [1 0; 1 1]))
    A_iso = ModuleIsomorphism(M, M, matrix(ZZ, [1 1; 0 1]))
    B_iso = ModuleIsomorphism(M, M, matrix(ZZ, [1 0; 1 1]))
    id_hom = ModuleHomomorphism(M, M, identity_matrix(ZZ, 2))
    id_iso = ModuleIsomorphism(M, M, identity_matrix(ZZ, 2))

    C4 = GAP.Globals.CyclicGroup(4)
    c4_context = MackeyContext(C4)
    trivial = findfirst(H -> Int(GAP.Globals.Size(H)) == 1, c4_context.subgroups)
    whole = findfirst(
        H -> Int(GAP.Globals.Size(H)) == Int(GAP.Globals.Size(C4)),
        c4_context.subgroups,
    )
    path = c4_context.paths[(trivial, whole)]
    @test length(path) == 2

    values = [M for _ in c4_context.subgroups]
    restrictions = [
        i == path[1] ? A_hom : i == path[2] ? B_hom : id_hom
        for i in eachindex(c4_context.covers)
    ]
    transfers = copy(restrictions)
    conjugations = [id_iso for _ in c4_context.generators, _ in c4_context.subgroups]
    c4_mackey_functor = MackeyFunctor(
        c4_context,
        values,
        restrictions,
        transfers,
        conjugations;
        verify = false,
    )

    res = MackeyFunctors.restriction(c4_mackey_functor, trivial, whole)
    tr = MackeyFunctors.transfer(c4_mackey_functor, trivial, whole)

    @test all(gens(M)) do x
        res(x) == restrictions[path[1]](restrictions[path[2]](x))
    end
    @test all(gens(M)) do x
        tr(x) == transfers[path[2]](transfers[path[1]](x))
    end

    S3 = GAP.Globals.SymmetricGroup(3)
    s3_context = MackeyContext(S3)
    s3_restrictions = [id_hom for _ in s3_context.covers]
    s3_transfers = [id_hom for _ in s3_context.covers]
    s3_conjugations = [
        i == 1 ? A_iso : i == 2 ? B_iso : id_iso
        for i in eachindex(s3_context.generators), _ in s3_context.subgroups
    ]
    s3_mackey_functor = MackeyFunctor(
        s3_context,
        [M for _ in s3_context.subgroups],
        s3_restrictions,
        s3_transfers,
        s3_conjugations;
        verify = false,
    )

    conjugation_word = [(1, 1), (2, 1)]
    conj = MackeyFunctors.conjugation(s3_mackey_functor, conjugation_word, 1)
    @test all(gens(M)) do x
        conj(x) == A_iso(B_iso(x))
    end
end

@testset "Mackey functor homomorphisms" begin
    context = MackeyContext(GAP.Globals.CyclicGroup(2))
    M = free_module(ZZ, 1)
    mackey_functor = constant_mackey_functor(context, M)
    id_component = ModuleHomomorphism(M, M, identity_matrix(ZZ, 1))
    components = [id_component for _ in context.subgroups]

    homomorphism = MackeyFunctorHomomorphism(
        mackey_functor,
        mackey_functor,
        components,
    )

    @test homomorphism isa MackeyFunctorHomomorphism
    @test homomorphism.context == context
    @test homomorphism.domain === mackey_functor
    @test homomorphism.codomain === mackey_functor
    @test length(homomorphism.components) == length(context.subgroups)

    identity_homomorphism = MackeyFunctors.id_homomorphism(mackey_functor)
    @test identity_homomorphism isa MackeyFunctorHomomorphism
    @test identity_homomorphism.domain === mackey_functor
    @test identity_homomorphism.codomain === mackey_functor
    @test all(identity_homomorphism.components) do component
        MackeyFunctors.is_identity_module_homomorphism(component)
    end

    @test_throws ArgumentError MackeyFunctorHomomorphism(
        mackey_functor,
        mackey_functor,
        components[1:end-1],
    )

    wrong_domain = free_module(ZZ, 1)
    wrong_domain_component =
        ModuleHomomorphism(wrong_domain, M, zero_matrix(ZZ, 1, 1))
    wrong_domain_components = Generic.ModuleHomomorphism[components...]
    wrong_domain_components[1] = wrong_domain_component
    @test_throws ArgumentError MackeyFunctorHomomorphism(
        mackey_functor,
        mackey_functor,
        wrong_domain_components,
    )

    wrong_codomain = free_module(ZZ, 1)
    wrong_codomain_component =
        ModuleHomomorphism(M, wrong_codomain, zero_matrix(ZZ, 1, 1))
    wrong_codomain_components = Generic.ModuleHomomorphism[components...]
    wrong_codomain_components[1] = wrong_codomain_component
    @test_throws ArgumentError MackeyFunctorHomomorphism(
        mackey_functor,
        mackey_functor,
        wrong_codomain_components,
    )

    N = free_module(ZZ, 2)
    id_hom = ModuleHomomorphism(N, N, identity_matrix(ZZ, 2))
    id_iso = ModuleIsomorphism(N, N, identity_matrix(ZZ, 2))
    A_iso = ModuleIsomorphism(N, N, matrix(ZZ, [1 1; 0 1]))
    B_iso = ModuleIsomorphism(N, N, matrix(ZZ, [1 0; 1 1]))

    s3_context = MackeyContext(GAP.Globals.SymmetricGroup(3))
    s3_functor = MackeyFunctor(
        s3_context,
        [N for _ in s3_context.subgroups],
        [id_hom for _ in s3_context.covers],
        [id_hom for _ in s3_context.covers],
        [
            subgroup_index == 1 ? A_iso : subgroup_index == 2 ? B_iso : id_iso
            for _ in eachindex(s3_context.generators), subgroup_index in eachindex(s3_context.subgroups)
        ];
        verify = false,
    )
    @test MackeyFunctorHomomorphism(
        s3_functor,
        s3_functor,
        [id_hom for _ in s3_context.subgroups],
    ) isa MackeyFunctorHomomorphism
end

@testset "Matrices of Mackey functor homomorphisms" begin
    context = MackeyContext(GAP.Globals.CyclicGroup(2))
    constant_Z = constant_mackey_functor(context, ZZ)

    function scalar_hom(a)
        return MackeyFunctorHomomorphism(
            constant_Z,
            constant_Z,
            Generic.ModuleHomomorphism[
                ModuleHomomorphism(
                    value(constant_Z, H_index),
                    value(constant_Z, H_index),
                    matrix(ZZ, 1, 1, [ZZ(a)]),
                )
                for H_index in eachindex(context.subgroups)
            ],
        )
    end

    id_map = scalar_hom(1)
    double_map = scalar_hom(2)
    zero_map = scalar_hom(0)

    row_map = MackeyFunctorHomomorphism([id_map double_map])
    @test row_map isa MackeyFunctorHomomorphism
    @test all(eachindex(context.subgroups)) do H_index
        source_gen = gens(value(row_map.domain, H_index))[1]
        target_gens = gens(value(row_map.codomain, H_index))
        row_map.components[H_index](source_gen) ==
            target_gens[1] + ZZ(2) * target_gens[2]
    end

    column_map = MackeyFunctorHomomorphism([id_map, double_map])
    @test column_map isa MackeyFunctorHomomorphism
    @test all(eachindex(context.subgroups)) do H_index
        source_gens = gens(value(column_map.domain, H_index))
        target_gen = gens(value(column_map.codomain, H_index))[1]
        column_map.components[H_index](source_gens[1]) == target_gen &&
            column_map.components[H_index](source_gens[2]) == ZZ(2) * target_gen
    end

    row_from_vector = block_homomorphism(
        [id_map, double_map];
        orientation=:row,
    )
    @test all(eachindex(context.subgroups)) do H_index
        source_gen = gens(value(row_from_vector.domain, H_index))[1]
        target_gens = gens(value(row_from_vector.codomain, H_index))
        row_from_vector.components[H_index](source_gen) ==
            target_gens[1] + ZZ(2) * target_gens[2]
    end

    matrix_map = block_homomorphism([
        id_map zero_map
        double_map id_map
    ])
    @test matrix_map isa MackeyFunctorHomomorphism
    @test all(eachindex(context.subgroups)) do H_index
        source_gens = gens(value(matrix_map.domain, H_index))
        target_gens = gens(value(matrix_map.codomain, H_index))
        matrix_map.components[H_index](source_gens[1]) == target_gens[1] &&
            matrix_map.components[H_index](source_gens[2]) ==
                ZZ(2) * target_gens[1] + target_gens[2]
    end

    diagonal_map = MackeyFunctors.direct_sum(id_map, double_map)
    @test diagonal_map isa MackeyFunctorHomomorphism
    @test all(eachindex(context.subgroups)) do H_index
        source_gens = gens(value(diagonal_map.domain, H_index))
        target_gens = gens(value(diagonal_map.codomain, H_index))
        diagonal_map.components[H_index](source_gens[1]) == target_gens[1] &&
            diagonal_map.components[H_index](source_gens[2]) ==
                ZZ(2) * target_gens[2]
    end

    @test_throws ArgumentError MackeyFunctorHomomorphism(
        [id_map, double_map];
        orientation=:diagonal,
    )
end

@testset "Mackey functor Hom modules" begin
    context = MackeyContext(GAP.Globals.CyclicGroup(2))
    M = free_module(ZZ, 1)
    constant_Z = constant_mackey_functor(context, M)
    endomorphism_module = Hom(constant_Z, constant_Z)

    @test endomorphism_module isa MackeyFunctorHomModule
    @test underlying_module(endomorphism_module) isa AbstractAlgebra.FPModule
    @test ngens(endomorphism_module) == 1
    @test isempty(relations(endomorphism_module))

    identity_mackey_hom = MackeyFunctors.id_homomorphism(constant_Z)
    identity_element = as_hom_module_element(endomorphism_module, identity_mackey_hom)
    round_trip_identity = as_homomorphism(endomorphism_module, identity_element)
    @test round_trip_identity isa MackeyFunctorHomomorphism
    @test all(eachindex(context.subgroups)) do subgroup_index
        round_trip_identity.components[subgroup_index] ==
            identity_mackey_hom.components[subgroup_index]
    end

    precompose_identity = precomposition_map(
        endomorphism_module,
        endomorphism_module,
        identity_mackey_hom,
    )
    postcompose_identity = postcomposition_map(
        endomorphism_module,
        endomorphism_module,
        identity_mackey_hom,
    )
    @test all(gens(underlying_module(endomorphism_module))) do x
        precompose_identity(x) == x
    end
    @test all(gens(underlying_module(endomorphism_module))) do x
        postcompose_identity(x) == x
    end

    zero_functor = zero_mackey_functor(context, ZZ)
    zero_to_constant = _zero_mackey_functor_homomorphism(
        zero_functor,
        constant_Z,
    )
    constant_to_zero = _zero_mackey_functor_homomorphism(
        constant_Z,
        zero_functor,
    )
    hom_zero_constant = Hom(zero_functor, constant_Z)
    hom_constant_zero = Hom(constant_Z, zero_functor)

    precompose_zero = precomposition_map(
        endomorphism_module,
        hom_zero_constant,
        zero_to_constant,
    )
    postcompose_zero = postcomposition_map(
        endomorphism_module,
        hom_constant_zero,
        constant_to_zero,
    )
    @test all(gens(underlying_module(endomorphism_module))) do x
        iszero(precompose_zero(x))
    end
    @test all(gens(underlying_module(endomorphism_module))) do x
        iszero(postcompose_zero(x))
    end

    @test_throws ArgumentError precomposition_map(
        endomorphism_module,
        endomorphism_module,
        MackeyFunctors.id_homomorphism(zero_functor),
    )
    @test_throws ArgumentError postcomposition_map(
        endomorphism_module,
        endomorphism_module,
        MackeyFunctors.id_homomorphism(zero_functor),
    )

    F1 = free_module(ZZ, 1)
    twoF1, = sub(F1, [F1([ZZ(2)])])
    Z2, = quo(F1, twoF1)
    constant_Z2 = constant_mackey_functor(context, Z2)
    hom_Z2_Z = Hom(constant_Z2, constant_Z)
    @test ngens(hom_Z2_Z) == 0

    id_component = ModuleHomomorphism(M, M, identity_matrix(ZZ, 1))
    bad_component_element = as_hom_module_element(
        endomorphism_module.component_hom_modules[1],
        id_component,
    )
    bad_ambient_element =
        endomorphism_module.ambient_injections[1](bad_component_element)
    has_bad_lift, = has_preimage_with_preimage(
        endomorphism_module.inclusion,
        bad_ambient_element,
    )
    @test !has_bad_lift

    N = free_module(ZZ, 2)
    id_N_hom = ModuleHomomorphism(N, N, identity_matrix(ZZ, 2))
    id_N_iso = ModuleIsomorphism(N, N, identity_matrix(ZZ, 2))
    A_hom = ModuleHomomorphism(N, N, matrix(ZZ, [1 1; 0 1]))
    A_iso = ModuleIsomorphism(N, N, matrix(ZZ, [1 1; 0 1]))
    B_iso = ModuleIsomorphism(N, N, matrix(ZZ, [1 0; 1 1]))

    s3_context = MackeyContext(GAP.Globals.SymmetricGroup(3))
    s3_functor = MackeyFunctor(
        s3_context,
        [N for _ in s3_context.subgroups],
        [id_N_hom for _ in s3_context.covers],
        [id_N_hom for _ in s3_context.covers],
        [
            subgroup_index == 1 ? A_iso : subgroup_index == 2 ? B_iso : id_N_iso
            for _ in eachindex(s3_context.generators), subgroup_index in eachindex(s3_context.subgroups)
        ];
        verify = false,
    )
    s3_endomorphism_module = Hom(s3_functor, s3_functor)
    s3_identity_element = as_hom_module_element(
        s3_endomorphism_module,
        MackeyFunctors.id_homomorphism(s3_functor),
    )
    @test as_homomorphism(
        s3_endomorphism_module,
        s3_identity_element,
    ) isa MackeyFunctorHomomorphism

    bad_conjugation_ambient = zero(s3_endomorphism_module.ambient_module)
    for subgroup_index in eachindex(s3_context.subgroups)
        bad_component = as_hom_module_element(
            s3_endomorphism_module.component_hom_modules[subgroup_index],
            A_hom,
        )
        bad_conjugation_ambient +=
            s3_endomorphism_module.ambient_injections[subgroup_index](
                bad_component,
            )
    end
    has_bad_conjugation_lift, = has_preimage_with_preimage(
        s3_endomorphism_module.inclusion,
        bad_conjugation_ambient,
    )
    @test !has_bad_conjugation_lift
end
