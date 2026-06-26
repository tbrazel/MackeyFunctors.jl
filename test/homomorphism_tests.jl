using MackeyFunctors.AbstractAlgebraLocal: HomModule,
    as_hom_module_element,
    as_homomorphism,
    underlying_module

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

@testset "Mackey functor Hom modules" begin
    context = MackeyContext(GAP.Globals.CyclicGroup(2))
    M = free_module(ZZ, 1)
    constant_Z = constant_mackey_functor(context, M)
    endomorphism_module = HomModule(constant_Z, constant_Z)

    @test endomorphism_module isa MackeyFunctorHomModule
    @test underlying_module(endomorphism_module) isa AbstractAlgebra.FPModule
    @test ngens(endomorphism_module) == 1
    @test isempty(relations(endomorphism_module))

    identity_mackey_hom = MackeyFunctors.id_homomorphism(constant_Z)
    identity_element = as_hom_module_element(endomorphism_module, identity_mackey_hom)
    round_trip_identity = as_homomorphism(endomorphism_module, identity_element)
    @test round_trip_identity isa MackeyFunctorHomomorphism
    @test all(eachindex(context.subgroups)) do subgroup_index
        MackeyFunctors.map_eq(
            round_trip_identity.components[subgroup_index],
            identity_mackey_hom.components[subgroup_index],
        )
    end

    F1 = free_module(ZZ, 1)
    twoF1, = sub(F1, [F1([ZZ(2)])])
    Z2, = quo(F1, twoF1)
    constant_Z2 = constant_mackey_functor(context, Z2)
    hom_Z2_Z = HomModule(constant_Z2, constant_Z)
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
    s3_endomorphism_module = HomModule(s3_functor, s3_functor)
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
