using MackeyFunctors.AbstractAlgebraLocal: zero_homomorphism

@testset "Kernels and cokernels of Mackey functor homomorphisms" begin
    context = MackeyContext(GAP.Globals.CyclicGroup(2))
    constant_Z = constant_mackey_functor(context, ZZ)
    double_components = Generic.ModuleHomomorphism[
        ModuleHomomorphism(
            value(constant_Z, H_index),
            value(constant_Z, H_index),
            matrix(ZZ, 1, 1, [ZZ(2)]),
        )
        for H_index in eachindex(context.subgroups)
    ]
    double_map = MackeyFunctorHomomorphism(
        constant_Z,
        constant_Z,
        double_components,
    )

    K, inclusion = kernel(double_map)
    Q, projection = cokernel(double_map)

    @test K isa MackeyFunctor
    @test Q isa MackeyFunctor
    @test inclusion isa MackeyFunctorHomomorphism
    @test projection isa MackeyFunctorHomomorphism
    @test inclusion.domain === K
    @test inclusion.codomain === constant_Z
    @test projection.domain === constant_Z
    @test projection.codomain === Q

    @test all(eachindex(context.subgroups)) do H_index
        ngens(value(K, H_index)) == 0
    end
    @test all(eachindex(context.subgroups)) do H_index
        AbstractAlgebra.invariant_factors(value(Q, H_index)) == BigInt[2]
    end
    @test all(eachindex(context.subgroups)) do H_index
        all(gens(value(constant_Z, H_index))) do x
            iszero(projection.components[H_index](double_map.components[H_index](x)))
        end
    end

    identity_map = id_homomorphism(constant_Z)
    identity_kernel, identity_kernel_inclusion = kernel(identity_map)
    identity_cokernel, identity_cokernel_projection = cokernel(identity_map)

    @test all(eachindex(context.subgroups)) do H_index
        ngens(value(identity_kernel, H_index)) == 0 &&
            ngens(value(identity_cokernel, H_index)) == 0
    end
    @test identity_kernel_inclusion.domain === identity_kernel
    @test identity_cokernel_projection.codomain === identity_cokernel

    zero_components = Generic.ModuleHomomorphism[
        zero_homomorphism(value(constant_Z, H_index), value(constant_Z, H_index))
        for H_index in eachindex(context.subgroups)
    ]
    zero_map = MackeyFunctorHomomorphism(
        constant_Z,
        constant_Z,
        zero_components,
    )
    zero_kernel, zero_kernel_inclusion = kernel(zero_map)

    @test all(eachindex(context.subgroups)) do H_index
        all(gens(value(zero_kernel, H_index))) do x
            iszero(zero_map.components[H_index](zero_kernel_inclusion.components[H_index](x)))
        end
    end
end

@testset "Epimorphisms from free Mackey functors" begin
    context = MackeyContext(GAP.Globals.CyclicGroup(2))
    constant_Z = constant_mackey_functor(context, ZZ)
    epi = epimorphism_from_free(constant_Z)
    epi_cokernel, = cokernel(epi)

    @test epi isa MackeyFunctorHomomorphism
    @test epi.codomain === constant_Z
    @test all(eachindex(context.subgroups)) do H_index
        all(iszero, gens(value(epi_cokernel, H_index)))
    end

    reversed_epi = epimorphism_from_free(
        constant_Z;
        level_order=reverse(collect(eachindex(context.subgroups))),
    )
    reversed_cokernel, = cokernel(reversed_epi)
    @test all(eachindex(context.subgroups)) do H_index
        all(iszero, gens(value(reversed_cokernel, H_index)))
    end

    F1 = free_module(ZZ, 1)
    twoF1, = sub(F1, [F1([ZZ(2)])])
    Z2, = quo(F1, twoF1)
    constant_Z2 = constant_mackey_functor(context, Z2)
    quotient_epi = epimorphism_from_free(constant_Z2)
    quotient_cokernel, = cokernel(quotient_epi)

    @test quotient_epi.codomain === constant_Z2
    @test all(eachindex(context.subgroups)) do H_index
        all(iszero, gens(value(quotient_cokernel, H_index)))
    end

    zero_functor = zero_mackey_functor(context, ZZ)
    zero_epi = epimorphism_from_free(zero_functor)
    @test zero_epi.domain === zero_functor
    @test zero_epi.codomain === zero_functor
    @test all(eachindex(context.subgroups)) do H_index
        MackeyFunctors.is_identity_module_homomorphism(zero_epi.components[H_index])
    end
end
