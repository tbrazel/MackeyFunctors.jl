using MackeyFunctors.AbstractAlgebraLocal: HomModule, underlying_module, as_hom_module_element, as_homomorphism, TensorProduct, tensor_product, tensor_product_element, submodules_matrix

@testset "module morphisms" begin
    M = FreeModule(ZZ, 1)
    m1, m2, m3 = M([ZZ(1)]), M([ZZ(2)]), M([ZZ(3)])
    N2, _ = sub(M, [m2])
    N3, _ = sub(M, [m3])
    L2, g2 = quo(M, N2)
    L3, g3 = quo(M, N3)
    id = identity_matrix(ZZ, 1)

    @test_throws ArgumentError ModuleHomomorphism(L2, M, id)
    @test_throws ArgumentError ModuleHomomorphism(L2, L3, id)
    @test_throws ArgumentError ModuleIsomorphism(L2, L3, id)
    @test ModuleHomomorphism(L2, L3, ZZ(3)*id) == AbstractAlgebra.ModuleHomomorphism(L2, L3, ZZ(3)*id)
    @test ModuleHomomorphism(M, L3, id) == AbstractAlgebra.ModuleHomomorphism(M, L3, id)
    @test ModuleHomomorphism(L3, L3, ZZ(2)*id) == AbstractAlgebra.ModuleHomomorphism(L3, L3, ZZ(2)*id)

    @test_throws ArgumentError ModuleHomomorphism(L2, M, [m1])
    @test_throws ArgumentError ModuleHomomorphism(L2, L3, [g3(m1)])
    @test_throws ArgumentError ModuleIsomorphism(L2, L3, [g3(m1)])
    @test ModuleHomomorphism(L2, L3, [g3(m3)]) == AbstractAlgebra.ModuleHomomorphism(L2, L3, [g3(m3)])
    @test ModuleHomomorphism(M, L3, [g3(m1)]) == AbstractAlgebra.ModuleHomomorphism(M, L3, [g3(m1)])
    @test ModuleHomomorphism(L3, L3, [g3(m2)]) == AbstractAlgebra.ModuleHomomorphism(L3, L3, [g3(m2)])
end

@testset "Hom modules" begin
    F1 = free_module(ZZ, 1)
    twoF1, = sub(F1, [F1([ZZ(2)])])
    fourF1, = sub(F1, [F1([ZZ(4)])])
    Z2, = quo(F1, twoF1)
    Z4, = quo(F1, fourF1)

    hom_Z2_Z2 = HomModule(Z2, Z2)
    @test underlying_module(hom_Z2_Z2) isa AbstractAlgebra.FPModule
    @test ngens(hom_Z2_Z2) == 1
    @test length(relations(hom_Z2_Z2)) == 1
    @test relations(hom_Z2_Z2)[1][1, 1] == ZZ(2)

    id_Z2 = ModuleHomomorphism(Z2, Z2, matrix(ZZ, 1, 1, [ZZ(1)]))
    id_element = as_hom_module_element(hom_Z2_Z2, id_Z2)
    @test 2*id_element == zero(underlying_module(hom_Z2_Z2))
    @test as_homomorphism(hom_Z2_Z2, id_element) == id_Z2

    F2 = free_module(ZZ, 2)
    hom_Z2_power = HomModule(F2, Z2)
    @test ngens(hom_Z2_power) == 2
    @test length(relations(hom_Z2_power)) == 2
    @test Set([[relation[1, i] for i in 1:ncols(relation)] for relation in relations(hom_Z2_power)]) ==
        Set([[ZZ(2), ZZ(0)], [ZZ(0), ZZ(2)]])
    first_projection = ModuleHomomorphism(F2, Z2, matrix(ZZ, 2, 1, [ZZ(1), ZZ(0)]))
    @test as_homomorphism(hom_Z2_power, as_hom_module_element(hom_Z2_power, first_projection)) == first_projection

    hom_Z2_Z4 = HomModule(Z2, Z4)
    @test ngens(hom_Z2_Z4) == 1
    killed_by_two = ModuleHomomorphism(Z2, Z4, matrix(ZZ, 1, 1, [ZZ(2)]))
    killed_by_two_element = as_hom_module_element(hom_Z2_Z4, killed_by_two)
    @test as_homomorphism(hom_Z2_Z4, killed_by_two_element) == killed_by_two

    M = free_module(QQ, 2)
    N = free_module(QQ, 3)
    free_hom_module = HomModule(M, N)
    f = ModuleHomomorphism(M, N, matrix(QQ, [1 2 3; 4 5 6]))
    f_element = as_hom_module_element(free_hom_module, f)
    @test matrix(as_homomorphism(free_hom_module, f_element)) == matrix(f)
end

@testset "Tensor products of FPModules" begin
    F1 = free_module(ZZ, 1)
    F2 = free_module(ZZ, 2)
    twoF1, = sub(F1, [F1([ZZ(2)])])
    threeF1, = sub(F1, [F1([ZZ(3)])])
    fourF1, = sub(F1, [F1([ZZ(4)])])
    Z2, = quo(F1, twoF1)
    Z3, = quo(F1, threeF1)
    Z4, = quo(F1, fourF1)

    # Z^2 * Z
    free_tensor = TensorProduct(F2, F1)
    @test underlying_module(free_tensor) isa AbstractAlgebra.FPModule
    @test ngens(free_tensor) == 2
    @test isempty(relations(free_tensor))
    @test tensor_product_element(free_tensor, F2([ZZ(3), ZZ(5)]), gen(F1, 1)) ==
        underlying_module(free_tensor)([ZZ(3), ZZ(5)])

    #Z/2 * Z/4
    z2_tensor_z4 = tensor_product(Z2, Z4)
    @test underlying_module(z2_tensor_z4) isa AbstractAlgebra.FPModule
    @test ngens(z2_tensor_z4) == 1
    @test AbstractAlgebra.invariant_factors(underlying_module(z2_tensor_z4)) == BigInt[2]
    pure_tensor = tensor_product_element(z2_tensor_z4, gen(Z2, 1), gen(Z4, 1))
    @test 2*pure_tensor == zero(underlying_module(z2_tensor_z4))
    @test pure_tensor == gen(underlying_module(z2_tensor_z4), 1)

    # Z/2 * Z/3
    z2_tensor_z3 = TensorProduct(Z2, Z3)
    @test ngens(z2_tensor_z3) == 0
    @test tensor_product_element(z2_tensor_z3, gen(Z2, 1), gen(Z3, 1)) ==
        zero(underlying_module(z2_tensor_z3))

    # Z^2 * Z/4
    bilinear_tensor = TensorProduct(F2, Z4)
    @test tensor_product_element(bilinear_tensor, F2([ZZ(2), ZZ(3)]), gen(Z4, 1)) ==
        2*tensor_product_element(bilinear_tensor, gen(F2, 1), gen(Z4, 1)) +
        3*tensor_product_element(bilinear_tensor, gen(F2, 2), gen(Z4, 1))

    # Tensor two maps together
    left_map = ModuleHomomorphism(F2, F2, matrix(ZZ, [1 2; 3 4]))
    right_map = ModuleHomomorphism(F1, F1, matrix(ZZ, 1, 1, [ZZ(5)]))
    tensor_map = tensor_product(left_map, right_map)
    source_tensor = TensorProduct(F2, F1)
    target_tensor = TensorProduct(F2, F1)
    @test domain(tensor_map) == underlying_module(source_tensor)
    @test codomain(tensor_map) == underlying_module(target_tensor)
    @test tensor_map(gen(domain(tensor_map), 1)) == underlying_module(target_tensor)([ZZ(5), ZZ(10)])
    @test tensor_map(gen(domain(tensor_map), 2)) == underlying_module(target_tensor)([ZZ(15), ZZ(20)])

    # id_{Z/2} * (Z/2 -> Z/4) gives zero
    id_Z2 = ModuleHomomorphism(Z2, Z2, matrix(ZZ, 1, 1, [ZZ(1)]))
    double_into_Z4 = ModuleHomomorphism(Z2, Z4, matrix(ZZ, 1, 1, [ZZ(2)]))
    torsion_tensor_map = tensor_product(id_Z2, double_into_Z4)
    @test ngens(domain(torsion_tensor_map)) == 1
    @test ngens(codomain(torsion_tensor_map)) == 1
    @test iszero(torsion_tensor_map)

    # Can't tensor modules over different base rings
    F_QQ_1 = free_module(QQ, 1)
    @test_throws ArgumentError TensorProduct(F1, F_QQ_1)
end
