# tensor products of FP modules

struct TensorProduct{R <: RingElement, M <: AbstractAlgebra.FPModule{R}, F} <: AbstractAlgebra.FPModule{R}
    mod::M  # underlying module (so far, free module or quotient module)
    f::F    # multilinear structure map with codomain `mod`
end

underlying_module(M::TensorProduct) = M.mod

AbstractAlgebra.base_ring_type(::Type{<:TensorProduct{R}}) where R <: RingElement = parent_type(R)
AbstractAlgebra.coefficient_ring_type(T::Type{<:TensorProduct}) = base_ring_type(T)

for f in (:base_ring, :number_of_generators, :relations)
    @eval AbstractAlgebra.$f(M::TensorProduct) = $f(M.mod)
end
AbstractAlgebra.coefficient_ring(M::TensorProduct) = base_ring(M)

AbstractAlgebra.gen(M::TensorProduct, i::Integer) = M(gen(M.mod, i))
AbstractAlgebra.gens(M::TensorProduct) = map(M, gens(M.mod))

"""
    structure_map(t::TensorProduct{R}) where R <: RingElement

Return the canonical multilinear structure map of the tensor product `t`.
"""
function structure_map(t::TensorProduct{R}) where R <: RingElement
    (vs::AbstractAlgebra.FPModuleElem{R}...) -> t(t.f(vs...))
end

struct TensorProductElem{R <: RingElement, T <: TensorProduct{R}} <: AbstractAlgebra.FPModuleElem{R}
    parent::T
    v::Generic.MatSpaceElem{R}  # Generic._matrix uses this field to get coordinates
end

AbstractAlgebra.parent_type(::Type{TensorProductElem{R, T}}) where {R <: RingElement, T <: TensorProduct{R}} = T
AbstractAlgebra.elem_type(::Type{T}) where {R <: RingElement, T <: TensorProduct{R}} = TensorProductElem{R, T}
AbstractAlgebra.parent(v::TensorProductElem) = v.parent

(t::TensorProduct)(v::AbstractVector) = t(t.mod(v))
(t::TensorProduct{R})(v::Generic.MatSpaceElem{R}) where R = t(t.mod(v))

function (t::TensorProduct{R})(v::AbstractAlgebra.FPModuleElem{R}) where R
    parent(v) === t.mod || throw(ArgumentError("incompatible arguments"))
    TensorProductElem(t, Generic._matrix(v))
end

"""
    tensor_product(MS::FPModule{R}...) where R <: RingElement

Return a pair `(M, f)` where `M` is the tensor product of the modules `MS`
and `f` the canonical multilinear structure map from `MS` to `M`.
All argument modules must be defined over the same base ring.
If all arguments are of type `Generic.FreeModule`, then `underlying_module(M)` is too.

See also [`structure_map`](@ref).
"""
tensor_product() = throw(ArgumentError("empty tensor products are not supported"))

function tensor_product(ms::AbstractAlgebra.FPModule{R}...) where R <: RingElement
    isempty(ms) && throw(ArgumentError("empty tensor products are not supported"))
    allequal(base_ring, ms) || throw(ArgumentError("all modules must have the same base ring"))
    ranks = map(ngens, ms)
    M = FreeModule(base_ring(ms[1]), prod(ranks))

    # structure map for tensor product of free modules
    f = function(vs::AbstractAlgebra.FPModuleElem{R}...)
        length(vs) == length(ms) || throw(ArgumentError("wrong number of arguments"))
        all(((v, m),) -> parent(v) === m, zip(vs, ms)) || throw(ArgumentError("arguments must be elements of the tensor factors"))
        vm = map(prod, Iterators.product(map(Generic._matrix, vs)...))
        M(reshape(vm, :))
    end

    if all(m -> m isa Generic.FreeModule, ms)
        # free tensor product
        TM =  TensorProduct(M, f)
        return TM, structure_map(TM)
    end

    # non-free tensor product
    civ = reshape(CartesianIndices(ranks), :)
    tensor_relations = [M([@inbounds relation[ci[k]] for ci in civ]) for (k, m) in enumerate(ms) for relation in relations(m)]
    N, _ = sub(M, tensor_relations)
    Q, q = quo(M, N)
    TQ = TensorProduct(Q, q ∘ f)
    return TQ, structure_map(TQ)
end

function tensor_product(ms::AbstractAlgebra.FPModule...)
    isempty(ms) && throw(ArgumentError("empty tensor products are not supported"))
    allequal(base_ring, ms) || throw(ArgumentError("all modules must have the same base ring"))
    throw(ArgumentError("all modules must have compatible coefficient element types"))
end

"""
    tensor_product(fs::Union{Generic.ModuleHomomorphism{R}, Generic.ModuleIsomorphism{R}}...;
            [domain], [codomain]) where R <: RingElement

Return the linear morphism from the tensor product of the domains of the maps `fs`
to the tensor product of the codomains. If `domain` or `codomain` are not given as
keyword arguments, they are constructed anew.

If all morphisms are of type `ModuleIsomorphism`, then so is their tensor product.
Otherwise it is a `ModuleHomomorphism`.
"""
function tensor_product(fs::Union{Generic.ModuleHomomorphism{R}, Generic.ModuleIsomorphism{R}}...;
        domain = first(tensor_product(map(domain, fs)...)),
        codomain = first(tensor_product(map(codomain, fs)...)),
    ) where R <: RingElement
    factor_gens = map(gens ∘ AbstractAlgebraLocal.domain, fs)
    a = [Generic._matrix(codomain.f(map((f, v) -> f(v), fs, vs)...)) for vs in Iterators.product(factor_gens...)]
    b = matrix([v[i] for v in reshape(a, :), i in 1:ngens(codomain)])
    if all(f -> f isa Generic.ModuleIsomorphism, fs)
        ModuleIsomorphism(domain, codomain, b)  # if the constructor allowed, we could also provide the inverse map
    else
        ModuleHomomorphism(domain, codomain, b)
    end
end
