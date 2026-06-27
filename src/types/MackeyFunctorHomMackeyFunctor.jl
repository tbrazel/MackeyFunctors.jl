"""
    MackeyFunctorHomMackeyFunctor(M::MackeyFunctor, N::MackeyFunctor)
    InternalHom(M::MackeyFunctor, N::MackeyFunctor)

Represent the internal hom in the category of Mackey functors.

Mathematically, when M and N are G-Mackey functors, the internal hom is the
G-Mackey functor [M,N] whose values are
    [M,N](G/H) = Hom(M_{G/H}, N)
where M_{G/H} represents the G/H-shift of M. When H is a subgroup of K, the
map G/H -> G/K induces maps M_{G/H} <-> M_{G/K} and these yield the
restriction and transfer maps of [M,N]. Likewise for conjugation!
"""