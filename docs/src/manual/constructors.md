# Constructors
## `G`-Module Constructors

```@docs
GModule
permutation_module
```

## Mackey Functor Constructors

### Shifts
One of the major methods used in our constructors is the [`shift`](@ref) operation. Given a Mackey functor $M$ and a finite $G$-set $X$, the general shift operation yields a new Mackey functor $M_X$, defined by the formula $M_X(Y):= M(X\times Y)$. This allows us to take existing constructors and turn them into more general ones. For instance:

- the [`free_mackey_functor`](@ref) at level $H$ is the shift by $H$ of the [`burnside_mackey_functor`](@ref)
- the free cohomological Mackey functor (TODO: export/docstring needed) at level $H$ is the shift by $H$ of the [`constant_mackey_functor`](@ref) at $\mathbb{Z}$.

The shift operation is as follows:

```@docs
shift
```
We can also shift a [Mackey functor homomorphism](@ref MackeyFunctorHomomorphism) by a subgroup index.

### Fixed point Mackey functors

Let $M$ be a (left) $G$-module (of type [`GModule`](@ref)). Then we obtain the *fixed point Mackey functor* $\text{FP}(M)$, whose value at level $H$ is the fixed submodule $M^H$. For $H\le K$, restriction is given by the natural inclusion of submodules $M^K \subseteq M^H$, and transfer is obtained by fixing left coset representatives $K = \cup_i k_i H$, and defining the map $$M^H \to M^K; x \mapsto \sum_i k_i x$$ (this can be seen to be independent of the choice of coset representatives). Conjugation by $g\in G$ is given by the natural map
$$M^H \to M^{gHg^{-1}}; x \mapsto gx.$$
```@docs
fixedpoint_mackey_functor
```

### Constant Mackey functors
For any module $M$, we define the *constant Mackey functor valued at $M$*, denoted $\underline{M}$, to be the fixed point Mackey functor attached to $M$, where $M$ is considered as a $G$-module with trivial action.
```@docs
constant_mackey_functor
```
As an example we can compute the *zero Mackey functor*:
```@docs
zero_mackey_functor
```

### The Burnside Mackey functor
For a group $G$, the *Burnside ring* $A(G)$ is the ring of isomorphism classes of virtual finite $G$-sets. For any $G$, this promotes to a $G$-Mackey functor, defined by the data $H\mapsto A(H)$, with natural restrictions. Transfers are obtained by induction of finite $G$-sets.
```@docs
burnside_mackey_functor
```
### Free Mackey functors
For each $H\le G$, we define the *free Mackey functor at level $H$* to be the shift of the [Burnside Mackey functor](@ref burnside_mackey_functor) by $G/H$.
```@docs
free_mackey_functor
```
