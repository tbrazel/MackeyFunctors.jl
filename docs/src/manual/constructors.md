# Constructors

One of the major methods used in our constructors is the [`shift`](@ref) operation. Given a Mackey functor $M$ and a finite $G$-set $X$, the general shift operation yields a new Mackey functor $M_X$, defined by the formula $M_X(Y):= M(X\times Y)$. This allows us to take existing constructors and turn them into more general ones. For instance:

- the [`free_mackey_functor`](@ref) at level $H$ is the shift by $H$ of the [`burnside_mackey_functor`](@ref)
- the free cohomological Mackey functor (TODO: export/docstring needed) at level $H$ is the shift by $H$ of the [`constant_mackey_functor`](@ref) at $\mathbb{Z}$.

The shift operation is as follows:

```@docs
shift
```

Other constructors include

```@docs
constant_mackey_functor
burnside_mackey_functor
free_mackey_functor
```
