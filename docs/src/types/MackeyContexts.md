# Mackey contexts

One of our time-saving design choices is the idea of a [`MackeyContext`](@ref), which precomputes relevant information about a finite group before working with Mackey functors for that group. The idea is to make the construction of a Mackey functor a two-step process: step 1 is to build the context, and then step 2 is to build a Mackey functor for the context. This allows us to avoid repeating certain costly computations.

For example, any $C_2$-Mackey functor $M$ must satisfy $\operatorname{res}^{C_2}_e(\operatorname{tr}_e^{C_2}(x)) = x + \operatorname{cong}_{\gamma,e}(x)$ for all elements $x \in M(C_2/e)$ -- this is an example of a *double-coset formula*. For any group $G$, there is a finite list of specific double-coset formulae that a $G$-Mackey functor must satisfy, which we would like to verify when we attempt to construct a $G$-Mackey functor in software. It is easy to specify these formulae mathematically in one line, but it can be slow to enumerate them computationally. So, we precompute all of the double-coset formulae for a given group $G$ and store this information in the `MackeyContext`, allowing us to check the formulae more quickly for each individual Mackey functor we construct.

**Definition** Let $H\le K\le G$ be subgroups of a group $G$. We say that the tuple $(H,K)$ is a *cover* (and we will denote this by $H\lessdot K$) if $H$ is a proper subgroup of $K$, and if $H\le L \le K$ is a pair of subgroup inclusions, then either $H=L$ or $L=K$.

The MackeyContext data type stores the following information about the group:

- the group itself
- a list of all subgroups of the group
- a list of all the *covers* $H\lessdot K$
- a list of *paths* between two subgroups: for every tuple of subgroups $(H,K)$ with $H\le K$, we store a single sequence of covers $H = H_0 \lessdot H_1 \lessdot \cdots \lessdot H_n = K$
- a list of generators for the group
- a list of relation words among those generators
- a matrix `generator_left_conjugation_matrix`, whose rows are indexed by generators of the group and whose columns are indexed by subgroups, and where the $[g,H]$th entry is the index in the list of subgroups of the subgroup $gHg^{-1}$
- a similar matrix `generator_right_conjugation_matrix` telling us how to find $g^{-1}Hg$ for a generator $g$
- a dictionary `double_coset_info_cache`: the keys are triples of subgroups $(J,H,K)$, interpreted as double cosets $J \backslash H / K$, and the values are vectors of `DoubleCosetInfo` values.  Each entry stores a group-element representative $x$ together with the subgroup indices for both $J^x\cap K$ and $J\cap xKx^{-1}$, since these two conjugate intersections are used in different constructions.

**Proposition**: Let $M$ give the data of $M(H)$ for every subgroup $H\le G$, a conjugation map $c_g \colon M(H) \to M(gHg^{-1})$ for every generator $g\in G$, and restriction and transfer maps for every cover $H\lessdot K$. Then if the double coset formula holds along covers, it holds in general.

```@docs
MackeyContext
```
