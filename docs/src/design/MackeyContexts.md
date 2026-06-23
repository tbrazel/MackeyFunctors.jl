# Mackey contexts

One of our time-saving design choices is the idea of a [`MackeyContext`](@ref), which precompiles relevant information about a finite group before working with Mackey functors for that group. The idea is to make the construction of a Mackey functor a two-step process: step 1 is to build the context, and then step 2 is to build a Mackey functor for the context. This allows us to avoid certain costly computations, for instance computing what the formulas for the double coset formulas should look like.

**Definition** Let $H\le K\le G$ be subgroups of a group $G$. We say that the tuple $(H,K)$ is a *cover* (and we will denote this by $H\lesssim K$) if $H$ is a proper subgroup of $K$, and if $H\le L \le K$ is a pair of subgroup inclusions, then either $H=L$ or $L=K$.

The MackeyContext data type stores the following information about the group:
- the group itself
- a list of all subgroups of the group
- a list of all the *covers* $H\lessim K$
- a list of *paths* between two subgroups: for every tuple of subgroups $(H,K)$ with $H\le K$, we store a single sequence of covers $H = H_0 \lessim H_1 \lessim \cdots \lessim H_n = K$
- a list of generators for the group
- a matrix `generatorLeftConjugationMatrix`, whose rows are indexed by generators(?) of the group and whose columns are indexed by subgroups, and where the $[g,H]$th entry is the index in the list of subgroups of the subgroup $gHg^{-1}$
- a similar matrix `generatorRightConjugationMatrix` telling us how to find $g^{-1}Hg$ for a generator $g$
- a dictionary `doubleCosetRepresentatives`: the keys are triples of subgroups $(J,H,K)$, with $J\lessim H$ and $K\lessim H$, and the values are a vector of `DoubleCosetFormulaTerm` types: these are tuples of generator words and subgroup indices, indended to help us locate $J^x \cap K$.

