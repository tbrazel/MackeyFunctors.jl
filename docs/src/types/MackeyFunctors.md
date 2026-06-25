# Mackey functors

When asking a user to construct a Mackey functor for a group $G$, it seems reasonable to ask them to supply the values $M(H)$ for each subgroup $H\le G$, but asking for restriction/transfer along *every* subgroup inclusion and conjugation for *every* group element feels a bit excessive. A more concise collection of data is to ask them to supply is:

- the values $M(H)$ for each $H\le G$
- restriction and transfer along covers
- conjugation by generators of $G$

This is a compression of the data of a Mackey functor available in [Webb](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/WebbMF.pdf). Obviously steps are needed to verify that this does indeed satisfy all the axioms of a Mackey functor. Some checks that need to be verified include:

**Check 1**: If $h\in H$, then the conjugation map $c_h$ is the identity on $M(H)$.

**Check 2**: Conjugation by generators of $G$ respect the relations in $G$.

**Check 3**: Conjugation by generators commutes with restrictions and transfers along covers.

**Check 4**: For any two factorizations of $H \le K$ into covers, the supplied cover restriction maps and cover transfer maps yield a well-defined restriction/transfer between $M(H)$ and $M(K)$.

**Check 5**: The double coset formula holds for triples $(J,H,K)$ where $J\lesssim H$ and $K\lesssim H$ are covers.

**Theorem**: The data supplied above, satisfying all three checks, uniquely specifies a $G$-Mackey functor.

```@docs
MackeyFunctor
```

## Extracting data from a Mackey functor
Given a Mackey functor, we can access various parts of the data.

```@docs
coefficient_ring
transfer
restriction
conjugation
value
```