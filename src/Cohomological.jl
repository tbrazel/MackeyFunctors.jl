using GAP
using AbstractAlgebra
using MackeyFunctors


##Implement is_cohomological Boolean method
#25

#check each cover_index
#use map_equal
#compose res then transfer in julia (as matrices?)
"""
A Mackey functor M is cohomological if the composition of the transfer and restriction maps between any cover relation K < H
and K is equal to the index [K : H].

"""
##function is_cohomological(context::MackeyContext, H_idx::Int)
        ##subgroups  = context.subgroups
        ##H   = subgroups[H_idx]

        #now define cover restrictions and cover transfers over all indexes of H
        ##cover_restrictions = map(context.covers) do (J_idx, K_idx)
        ##K  = subgroups[K_idx]
        ##end

        ##for (H, (i, j)) in enumerate(context.covers)
            ##cover_restrictions

##end
