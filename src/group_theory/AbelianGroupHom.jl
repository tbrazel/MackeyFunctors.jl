# Check if a matrix defines a finitely generated abelian group homomorphism
function defines_fg_abelian_map(source_group::Vector{Int}, target_group::Vector{Int}, M::AbstractMatrix{<:Integer})
    # Number of rows and columns of the matrix M
    num_rows_of_M = size(M,1)
    num_cols_of_M = size(M,2)

    # Verify everything matches correctly
    length(source_group) == num_cols_of_M || throw(DimensionMismatch("Source group doesn't match the number of rows of the matrix"))
    length(target_group) == num_rows_of_M || throw(DimensionMismatch("Target group doesn't match the number of rows of the matrix"))

    # We now have to check the torsion relations hold
    for j in eachindex(source_group)

        # In the jth spot, the domain abelian group is Z/n
        n = source_group[j]
        
        # If n=0, we are mapping out of Z, so there is nothing to check
        n == 0 && continue   # source generator has infinite order

        # For each entry in the target
        for i in eachindex(target_group)

            # The ith spot is Z/m
            m = target_group[i]

            # The (i,j)th entry in M is a map a : Z/n -> Z/m
            a = M[i, j]

            # If we have a map Z/n -> Z
            if m == 0
                # Then we need a =0 or n=0
                n * a == 0 || return false
            
            # Otherwise we have some map Z/n -> Z/m
            else
                # And this only makes sense if a*n is 0 mod m
                mod(n * a, m) == 0 || return false
            end
        end
    end
    return true
end

# Check if two integral matrices define identical maps between fg abelian groups
function are_equal_abelian_group_homomorphisms(source_group::Vector{Int}, target_group::Vector{Int}, M::AbstractMatrix{<:Integer},N::AbstractMatrix{<:Integer})
    # First check the two inputted maps are each well-defined
    defines_fg_abelian_map(source_group,target_group,M) || throw(DimensionMismatch("First inputted matrix is not a well-defined map between the two groups"))
    defines_fg_abelian_map(source_group,target_group,N) || throw(DimensionMismatch("Second inputted matrix is not a well-defined map between the two groups"))

    for j in eachindex(source_group)
        # In the jth spot, the domain abelian group is Z/n
        n = source_group[j]
        # For each entry in the target
        for i in eachindex(target_group)
            # The ith spot is Z/m
            m = target_group[i]
            if m==0
                a = M[i,j]
                b = N[i,j]
                a == b || return false
            else
                a = M[i,j]
                b = N[i,j]
                mod(a,m) == mod(b,m) || return false
            end
        end
    end
    return true
end

# Checks if M:group->group is the identity
function is_identity_abelian_group_homomorphism(group::Vector{Int},M::AbstractMatrix{<:Integer})
    defines_fg_abelian_map(group,group,M) || throw(DimensionMismatch("Matrix is not a well-defined endomorphism of the group"))

    for j in eachindex(group)
        for i in eachindex(group)
            if i!=j
                mod(M[i,j],group[j]) == 0 || return false
            else
                mod(M[i,j],group[j]) == 1 || return false
            end
        end
    end
    return true
end
