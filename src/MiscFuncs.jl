# "haskey" replacement for IdDicts
function find_gap_pair_key(dict, A, B)
    for key in keys(dict)
        (A0, B0) = key
        if A0 === A && B0 == B
            return key
        end
    end

    return nothing
end