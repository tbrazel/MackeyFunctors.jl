# attempt at pretty printing for Mackey functors

# produces string of group description with SmallGroup ID #
function _desc_id(G::Group)
    GAP.Globals.Size(G) == 1 && return "e"

    (n,k) = GAP.Globals.IdGroup(G)
    name = string(GAP.Globals.StructureDescription(G))
    return join([name, " (", string(n), ",", string(k), ")"])
end

# produces a string of length w which contains s as a substring in the middle 
function _center(s::AbstractString, w::Int; spacer = " ")
    pad = max(w - length(s), 0)
    left = pad ÷ 2
    right = pad - left
    string(repeat(spacer, left), s, repeat(spacer, right))
end

# converts matrix to list of strings
function _matrix_lines_textplain(A)
    sprint(io -> begin
        ctx = IOContext(io, :limit => false)
        show(ctx, "text/plain", A)
    end) |> x -> split(x, '\n')
end

function Base.show(io::IO, obj::MackeyFunctor)
    println(io, "MackeyFunctor for group ", _desc_id(obj.context.group), " over base ring ", coefficient_ring(obj))
    
    # TODO: find way to toggle this limit?
    # length(obj.context.covers) <= 10 || return

    for (i,(h,k)) in enumerate(obj.context.covers)
        H = obj.context.subgroups[h]
        K = obj.context.subgroups[k]
        kname = String(GAP.Globals.StructureDescription(obj.context.subgroups[k]))
        
        # restriction and transfer matrices
        R = matrix(obj.cover_restrictions[i])
        T = matrix(obj.cover_transfers[i])

        # lines of matrices and their widths
        linesR = _matrix_lines_textplain(R)
        linesT = _matrix_lines_textplain(T)
        wR = maximum(length, linesR; init=0)
        wT = maximum(length, linesT; init=0)

        # space between matrices
        gap = 4
        spacer = repeat(" ", gap)

        println("\n", _center(join([" ", _desc_id(H), " < ", _desc_id(K), " "]), wR + wT + gap; spacer = "-"))
        println(io, _center("res", wR), spacer, _center("tr", wT))

        # print lines of matrices one by one, padding if one matrix runs out of lines
        n = max(length(linesR), length(linesT))
        for i in 1:n
            r = i <= length(linesR) ? linesR[i] : ""
            t = i <= length(linesT) ? linesT[i] : ""
            println(io, rpad(r, wR), spacer, rpad(t, wT))
        end
    end
end