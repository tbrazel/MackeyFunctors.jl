# attempt at pretty printing for Mackey functors

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
    (n,k) = GAP.Globals.IdGroup(obj.context.group)
    name = GAP.Globals.StructureDescription(obj.context.group)
    println(io, "MackeyFunctor for group ", String(name), " (", n, ",", k, ") over base ring ", coefficient_ring(obj))
    
    length(obj.context.covers) <= 10 || return

    for (i,(h,k)) in enumerate(obj.context.covers)
        hname = String(GAP.Globals.StructureDescription(obj.context.subgroups[h]))
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

        println("\n", _center(join([" ", hname, " < ", kname, " "]), wR + wT + gap; spacer = "-"))
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