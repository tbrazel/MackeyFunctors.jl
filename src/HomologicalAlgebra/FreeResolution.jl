function free_resolution(
    mf::MackeyFunctor,
    max_deg::Int
)::Vector{MackeyFunctorHomomorphism}
    # begin resolution with epimorphism onto mf from free
    d = MackeyFunctorHomomorphism[epimorphism_from_free(mf)]

    # extend resolution out to max_deg
    for n in 1:max_deg 
        K, i = kernel(d[end])
        push!(d,compose(epimorphism_from_free(K),i))
    end

    return d
end 