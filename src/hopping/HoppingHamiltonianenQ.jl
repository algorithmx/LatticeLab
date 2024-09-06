mutable struct HoppingHamiltonianenQ
    LATT::Lattice
    MATS::Vector
    TS::Vector{HoppingParameter}
    COEFFS::Vector
    VS::Vector
end


@inline is_HoppingHamiltonianenQ(hhqs) = ( Set(fieldnames(HoppingHamiltonianenQ))==Set(fieldnames(typeof(hhqs)))
                                             && all(isa(MAT,Dict) for MAT ∈ hhqs.MATS)
                                             && all(is_HoppingParameter(T) for T ∈ hhqs.TS)
                                             && all(isa(V,Dict) for V ∈ hhqs.VS) )


check_compat(HQS::HoppingHamiltonianenQ) = ( check_compat(HQS.LATT)
                                           && all( size(MAT[r],1)==size(MAT[r],2)==total_num_orbits(HQS.LATT)
                                                   for MAT ∈ HQS.MATS for r ∈ keys(MAT) )
                                           && all(length(r)==dimensions(HQS.LATT) for MAT ∈ HQS.MATS for r ∈ keys(MAT))
                                           && all(issubset(Set(sck(T.HPLN)),Set(sck(HQS.LATT.LN.SPNB))) for T ∈ HQS.TS)
                                           && all(T.UC == HQS.LATT.UC for T ∈ HQS.TS)
                                           && all(check_compat(T) for T ∈ HQS.TS) )


copy(HQS::HoppingHamiltonianenQ) = HoppingHamiltonianenQ( copy(HQS.LATT), copy.(HQS.MATS),
                                                          copy.(HQS.TS), copy(HQS.COEFFS),
                                                          copy.(HQS.VS) )


@inline convert(::Type{T}, x::T) where {T<:HoppingHamiltonianenQ} = x

function convert(::Type{T}, hhq) where {T<:HoppingHamiltonianenQ}
    if is_HoppingHamiltonianQ(hhq)
        return HoppingHamiltonianenQ( convert(Lattice, hhq.LATT),
                                      [hhq.MAT,], [hhq.T,], [1.0,], [hhq.V,] )
    elseif is_HoppingHamiltonianenQ(hhq)
        return HoppingHamiltonianenQ(  convert(Lattice,        hhq.LATT),
                                       convert(Vector{Dict},   hhq.MATS),
                                       convert(Vector{HoppingParameter},    hhq.TS),
                                       convert(Vector{Number}, hhq.COEFFS),
                                       convert(Vector{Dict},   hhq.VS) )
    else
        throw("Typef error in function convert(::Type{T}, hhq) where {T<:HoppingHamiltonianenQ}")
        return hhq
    end
end


@inline consistency_check(hqs::HoppingHamiltonianenQ; TOL=1e-8) = all(consistency_check(MAT,TOL=TOL) for MAT ∈ hqs.MATS)


# favourite function
function HoppingHamiltonian(
    q::Vector,
    HQS::HoppingHamiltonianenQ
    )
    return qMATS(q, HQS.LATT.UC.a, HQS.MATS, HQS.COEFFS)
end

# ----------------------------------------------------------------------------

function *(c::Number, HQS::HoppingHamiltonianenQ)
    HQS1 = copy(HQS)
    HQS1.COEFFS = (c.*copy(HQS1.COEFFS))
    return HQS1
end

function *(c::Number, HQ::HoppingHamiltonianQ)
    return *(c, convert(HoppingHamiltonianenQ,HQ))
end


function +(HQS1::HoppingHamiltonianenQ, HQS2::HoppingHamiltonianenQ)
    @assert HQS1.LATT == HQS2.LATT
    return HoppingHamiltonianenQ(  HQS1.LATT,
                                   vcat(HQS1.MATS,   HQS2.MATS),
                                   vcat(HQS1.TS,     HQS2.TS),
                                   vcat(HQS1.COEFFS, HQS2.COEFFS),
                                   vcat(HQS1.VS, HQS2.VS) )
end

function +(HQS::HoppingHamiltonianenQ,HQ::HoppingHamiltonianQ)
    return +(HQS, convert(HoppingHamiltonianenQ,HQ))
end

function +(HQ::HoppingHamiltonianQ,HQS::HoppingHamiltonianenQ)
    return +(convert(HoppingHamiltonianenQ,HQ), HQS)
end

function +(HQ1::HoppingHamiltonianQ, HQ2::HoppingHamiltonianQ)
    return +(convert(HoppingHamiltonianenQ,HQ1), convert(HoppingHamiltonianenQ,HQ2))
end
