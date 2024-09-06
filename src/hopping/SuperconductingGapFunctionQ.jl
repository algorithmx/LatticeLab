
mutable struct SuperconductingGapFunctionQ
    LATT::Lattice
    MAT::Dict
end

mutable struct SuperconductingGapFunctionenQ
    LATT::Lattice
    MATS::Vector
    COEFFS::Vector
end

# ------------------------------------------------------------

copy(sg::SuperconductingGapFunctionQ) = SuperconductingGapFunctionQ(copy(sg.LATT),copy(sg.MAT))

copy(sg::SuperconductingGapFunctionenQ) = SuperconductingGapFunctionenQ(copy(sg.LATT),copy.(sg.MATS),copy(sg.COEFFS))

@inline is_SuperconductingGapFunction(dq) = (
    Set(fieldnames(SuperconductingGapFunctionQ))==Set(fieldnames(typeof(dq)))
    && is_Lattice( dq.LATT ) && isa( dq.MAT, Dict ) )

@inline is_SuperconductingGapFunctionen(dq) = (
    Set(fieldnames(SuperconductingGapFunctionenQ))==Set(fieldnames(typeof(dq)))
    && is_Lattice( dq.LATT ) && all(isa(M,Dict) for M ∈ dq.MATS) && all(isa(c,Number) for c ∈ dq.COEFFS) )


check_compat(dq::SuperconductingGapFunctionQ) = true #TODO

check_compat(dq::SuperconductingGapFunctionenQ) = true #TODO


@inline consistency_check(dq::SuperconductingGapFunctionQ; TOL=1e-8) = consistency_check(dq.MAT,TOL=TOL)

@inline consistency_check(dq::SuperconductingGapFunctionenQ; TOL=1e-8) = all(consistency_check(MAT,TOL=TOL) for MAT ∈ dq.MATS)


@inline convert(::Type{T}, x::T) where {T<:SuperconductingGapFunctionQ} = x
@inline convert(::Type{T}, x::T) where {T<:SuperconductingGapFunctionenQ} = x


function convert(::Type{T}, dq) where {T<:SuperconductingGapFunctionQ}
    @assert is_SuperconductingGapFunctionQ(dq)
    return SuperconductingGapFunctionQ( convert( Lattice, dq.LATT ),
                                        convert( Dict, dq.MAT ) )
end

function convert(::Type{T}, dq) where {T<:SuperconductingGapFunctionenQ}
    @assert is_SuperconductingGapFunctionenQ(dq)
    TN = typeof(sum(dq.COEFFS))
    return SuperconductingGapFunctionQ( convert( Lattice, dq.LATT ),
                                        convert( Vector{Dict}, dq.MATS ),
                                        convert( Vector{TN}, dq.COEFFS ) )
end

# ------------------------------------------------------------

function SuperconductingGapFunction(
    q::Vector,
    DQ::SuperconductingGapFunctionQ
    )
    return qMAT(q, DQ.MAT, DQ.LATT.UC.a)
end

function SuperconductingGapFunction(
    q::Vector,
    DQS::SuperconductingGapFunctionenQ
    )
    return qMATS(q, DQS.LATT.UC.a, DQS.MATS, DQS.COEFFS)
end


function kspace_s_wave_Δ(latt::Lattice)
    return SuperconductingGapFunctionQ( latt, Dict((0 ⨰ dimensions(latt))=>speye(total_num_orbits(latt))) )
end


# ------------------------------------------------------------

function +(D1::SuperconductingGapFunctionQ,D2::SuperconductingGapFunctionQ)
    @assert D1.LATT == D2.LATT
    return SuperconductingGapFunctionenQ(D1.LATT,[D1.MAT,D2.MAT],[one(ComplexF64),one(ComplexF64)])
end

function +(D1::SuperconductingGapFunctionQ,DS::SuperconductingGapFunctionenQ)
    @assert D1.LATT == DS.LATT
    return SuperconductingGapFunctionenQ( DS.LATT,
                                          vcat([D1.MAT,],DS.MATS),
                                          vcat([one(ComplexF64),],DS.COEFFS) )
end

function +(DS::SuperconductingGapFunctionenQ,D2::SuperconductingGapFunctionQ)
    return +(D2,DS)
end

function +(DS1::SuperconductingGapFunctionenQ,DS2::SuperconductingGapFunctionenQ)
    @assert DS1.LATT == DS2.LATT
    return SuperconductingGapFunctionenQ( DS1.LATT,
                                          vcat(DS1.MATS,DS2.MATS),
                                          vcat(DS1.COEFFS,DS2.COEFFS) )
end

function *(c::N,D::SuperconductingGapFunctionQ) where {N<:Number}
    return SuperconductingGapFunctionenQ(D.LATT,[D.MAT,],[c*one(ComplexF64),])
end

function *(c::N,DS::SuperconductingGapFunctionenQ) where {N<:Number}
    return SuperconductingGapFunctionenQ(DS.LATT,DS.MATS,c.*DS.COEFFS)
end
