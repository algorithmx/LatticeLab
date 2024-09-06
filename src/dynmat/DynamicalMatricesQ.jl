mutable struct DynamicalMatricesQ
    LATT::Lattice
    MATS::Vector       # mn => Φ
    SPS::Vector        # spring constants Dict(:f1=>1.0, :f2=>1.2)
    COEFFS::Vector
    M::Dict            # mass Dict(:mA=>1.0, :mB=>1.2)
end

@inline is_DynamicalMatricesQ(dmq) = (   Set(fieldnames(DynamicalMatricesQ))==Set(fieldnames(typeof(dmq)))
                                        && is_Lattice(dmq.LATT)
                                        && eltype(dmq.MATS  )<:Dict
                                        && eltype(dmq.SPS   )<:Dict
                                        && eltype(dmq.COEFFS)<:Number
                                        && isa(dmq.M,   Dict) )

check_compat(DQS::DynamicalMatricesQ) = ( check_compat(DQS.LATT)
                                        && all(length(r)==dimensions(DQS.LATT) for MAT ∈ DQS.MATS for r ∈ keys(MAT))
                                        && all( size(MAT[r],1)==size(MAT[r],2)==dimensions(DQS.LATT)*num_sublattice(DQS.LATT)
                                                for MAT ∈ DQS.MATS for r ∈ keys(MAT) )
                                        && all(issubset(Set(collect(keys(SP))),Set(keys(DQS.LATT.LN.SPNB))) for SP ∈ DQS.SPS)
                                        && issubset(Set(collect(keys(DQS.M))),Set(DQS.LATT.UC.m)) )

copy(DQS::DynamicalMatricesQ) = DynamicalMatricesQ( copy(DQS.LATT),
                                                    copy.(DQS.MATS), copy.(DQS.SPS), copy(DQS.COEFFS),
                                                    copy(DQS.M) )

@inline convert(::Type{T}, x::T) where {T<:DynamicalMatricesQ} = x

function convert(::Type{T}, dmq) where {T<:DynamicalMatricesQ}
    if is_DynamicalMatrixQ(dmq)
        return DynamicalMatricesQ( convert(Lattice,dmq.LATT),
                                    [dmq.MAT,], [dmq.SP,], [1.0,], convert(Dict,dmq.M) )
    elseif is_DynamicalMatricesQ(dmq)
        return DynamicalMatricesQ( convert(Lattice,dmq.LATT),
                                    convert(Vector{Dict},dmq.MATS),
                                    convert(Vector{Dict},dmq.SPS),
                                    convert(Vector{Number},dmq.COEFFS),
                                    convert(Dict,dmq.M) )
    else
        throw("Type error in function convert(::Type{T}, dmq) where {T<:DynamicalMatricesQ}")
    end
end

@inline consistency_check(DQS::DynamicalMatricesQ; TOL=1e-8) = all(consistency_check(MAT,TOL=TOL) for MAT ∈ DQS.MATS)

# favourite function
function DynamicalMatrix(
    q,
    DQS::DynamicalMatricesQ
    )
    # the Fourier sum
    return qMATS(q, DQS.LATT.UC.a, DQS.MATS, DQS.COEFFS)
end


# ----------------------------------------------------------------------------

function *(c::Number, DQS::DynamicalMatricesQ)
    DQS1 = copy(DQS)
    DQS1.COEFFS = (c.*copy(DQS1.COEFFS))
    return DQS1
end

function *(c::Number, DQ::DynamicalMatrixQ)
    return *(c, convert(DynamicalMatricesQ,DQ))
end


function +(DQS1::DynamicalMatricesQ, DQS2::DynamicalMatricesQ)
    @assert DQS1.M == DQS2.M
    @assert DQS1.LATT == DQS2.LATT
    return DynamicalMatricesQ( DQS1.LATT,
                               vcat(DQS1.MATS,   DQS2.MATS),
                               vcat(DQS1.SPS,    DQS2.SPS),
                               vcat(DQS1.COEFFS, DQS2.COEFFS),
                               DQS1.M )
end


function +(DQS::DynamicalMatricesQ,DQ::DynamicalMatrixQ)
    return +(DQS, convert(DynamicalMatricesQ,DQ))
end

function +(DQ::DynamicalMatrixQ,DQS::DynamicalMatricesQ)
    return +(convert(DynamicalMatricesQ,DQ), DQS)
end

function +(DQ1::DynamicalMatrixQ, DQ2::DynamicalMatrixQ)
    return +(convert(DynamicalMatricesQ,DQ1),convert(DynamicalMatricesQ,DQ2))
end
