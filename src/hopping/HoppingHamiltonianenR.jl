mutable struct HoppingHamiltonianenR
    LATT::Lattice
    HS::Vector
    EIG::EigenModes
    TS::Vector{HoppingParameter}
    VS::Vector
    COEFFS::Vector
end

copy(HRS::HoppingHamiltonianenR) = HoppingHamiltonianenR(   copy(HRS.LATT), copy.(HRS.HS),
                                                            copy(HRS.EIG),  copy.(HRS.TS),
                                                            copy.(HRS.VS),  copy(HRS.COEFFS)  )


@inline is_HoppingHamiltonianenR(hhr) = (  Set(fieldnames(HoppingHamiltonianenR))==Set(fieldnames(typeof(hhr)))
                                          && is_Lattice(hhr.LATT)
                                          && eltypes(hhr.HS)<:AbstractMatrix
                                          && is_EigenModes(hhr.EIG)
                                          && (length(hhr.TS)==0 || is_HoppingParameter(first(hhr.TS)))
                                          && eltype(hhr.VS)<:Dict
                                          && eltype(hhr.COEFFS)<:Number )


check_compat(HR::HoppingHamiltonianenR) = ( check_compat(HR.LATT)
                                         && check_compat(HR.EIG)
                                         && all(check_compat(T) for T ∈ HR.TS)
                                         && all(size(H,1)==size(H,2)==length(HR.EIG.ω) for H ∈ HR.HS) )


@inline convert(::Type{T}, x::T) where {T<:HoppingHamiltonianenR} = x


function convert(::Type{T}, hhr) where {T<:HoppingHamiltonianenR}
    if is_HoppingHamiltonianR(hhr)
        return HoppingHamiltonianenR( convert(Lattice,            hhr.LATT  ),
                                      [convert(AbstractMatrix,hhr.H),],
                                      convert(EigenModes,         hhr.EIG   ),
                                      [convert(HoppingParameter,hhr.T),],
                                      [convert(Dict,hhr.V),],
                                      [1.0,] )
    elseif is_HoppingHamiltonianenR(hhr)
        return HoppingHamiltonianenR( convert(Lattice,            hhr.LATT  ),
                                      [convert(AbstractMatrix,H) for H ∈ hhr.HS],
                                      convert(EigenModes,         hhr.EIG   ),
                                      [convert(HoppingParameter,T) for T ∈ hhr.TS],
                                      [convert(Dict,V) for V ∈ hhr.VS],
                                      [convert(Number,c) for c ∈ hhr.COEFFS] )
    else
        throw("Type error in function convert(::Type{T}, hhr) where {T<:HoppingHamiltonianenR}")
    end
end

# ----------------------------------------------------------------------------

function *(c::Number, HRS::HoppingHamiltonianenR)
    HRS1 = copy(HRS)
    HRS1.COEFFS = (c.*copy(HRS1.COEFFS))
    return HRS1
end

function *(c::Number, HR::HoppingHamiltonianR)
    return *(c, convert(HoppingHamiltonianenR,HR))
end


function +(HRS1::HoppingHamiltonianenR, HRS2::HoppingHamiltonianenR)
    @assert HRS1.LATT == HRS2.LATT
    return HoppingHamiltonianenR(  HRS1.LATT,
                                   vcat(HRS1.HS,     HRS2.HS),
                                   HRS1.EIG,
                                   vcat(HRS1.TS,     HRS2.TS),
                                   vcat(HRS1.VS,     HRS2.VS),
                                   vcat(HRS1.COEFFS, HRS2.COEFFS) )
end

function +(HRS::HoppingHamiltonianenR,HR::HoppingHamiltonianR)
    return +(HRS, convert(HoppingHamiltonianenR,HR))
end

function +(HR::HoppingHamiltonianR,HRS::HoppingHamiltonianenR)
    return +(convert(HoppingHamiltonianenR,HR), HRS)
end

function +(HR1::HoppingHamiltonianR, HR2::HoppingHamiltonianR)
    return +(convert(HoppingHamiltonianenR,HR1),convert(HoppingHamiltonianenR,HR2))
end


# -------------------------------------------------------------------------

function solve_eigenmodes!(HRS::HoppingHamiltonianenR; EPS=1e-8, test=true)
    @assert length(HRS.COEFFS) == length(HRS.HS)
    H = sum(HRS.COEFFS[i].*HRS.HS[i] for i=1:length(HRS.COEFFS))
    test_eigen_HR( e ) = all_real_v1( e, EPS ) || ( ! test )
    HRS.EIG = solve_eigenmodes( H, EPS=EPS, eigen_test=test_eigen_HR, eigen_process=(x->x) )
    return nothing
end
