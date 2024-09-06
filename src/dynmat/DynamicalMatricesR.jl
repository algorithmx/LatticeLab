mutable struct DynamicalMatricesR
    LATT::Lattice
    DS::Vector
    EIG::EigenModes
    SPS::Vector            # spring constants Dict(:f1=>1.0, :f2=>1.2)
    COEFFS::Vector
    MS::Vector             # mass Dict(:mA=>1.0, :mB=>1.2)
end

@inline is_DynamicalMatricesR(dmr) = (   Set(fieldnames(DynamicalMatricesR))==Set(fieldnames(typeof(dmr)))
                                        && is_Lattice(dmr.LATT)
                                        && eltype(dmr.DS    )<:AbstractMatrix
                                        && eltype(dmr.SPS   )<:Dict
                                        && eltype(dmr.COEFFS)<:Number
                                        && eltype(dmr.MS    )<:Dict   )

check_compat(DRS::DynamicalMatricesR) = ( check_compat(DRS.LATT)
                                        && all( size(D,1)==size(D,2)==dimensions(DRS.LATT)*num_inner_sites(DRS.LATT)
                                                for D ∈ DRS.DS )
                                        && all(issubset(Set(collect(keys(SP))),Set(keys(DRS.LATT.LN.SPNB))) for SP ∈ DRS.SPS)
                                        && all(issubset(Set(collect(keys(M))),Set(DRS.LATT.UC.m)) for M ∈ DRS.MS) )

copy(DRS::DynamicalMatricesR) = DynamicalMatricesR( copy(DRS.LATT),
                                                    copy.(DRS.DS), copy(DRS.EIG), copy.(DRS.SPS), copy(DRS.COEFFS),
                                                    copy.(DRS.MS) )

@inline convert(::Type{T}, x::T) where {T<:DynamicalMatricesR} = x

function convert(::Type{T}, dmr) where {T<:DynamicalMatricesR}
    if is_DynamicalMatrixR(dmr)
        return DynamicalMatricesR(  convert(Lattice,dmr.LATT),
                                    [dmr.D,], convert(EigenModes,dmr.EIG),
                                    [dmr.SP,], [1.0,], [dmr.M,]  )

    elseif is_DynamicalMatricesR(dmr)
        return DynamicalMatricesR(  convert(Lattice,dmr.LATT),
                                    [convert(AbstractMatrix,D) for D ∈ dmr.DS],
                                    convert(EigenModes, dmr.EIG),
                                    [convert(Dict,SP) for SP ∈ dmr.SPS],
                                    convert(Vector{Number},dmr.COEFFS),
                                    [convert(Dict,M) for M ∈ dmr.MS]  )
    else
        throw("Type error in function convert(::Type{T}, dmr) where {T<:DynamicalMatricesR}")
    end
end


## solve

function solve_eigenmodes(
    dmr::DynamicalMatricesR;
    EPS = 1e-10,
    test = true
    ) where { T<:AbstractMatrix }
    sqrt_eigen( x ) = sqrt(max(0.0,real(x)))
    test_eigen_DR( e ) = all_non_negative_real_v1( e, EPS ) || ( ! test )
    solve_eigenmodes( sum(coeff.*ds for (coeff,ds) ∈ zip(dmr.COEFFS,dmr.DS)),
                      EPS=EPS, eigen_test=test_eigen_DR, eigen_process=sqrt_eigen )
end


## arithmetic opreations

function *(c::Number, DRS::DynamicalMatricesR)
    DRS1 = copy(DRS)
    DRS1.COEFFS = (c.*copy(DRS1.COEFFS))
    return DRS1
end


function *(c::Number, DR::DynamicalMatrixR)
    return *(c, convert(DynamicalMatricesR,DR))
end


function +(DRS1::DynamicalMatricesR, DRS2::DynamicalMatricesR)
    @assert DRS1.LATT == DRS2.LATT
    return DynamicalMatricesR( DRS1.LATT,
                               vcat( DRS1.DS,     DRS2.DS     ),
                               DRS1.EIG,
                               vcat( DRS1.SPS,    DRS2.SPS    ),
                               vcat( DRS1.COEFFS, DRS2.COEFFS ),
                               vcat( DRS1.MS,     DRS2.MS     ) )
end


function +(DRS::DynamicalMatricesR, DR::DynamicalMatrixR)
    return +(DRS, convert(DynamicalMatricesR,DR))
end


function +(DR::DynamicalMatrixR, DRS::DynamicalMatricesR)
    return +(convert(DynamicalMatricesR,DR), DRS)
end


function +(DR1::DynamicalMatrixR, DR2::DynamicalMatrixR)
    return +(convert(DynamicalMatricesR,DR1), convert(DynamicalMatricesR,DR2))
end
