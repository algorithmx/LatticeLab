mutable struct HoppingParameter
    UC::LatticeInfo
    HPLN::Dict
end

#NOTE Example:
# Trivial hopping
# Dict( :f1 => 1.0, :f2 => 0.5 )
# Slater-Koster
# A = [[1.0,0.0] [0.0,1.0]]
# Dict( :f1 => Dict((1,2)=>A,(2,1)=>A'), :f2 => ... )

function isequal(hp1::HoppingParameter,hp2::HoppingParameter)
    return (hp1.UC==hp2.UC) && (hp1.HPLN==hp2.HPLN)
end

function ==(hp1::HoppingParameter,hp2::HoppingParameter)
    return isequal(hp1,hp2)
end

@inline is_HoppingParameter(hp) = ( Set(fieldnames(HoppingParameter))==Set(fieldnames(typeof(hp)))
                                && ( is_UnitCell(hp.UC)
                                && ( isa(hp.HPLN, Dict) && keytype(hp.HPLN)<:Symbol ) ) )


copy(hp::HoppingParameter) = HoppingParameter( copy(hp.UC), copy(hp.HPLN) )


check_compat(hp::HoppingParameter) = true #XXX


@inline convert(::Type{T}, hp::T) where {T<:HoppingParameter} = hp


function convert(::Type{T}, hp) where {T<:HoppingParameter}
    #XXX other kinds of LatticeInfo ?
    TUC = is_UnitCell(hp.UC) ? UnitCell : LatticeInfo
    HoppingParameter( convert(TUC, hp.UC), Dict(hp.HPLN)  )
end


#: ----------------- bond force models ------------------
#NOTE the common parameters are Δ::Vector sometimes it is useless...

function hopping_default(
    dim1::Int64, dim2::Int64
    )
    return zeros(Float64, dim1, dim2)
end


function hopping_t_for_all_orbits(
    dim1::Int64, dim2::Int64,
    t::T
    ) where {T<:Number}
    return Matrix( t*I, dim1, dim2)
end

#NOTE CONVENTION 1<--2
# H(1,2) C†(1) C(2)

@inline expiφ_Landau_gauge_0_x(Ri,Rj,Γ) = cis(Γ*((Ri[2]-Rj[2])*Rj[1]+0.5*(Ri[1]-Rj[1])*(Ri[2]-Rj[2])))

@inline expiφ_Landau_gauge_y_0(Ri,Rj,Γ) = cis(-Γ*((Ri[1]-Rj[1])*Rj[2]+0.5*(Ri[1]-Rj[1])*(Ri[2]-Rj[2])))


function hopping_t_for_all_orbits_Hofstadter_expiφ_Landaux(
    Ri::Vector, Rj::Vector, 
    dim1::Int64, dim2::Int64, 
    t::T,
    Γ::S
    ) where {T<:Number, S<:Real}
    return Matrix( (expiφ_Landau_gauge_0_x(Ri,Rj,Γ)*t)*I, dim1, dim2 )
end


function hopping_t_for_all_orbits_Hofstadter_expiφ_Landauy(
    Ri::Vector, Rj::Vector, 
    dim1::Int64, dim2::Int64, 
    t::T,
    Γ::S
    ) where {T<:Number, S<:Real}
    return Matrix( (expiφ_Landau_gauge_y_0(Ri,Rj,Γ)*t)*I, dim1, dim2 )
end


function hopping_matrix_direct_inject(
    dim1::Int64, dim2::Int64,
    mat::T
    ) where {T<:AbstractArray}
    @assert size(mat)==(dim1,dim2)
    return mat
end


# -----------------------------------------------------------
#TODO test Slater-Koster
function hopping_SlaterKoster(
    Δ::Vector, 
    dim1::Int64, dim2::Int64,
    V::Dict{Symbol,T},
    orbits1::Vector{Symbol},
    orbits2::Vector{Symbol}
    ) where {T<:Real}
    lmn = normalize(Δ)
    mat = Float64[ dot( getv(V, SlaterKosterIntegralNames[trim_ud(o1,o2)]),
                        SlaterKosterDict[trim_ud(o1,o2)](lmn...) )
                            for o1 ∈ orbits1, o2 ∈ orbits2 ]
    @assert size(mat)==(dim1,dim2)
    return mat
end

# ----------------

function on_site_potential(
    dim1::Int64,
    mat_or_number
    )
    if mat_or_number isa Number
        return Matrix( mat_or_number*I, dim1, dim1)
    elseif mat_or_number isa AbstractMatrix
        @assert size(mat_or_number)==(dim1,dim1)
        return mat_or_number
    else
        return zeros(Float64, dim1, dim1)
    end
end

# ----------------


function hopping_parameter_matrix(
    type::Symbol,
    ΔRIRJ::Tuple{Array{N,1},Array{N,1},Array{N,1},Int64,Int64},
    params
    ) where { N<:Real }
    (Δ, Ri, Rj, dim1, dim2) = ΔRIRJ
    if type == :ALL
        return hopping_t_for_all_orbits(dim1,dim2,params[1])
    elseif type == :MAT
        return hopping_matrix_direct_inject(dim1,dim2,params[1])
    elseif type == :SK
        return hopping_SlaterKoster(Δ,dim1,dim2,params...)
    elseif type == :HOFALLX
        return hopping_t_for_all_orbits_Hofstadter_expiφ_Landaux(Ri, Rj, dim1, dim2, params...)
    elseif type == :HOFALLY || type == :HOFALL
        return hopping_t_for_all_orbits_Hofstadter_expiφ_Landauy(Ri, Rj, dim1, dim2, params...)
    else
        # the first parameter is always the bond vector
        return hopping_default(dim1,dim2)
    end
end


function zero_HoppingParameter(latt::Lattice)
    A,B,F = findnz(latt.f)
    return HoppingParameter( latt.UC, Dict(k=>(:DEFAULT,(0.0,)) for k ∈ sort(unique(F))) )
end

function zero_onsite_potential(latt::Lattice)
    return Dict(k=>0.0 for k ∈ sort(unique(latt.UC.m)))
end

function one_onsite_potential(latt::Lattice)
    return Dict(k=>1.0 for k ∈ sort(unique(latt.UC.m)))
end

# ------------------------------------------------------------------
