function interatomic_force_matrix(
    type::Symbol,
    Δij::Tuple,
    params
    )
    (Δ::Vector,si::Int64,sj::Int64) = Δij
    if type == :FC3X3
        return params[1]
    end
    if type == :L
        return _L_(Δ,si,sj,params...)
    elseif type == :LT
        return _LT_(Δ,si,sj,params...)
    elseif type == :LTZ
        return _LTZ_(Δ,si,sj,params...)
    elseif type == :LTR
        return _LTR_(Δ,si,sj,params...)
    #elseif type == :BondBending
    #    return bond_bending(Δ,si,sj,params...)
    #elseif type == :StraightBondBending
    #    return straight_bond_bending(Δ,si,sj,params...)
    else
        # the first parameter is always the bond vector
        return default_force(Δ)
    end
end

##: -------------------------------------------------------

@inline dgn(m,v) = LinearAlgebra.Matrix(v*I,m,m)

@inline nXnXT(x) = (one(eltype(x))/dot(x,x)).*(transpose(x).*x)

@inline I_nXnXT(x) = eye(length(x)) .- nXnXT(x)



##: ------------------ bond force models ------------------

#: for both 2D and 3D 

default_force(Δ) = zeros(Float64,length(Δ),length(Δ))

function _L_(
    Δ::Vector{T}, si::Int64, sj::Int64,
    kbs::Real
    ) where {T<:Real}
    (-kbs) .* nXnXT(Δ)
end

function _LZ_(
    Δ::Vector{T1}, si::Int64, sj::Int64,
    nz::Vector{T2},
    kbs::Real, kz::Real
    ) where {T1<:Real,T2<:Real}
    (-kbs) .* nXnXT(Δ) .+ (-kz) .* nXnXT(nz)
end

function _LT_(
    Δ::Vector{T1}, si::Int64, sj::Int64,
    kl::Real, kt::Real
    ) where {T1<:Real}
    #(-kl) .* nXnXT(Δ) .+ (-kt) .* I_nXnXT(Δ)
    (kt-kl) .* nXnXT(Δ) .+ dgn(length(Δ),-kt)
end

function _LTZ_(
    Δ::Vector{T1}, si::Int64, sj::Int64,
    nz::Vector{T2},
    kl::Real, kt::Real, kz::Real
    ) where {T1<:Real,T2<:Real}
    (-kl) .* nXnXT(Δ) .+ (-kt) .* nXnXT(cross3(Δ,nz)) .+ (-kz) .* nXnXT(nz)
end

#: 2D ONLY

function _LTR_(
    Δ::Vector{T1}, si::Int64, sj::Int64,
    kl::Real, kt::Real, kr::Real
    ) where {T1<:Real}
    @assert length(Δ)==2
    n  = one(eltype(Δ))/dot(Δ,Δ)
    Δt = [-Δ[2],Δ[1]]
    return (-n).*(kl.*(transpose(Δ).*Δ) .+ kt.*(transpose(Δt).*Δt) .+ kr.*(transpose(Δt).*Δ.+transpose(Δ).*Δt))
end
