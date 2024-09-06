abstract type LatticeInfo end
(::Type{T})(x) where {T<:LatticeInfo} = x  # WTF !!!
copy(x::LatticeInfo) = copy(x)

Var = Union{Number,Symbol}
is_Var(c) = isa(c,Var)

Coordinates{TR} = Array{TR,2}
is_Coordinates(c) = ( isa(c,Array) && ndims(c)==2 && eltype(c)<:Real )

Index = Vector{Int64}
is_Index(c) = ( isa(c,Vector) && eltype(c)<:Integer )

Indices = Array{Int64,2}
is_Indices(c) = ( isa(c,Array) && ndims(c)==2 && eltype(c)<:Integer )

Mode = Vector{ComplexF64}
is_Mode(c) = ( isa(c,Vector) && eltype(c)<:Complex )

Modes = Array{ComplexF64,2}
is_Modes(c) = ( isa(c,Array) && ndims(c)==2 && eltype(c)<:Complex )

Masses = Vector{Var}
is_Masses(c) = ( isa(c,Vector) && eltype(c)<:Var )

@inline function position_index(M::Masses)
    Mu = unique(M)
    return Dict(i=>findfirst(Mu.==M[i]) for i ∈ 1:length(M))
end

Orbits = Vector{Symbol}
is_Orbits(o) = ( isa(o,Vector) && eltype(o)<:Symbol )

@inline isequal(OA::Orbits,OB::Orbits) = all(OA.==OB)

# a mode is a "flattened" displacement vector field
check_compat(coo::Coordinates, moo::Modes) = ( size(coo,1)*size(coo,2)==size(moo,1) )
check_compat(moo::Modes, coo::Coordinates) = ( size(coo,1)*size(coo,2)==size(moo,1) )
