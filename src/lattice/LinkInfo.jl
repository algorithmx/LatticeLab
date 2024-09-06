# direction
Direction = Vector


"""
    Spring = Dict{Tuple{Int64,Int64},Vector{Direction}}

Type `String` is a `Dict` with keys of sublattice index pairs (a,b) 
and values of bond directions point from a to b. For a pair of sublattice, 
there can be multiple directions.

# Examples

```
Dict( (1, 1) => [[0.0, 0.0, 1.0],] )
Dict( (subl_1,subl_2) => [δ_1, δ_1, ..., δ_N] )
```
"""
Spring = Dict{Tuple{Int64,Int64},Vector{Direction}}



is_Spring(sp) = ( keytype(sp)<:Tuple{Int64,Int64}
               && valtype(sp)<:Vector
               && (eltype(values(sp))<:Vector || eltype(values(sp))<:Direction) )


copy(sp::Spring) = Dict(copy(k)=>copy(sp[k]) for k ∈ keys(sp))

check_compat( sp::Spring,
              Nsubl::Int64) = all([ all(t.>0) && all(t.<=Nsubl) for t ∈ keys(sp) ])

function isequal(SP1::Spring,SP2::Spring)
    return SP1==SP2
end




"""
    mutable struct LinkInfo
        UC::LatticeInfo
        SPNB::Dict
    end

# Note 

LinkInfo.SPNB is the `Dict` with spring constant label as keys 
and `Spring` as values, i.e. 

```LinkInfo.SPNB[:f1] = Dict((1,2)=>[[0.5,0.5,0.0],], (2,1)=>[[-0.5,-0.5,0.0],],)```


"""
mutable struct LinkInfo
    UC::LatticeInfo
    SPNB::Dict{Symbol,Spring}
end



function isequal(LN1::LinkInfo,LN2::LinkInfo)
    return (LN1.UC==LN2.UC) && (LN1.SPNB==LN2.SPNB)
end

function ==(LN1::LinkInfo,LN2::LinkInfo)
    return isequal(LN1,LN2)
end

@inline is_LinkInfo(fc) = ( Set(fieldnames(LinkInfo))==Set(fieldnames(typeof(fc)))
                              && ( is_UnitCell(fc.UC) || false ) #XXX other kinds of LatticeInfo ?
                              && isa(fc.SPNB, Dict) && valtype(fc.SPNB)<:Spring )

copy(fc::LinkInfo) = LinkInfo( copy(fc.UC), copy(fc.SPNB) )

check_compat(fc::LinkInfo) = ( all([check_compat(nb,fc.UC.nsubl) for nb ∈ values(fc.SPNB)])
                                 && check_compat(fc.UC) )


@inline convert(::Type{T}, fc::T) where {T<:LinkInfo} = fc

function convert(::Type{T}, fc) where {T<:LinkInfo}
    @assert is_LinkInfo(fc)
    #+ other kinds of LatticeInfo ?
    TUC = is_UnitCell(fc.UC) ? UnitCell : LatticeInfo
    KSP = keytype(fc.SPNB)
    LinkInfo( convert( TUC, fc.UC ), Dict{KSP,Spring}( fc.SPNB ) )
end

