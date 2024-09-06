#NOTE meaning of the components
#(1) origin, usually slightly different from [0,0,0] to guarantee 
#    all sites to be in the unit cell
#(2) translation unit to define the bulk, in terms of Bavais basis
#    example t1=1*a1+2*a2, t2=2*a1+1*a2, trans = [[1,2] [2,1]]
#(3) unit of translation in each side of the bulk, counted by t1,t2
#    example: shift[:,1] = n*t1+m*t2

"""
    BoundingBox = Tuple{Vector,Matrix{Int64},Vector{Int64},Vector{Bool}}


    (origin, trans, shift_units, pbc) = bbox
"""
BoundingBox = Tuple{Vector,Matrix{Int64},Vector{Int64},Vector{Bool}}

is_BoundingBox(bb) = ( isa(bb,Tuple)
                    && length(bb)==4
                    && isa(bb[1],Vector)
                    && isa(bb[2],Matrix) && eltype(bb[2])<:Integer
                    && isa(bb[3],Vector) && eltype(bb[3])<:Integer
                    && ((isa(bb[4],Vector) && eltype(bb[4])<:Bool) || isa(bb[4],BitArray{1})) )

copy(bb::BoundingBox) = BoundingBox((copy(bb[1]),copy(bb[2]),copy(bb[3]),copy(bb[4])))

check_compat(bb::BoundingBox) = ( length(bb[1])==size(bb[2],1)==size(bb[2],2)==length(bb[3])==length(bb[4]) )

@inline convert(::Type{T}, bbx::T) where {T<:BoundingBox} = bbx

function convert(::Type{T}, bbx::TB) where {T<:BoundingBox, TB<:Tuple}
    @assert is_BoundingBox(bbx)
    BoundingBox( Vector(bbx[1]),
                 Matrix{Int64}(bbx[2]),
                 Vector{Int64}(bbx[3]),
                 Vector{Bool}(bbx[4]) )
end


@inline function construct_parallelgon(origin,shift)
    return [ vec(origin).+vec(shift*collect(t))  # Cartesian
             for t ∈ Iterators.product([0:1 for k=1:size(shift,1)]...) ]
end


#NOTE bbox[2] is the basis of a supercell.
# It should be identity for conventional construction,
# but can also be costomized to construct special shape of bulk
@inline function origin_shift(bbox::BoundingBox, a)
    return ( bbox[1], a*(bbox[2]*diagm(0=>bbox[3])) )
end


# (origin,trans,units,pbc) = bbox
#p1 = inv(basis*(trans*diagm(0=>units)))*(point.-origin)

@inline inbbox_old(point_wrt_origin,invS) =  all(x->(x<1.0 && x>=0.0),(invS*point_wrt_origin))

#+  MORE OPTIMIZATION !!!
function inbbox(point_wrt_origin,invS)
    A = invS*point_wrt_origin
    if any([A[1]<0.0, A[1]>=1.0, A[2]<0.0, A[2]>=1.0])
        return false
    else
        return (length(A)<=2) || (A[3]>=0.0 && A[3]<1.0)
    end
end


@inline inbboxt(point_wrt_origin,invS,t) = length(t)==0 || (all((invS*point_wrt_origin)[t].<1.0) && all((invS*point_wrt_origin)[t].>=0.0))
@inline inbbox3(point_wrt_origin,invS) = ( ((invS*point_wrt_origin)[3]<1.0) && ((invS*point_wrt_origin)[3]>=0.0) )


@inline Vceil(Av)  = Int64.( ceil.(Av))
@inline Vfloor(Av) = Int64.(floor.(Av))


function bounding_box_Nmin_Nmax_old_version(
    basis,
    bbox::BoundingBox,
    margin
    )::Tuple
    (origin, shift) = origin_shift(bbox,basis)
    verts = construct_parallelgon(origin,shift) # Cartesian
    invA  = inv(basis)
    oNf   = hcat(map(x->Vfloor(invA*x), verts)...) # frac
    oNc   = hcat(map(x-> Vceil(invA*x), verts)...) # frac
    oN1   = minimum(hcat(minimum(oNf,dims=2)|>vec,minimum(oNc,dims=2)|>vec),dims=2) |> vec
    oN2   = maximum(hcat(maximum(oNf,dims=2)|>vec,maximum(oNc,dims=2)|>vec),dims=2) |> vec
    #: it was (oN1.-margin,oN2.+margin) without +1
    return (oN1.-(margin), oN2.+(margin+1))
end


function bounding_box_Nmin_Nmax(
    basis,
    bbox::BoundingBox,
    margin
    )::Tuple
    (origin, shift) = origin_shift(bbox,basis)
    verts = construct_parallelgon(origin,shift) # Cartesian
    invA  = inv(basis)
    oNf   = hcat(map(x->Vfloor(invA*x), verts)...) # frac
    oNc   = hcat(map(x-> Vceil(invA*x), verts)...) # frac
    oN1   = minimum(hcat(minimum(oNf,dims=2)|>vec,minimum(oNc,dims=2)|>vec),dims=2) |> vec
    oN2   = maximum(hcat(maximum(oNf,dims=2)|>vec,maximum(oNc,dims=2)|>vec),dims=2) |> vec
    #: it was (oN1.-margin,oN2.+margin) without +1
    pbc = bbox[4] #% new 
    if length(pbc)==3 
        if pbc[3]
            return (oN1.-(margin), oN2.+(margin+1))
        else # slab perpendicular to z-axis
            return ([(oN1[1:2].-(margin))..., 0], [(oN2[1:2].+(margin+1))..., 0]) #% new 
        end
    else
        return ([(oN1[1:2].-(margin))...,], [(oN2[1:2].+(margin+1))...,]) #% new 
    end
end
