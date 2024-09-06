#NOTE has been updated to work for dimension = 1,2,3,4,...

##

using LinearAlgebra
using SparseArrays

@inline function delete_empty(
    dic::Dict{K,V}
    ) where { K<:Any, V<:SparseMatrixCSC }
    return Dict( k=>dic[k] for k ∈ keys(dic)
                            if length(dic[k].nzval)>0 && sum(abs.(dic[k]))>1e-20 )
end

@inline function delete_empty(
    dic::Dict{K,V}
    ) where { K<:Any, V<:Vector }
    return Dict( k=>unique(v) for (k,v) ∈ dic if length(v)>0 )
end


# ----------------------------------------------------------
_GLOBAL_MINIMAL_PARTITION_SIZE_ = 8

#NOTE find the min/max box for each rows of rr
@inline function get_BOX(rr)
    return map( x->[x[1],x[2]],
                zip(vec(minimum(rr,dims=2)),vec(maximum(rr,dims=2))) )
end


#NOTE meaning of radius :
# the radius of the maximal sphere that
# covers each region in the bisection of R0
# R(i) = (x(i)_max - x(i)_min)/4
# R = √ ∑ᵢ R(i)^2
@inline radius_bisection(bx) = sqrt( sum( map(x->(1.05*(x[2]-x[1])/4)^2, bx) ) )


#NOTE if the "thickness" of the point set along a direction
#> is less than ϵ0 : take center
#> else : take 1/4 and 3/4 (bisect)
@inline function bisect_each_dimension(
    x::Vector,
    ϵ0::Float64
    )
    return ((abs(x[2]-x[1])<ϵ0) ? [[0.5,0.5]' *x,] : [[0.75,0.25]' *x, [0.25,0.75]' *x])
end


#NOTE bisect the region : return centers of the regions
# the argument `box` is the result of `get_BOX(rr)`
function get_origin(
    box::Vector,
    ϵ0::Float64
    )
    return map(collect, Iterators.product(map(x->bisect_each_dimension(x,ϵ0),box)...))
end


#NOTE index_R0 is the index for R0
# in the collection of all lattice positions R0
function refine(
    R0::Array{T,2},
    index_R0::Vector{Int64},
    ϵ0::Float64
    )::Dict  where { T<:Real }

    (dim,Nsites) = size(R0)
    box = get_BOX(R0)

    origins = get_origin(box,ϵ0)
    radius = radius_bisection(box)

    @inline pick(O,R) = [index_R0[i] for i=1:Nsites if norm(R0[:,i].-O)<=R]

    return Dict( (origin,radius) => pick(origin,radius) for origin ∈ origins ) |> delete_empty
end


function partition(
    indices::Vector{Int},
    R0::Array{T,2},
    r_min::Float64;
    drmin = 1e-3
    ) where {T<:Real}
    if length(indices) > _GLOBAL_MINIMAL_PARTITION_SIZE_
        return partition( refine(R0[:,indices], indices, drmin), R0, r_min, drmin=drmin )
    else
        if length(indices)==0
            println(length(indices))
        end
        return copy(indices)
    end
end


function partition(
    R0::Array{T,2},
    r_min::Float64;
    drmin = 1e-3
    ) where {T<:Real}
    partition(collect(1:size(R0,2)), R0, r_min; drmin = drmin)
end


function partition(
    dic::Dict,
    R0::Array{T,2},
    r_min::Float64;
    drmin = 1e-3
    ) where { T<:Real }
    return Dict( k=>((k[2]<=r_min) ? copy(v) : partition(v,R0,r_min,drmin=drmin)) for (k,v) ∈ dic )
end

# ------------------------------------------------------------------------------

# find_index_of_v_in_R0

function index_R0(
    v::Vector{T1},              # the coordinate of the point to search
    R0::Array{T2,2},            # the list of points
    Part::Vector{Int64},        # list of indices of R0, specifies a subset
    EPS::Float64                # precision control
    )::Int64  where { T1<:Real, T2<:Real }
    # @info "linear search"
    # find first index of v in R0[:,Part] w.r.t. `Part`
    ind = findfirst( vec(sum(abs.(R0[:,Part].-v),dims=1)).<EPS )
    # return the ACTUAL index  Part[ind]
    return ind===nothing ? -1 : Part[ind]
end


# the brute-force version
function index_R0(
    v::Vector{T1},              # the coordinate of the point to search
    R0::Array{T2,2},            # the list of points
    EPS::Float64                # precision control
    )::Int64  where { T1<:Real, T2<:Real }
    index_R0( v, R0, collect(1:size(R0,2)), EPS )
end


# O(log(n)) version
#+ SLOW !!! #+ optimize
function index_R0(
    v::Vector{T1},              # the coordinate of the point to search
    R0::Array{T2,2},            # the list of points
    Part::Dict,                 # the result of partition() and refine()
                                # Dict( (orig,rad) => pick(orig,rad) for orig ∈ origins )
    EPS::Float64                # precision control
    )::Int64  where { T1<:Real, T2<:Real }
    for ((origin,radius),part) ∈ Part
        if norm(v.-origin)<=radius
            # @info "new part"    
            return index_R0(v, R0, part, EPS)
        end
    end
    # for loop has been completed ==> not found
    return -1
end


##

#using Profile
#Profile.init(n = 10^9, delay = 0.0001)
#Profile.clear()

##

Coords = rand(3,25000) ;
p = partition(Coords, 0.02, drmin=0.001) ;

##

randk = rand(1:25000, 250000) ;
@time  A = [index_R0(Coords[:,k],Coords,p,1e-9) for k in randk] ;

##
