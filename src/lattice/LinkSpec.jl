mutable struct LinkSpec{R<:Real}
    Name::Symbol
    Directions::Vector{Vector{R}}
    Distances::Vector{R}
    Sublattices::Vector{Tuple{Int,Int}}
    Nth::Int
end

##

#Spring = Dict{Tuple{Int64,Int64},Vector{Direction}}

#Dict{Symbol,Spring}



##

