mutable struct MillerIndex
    dim::Int64
    hkl::Vector{Int64}
end

@inline is_MillerIndex(M) = ( Set(fieldnames(MillerIndex))==Set(fieldnames(typeof(M))) )

function to_projector( M::MillerIndex, a )
    @assert M.dim==3
    @assert size(a)==(3,3)
    V = normalize( a * M.hkl )
    if abs(V[1])<1e-10 && abs(V[2])<1e-10
        return Float64[[1,0] [0,1] [0,0]]
    else
        zhat = Float64[0,0,1]
        p1 = normalize(cross3(V,zhat))
        p2 = normalize(cross3(V,p1))
        return vcat(p1',p2')
    end
end

# A = [[1,0,0] [0,1,0] [0,0,1]]
# to_projector(MillerIndex(3,[1,1,1]),A)
