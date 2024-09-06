sortperm_by_vec(A::Array{T,2}) where {T<:Real} = sortperm(collect(1:size(A,2)), by=i->A[:,i])

function index_R0(
    v0::Vector{T1},              # the coordinate of the point to search
    R0::Array{T2,2},            # the list of points
    PosIndx,
    EPS::Float64,
    prec_control::Int
    )::Int64  where { T1<:Real, T2<:Real }
    DIM = size(R0,1)
    v  = round.(v0,digits=prec_control) .- fill(0.5EPS,DIM)
    L  = 1
    R  = size(R0,2)
    M  = (L+R) ÷ 2
    PM = PosIndx[M]
    vm = zeros(T1,DIM)
    while L <= R
        @inbounds  vm = round.(R0[:,PM],digits=prec_control)
        if norm(v.-vm) < EPS
            return PM
        elseif v > vm  # larger for sure
            L = M+1
            M = (L+R) ÷ 2
            @inbounds PM = PosIndx[M]
        else # v < vm
            R = M-1
            M = (L+R) ÷ 2
            @inbounds PM = PosIndx[M]
        end
    end
    return -1
end

##
#=

using LinearAlgebra

#using Profile
#Profile.init(n = 10^9, delay = 0.0001)
#Profile.clear()

Coords = 10rand(3,25000) ;
PV = sortperm_by_vec(Coords) ;
randk = rand(1:25000, 250000) ;

# 0.398138 seconds ~ 0.665311 seconds
@time res = [ index_R0(Coords[:,k], Coords, PV, 1e-8,9) for k in randk ] ;

res==randk

=#