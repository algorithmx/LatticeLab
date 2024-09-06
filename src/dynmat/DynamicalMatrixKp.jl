mutable struct DynamicalMatrixKp
    K::Vector
    LATT::Lattice
    KMAT::Dict               # (kx,) => Φ ; (kx,ky) => Φ ; ...
    SP::Dict        # spring constants Dict(:f1=>1.0, :f2=>1.2)
    M::Dict         # mass Dict(:mA=>1.0, :mB=>1.2)
end


mutable struct DynamicalMatrixSymbolicKp
    K::Vector
    LATT::Lattice
    MATS::Dict              # var => ( mn => Φ )
    SP::Dict                # spring constants Dict(:f1=>1.0, :f2=>1.2)
    SPparams::Dict          # parameters of functions in the above dict :f1 => (:x,:y,)
    M::Dict                 # mass Dict(:mA=>1.0, :mB=>1.2)
    Mparams::Dict           # parameters of functions in the above dict :m => (:x,)
end


# ------------------------------------------------------

# favourite function
function DynamicalMatrix(
    q::Vector{Float64},
    DQKP::DynamicalMatrixKp
    )
    # the Fourier sum
    return sum( Δk(q,DQKP.K,dk) .* DQKP.KMAT[dk] for dk ∈ keys(DQKP.KMAT) )
end


function Kp(
    KVEC::Vector{T},
    DQ::DynamicalMatrixQ;
    order=4,
    eps=1e-8
    )::DynamicalMatrixKp where { T<:Number }
    dim = dimensions(DQ.LATT)
    @assert dim==3 || dim==2
    KL = dim==3 ? (MakeKL3d(order)) : (MakeKL2d(order))
    @inline ce(M) = collect_expansion(M,eps)
    MAT = Dict( s => ∇(DQ.MAT, KVEC, klist, DQ.LATT.UC.a)
                for (s,klist) ∈ KL ) |> ce
    return DynamicalMatrixKp( KVEC, copy(DQ.LATT), MAT,
                              copy(DQ.SP), copy(DQ.M) )
end


function Kp(
    KVEC::Vector{T},
    DSQ::DynamicalMatrixSymbolicQ;
    order=4,
    eps=1e-8
    )::DynamicalMatrixSymbolicKp where { T<:Number }
    dim = dimensions(DSQ.LATT)
    @assert dim==3 || dim==2
    KL = dim==3 ? (MakeKL3d(order)) : (MakeKL2d(order))
    @inline ce(M) = collect_expansion(M,eps)
    MATS = Dict( sym => ( Dict( s => ∇(DSQ.MATS[sym], KVEC, klist, DSQ.LATT.UC.a)
                                for (s,klist) ∈ KL ) |> ce )
                 for sym ∈ keys(DSQ.MATS) )
    return DynamicalMatrixSymbolicKp( KVEC, copy(DSQ.LATT), MATS,
                                      copy(DSQ.SP), copy(DSQ.SPparams),
                                      copy(DSQ.M),  copy(DSQ.Mparams) )
end
