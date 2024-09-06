mutable struct HoppingHamiltonianKp
    K::Vector
    LATT::Lattice
    KMAT::Dict               # (kx,) => Φ ; (kx,ky) => Φ ; ...
    T::HoppingParameter     # spring constants Dict(:f1=>1.0, :f2=>1.2)
    V::Dict                 # on-site potential
end


mutable struct HoppingHamiltonianSymbolicKp
    K::Vector
    LATT::Lattice
    MATS::Dict              # var => ( (kx,) => Φ ; (kx,ky) => Φ ; ... )
    T::HoppingParameter     # spring constants Dict(:f1=>f1(x1,x2), :f2=>f2(x3,x4))
    Tparams::Dict           # parameters of functions in the above dict :f1 => (:x,:y,)
    V::Dict                 # on-site potential: :m => f(x)
    Vparams::Dict           # parameters of functions in the above dict :m => (:x,)
end

# ------------------------------------------------------

function HoppingHamiltonian(
    q::Vector{Float64},
    HQKP::HoppingHamiltonianKp
    )
    return sum( Δk(q,HQKP.K,dk) .* HQKP.KMAT[dk] for dk ∈ keys(HQKP.KMAT) )
end


function Kp(
    KVEC::Vector{T},
    HQ::HoppingHamiltonianQ;
    order=2,
    eps=1e-8
    )::HoppingHamiltonianKp where { T<:Number }
    dim = dimensions(HQ.LATT)
    @assert dim==3 || dim==2
    KL = dim==3 ? MakeKL3d(order) : MakeKL2d(order)
    @inline ce(M) = collect_expansion(M,eps)
    MAT = Dict( s => ∇(HQ.MAT, KVEC, klist, HQ.LATT.UC.a)
                for (s,klist) ∈ KL ) |> ce
    return HoppingHamiltonianKp( KVEC, copy(HQ.LATT), MAT,
                                 copy(HQ.T), copy(HQ.V) )
end


function Kp(
    KVEC::Vector{T},
    HSQ::HoppingHamiltonianSymbolicQ;
    order=2,
    eps=1e-8
    )::HoppingHamiltonianSymbolicKp where { T<:Number }
    dim = dimensions(HSQ.LATT)
    @assert dim==3 || dim==2
    KL = dim==3 ? MakeKL3d(order) : MakeKL2d(order)
    @inline ce(M) = collect_expansion(M,eps)
    MATS = Dict( sym => ( Dict( s => ∇(HSQ.MATS[sym], KVEC, klist, HSQ.LATT.UC.a)
                                for (s,klist) ∈ KL ) |> ce )
                 for sym ∈ keys(HSQ.MATS) )
    return HoppingHamiltonianSymbolicKp( KVEC, copy(HSQ.LATT), MATS,
                                         copy(HSQ.T), copy(HSQ.Tparams),
                                         copy(HSQ.V), copy(HSQ.Vparams) )
end
