mutable struct HoppingHamiltonianenKp
    K::Vector
    LATT::Lattice
    KMATS::Vector                    # (kx,) => Φ ; (kx,ky) => Φ ; ...
    TS::Vector{HoppingParameter}     # spring constants Dict(:f1=>1.0, :f2=>1.2)
    COEFFS::Vector
    VS::Vector
end

# ------------------------------------------------------

function HoppingHamiltonian(
    q::Vector{Float64},
    HQSKP::HoppingHamiltonianenKp
    )
    return sum( HQSKP.COEFFS[i] .* sum(Δk(q,HQSKP.K,dk) .* KMAT[dk] for dk ∈ keys(KMAT))
                for (i,KMAT) ∈ enumerate(HQSKP.KMATS) )
end


function Kp(
    KVEC::Vector{T},
    HQS::HoppingHamiltonianenQ;
    order=2,
    eps=1e-8
    )::HoppingHamiltonianenKp where { T<:Number }
    dim = dimensions(HQS.LATT)
    @assert dim==3 || dim==2
    KL = dim==3 ? MakeKL3d(order) : MakeKL2d(order)
    @inline ce(M) = collect_expansion(M,eps)
    a = HQS.LATT.UC.a
    MATS = [ ce(Dict(s => ∇(MAT, KVEC, klist, a) for (s,klist) ∈ KL))
             for MAT ∈ HQS.MATS ]
    return HoppingHamiltonianenKp( KVEC, copy(HQS.LATT),
                                   MATS, copy.(HQS.TS), copy(HQS.COEFFS),
                                   copy.(HQS.VS) )
end
