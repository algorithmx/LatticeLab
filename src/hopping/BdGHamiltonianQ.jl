#TODO trim and merge elements in BdGHamiltonianenQ

#TODO how to record the action of iσy on sublattice orbits? (it permutes the orbits...)

# BdG Hamiltonian is defined with respect to
# the Nambu basis
# Ψ = [ψ  T*ψ] = [ψ  iσy*K(ψ)]
# The diagonal part of the extended Hamiltonian is
# H₁₁ = H
# H₂₂ = T * H * T⁻¹ = (iσy*K) * H * (iσy*K)⁻¹ = - iσy * (H*) * iσy⁻¹
# iσy stands for the unitary transformation that rearranges
# the components of ψ (due to time-reversal rules), for example,
# ( ↑1k,↓1k,↑1-k,↓1-k, ↑2k,↓2k,↑2-k,↓2-k, ⋯ )

# ------------------------------------------------------------

mutable struct BdGHamiltonianQ{T1<:AbstractMatrix}
    H::HoppingHamiltonianQ
    Δ::SuperconductingGapFunctionQ
    iσy::T1
end


mutable struct BdGHamiltonianenQ{T1<:AbstractMatrix}
    HS::HoppingHamiltonianenQ
    ΔS::SuperconductingGapFunctionenQ
    iσy::T1
end


function iden_transform(hhq::HoppingHamiltonianQ)
    return speye(size(first(values(hhq.MAT)),1))
end

function iden_transform(hhq::HoppingHamiltonianenQ)
    return speye(size(first(values(hhq.MATS[1])),1))
end


function zero_gap(hhq::HoppingHamiltonianQ)
    ZK = zeros(Int64,  size(first(keys(hhq.MAT))))
    ZM = zeros(Float64,size(first(values(hhq.MAT))))
    return SuperconductingGapFunctionQ(hhq.LATT,Dict(ZK=>ZM))
end

function zero_gap(hhq::HoppingHamiltonianenQ)
    ZK = zeros(Int64,  size(first(keys(hhq.MATS[1]))))
    ZM = zeros(Float64,size(first(values(hhq.MATS[1]))))
    return SuperconductingGapFunctionenQ(hhq.LATT, [Dict(ZK=>ZM),], [zero(ComplexF64),])
end


@inline is_BdGHamiltonianQ(bdg) = ( Set(fieldnames(BdGHamiltonianQ))==Set(fieldnames(typeof(bdg)))
                                    && is_SuperconductingGapFunction(bdg.Δ) && is_HoppingHamiltonianQ(bdg.H) )

@inline is_BdGHamiltonianenQ(bdg) = ( Set(fieldnames(BdGHamiltonianenQ))==Set(fieldnames(typeof(bdg)))
                                      && is_SuperconductingGapFunctionen(bdg.ΔS) && is_HoppingHamiltonianenQ(bdg.HS) )


check_compat(bdg::BdGHamiltonianQ) = ( (bdg.H.LATT==bdg.Δ.LATT) && check_compat(bdg.H) && check_compat(bdg.Δ) )

check_compat(bdg::BdGHamiltonianenQ) = ( (bdg.HS.LATT==bdg.ΔS.LATT) && check_compat(bdg.HS) && check_compat(bdg.ΔS) )


@inline consistency_check(bdg::BdGHamiltonianQ; TOL=1e-8) = ( consistency_check(bdg.H,TOL=TOL) && consistency_check(bdg.Δ,TOL=TOL) )

@inline consistency_check(bdg::BdGHamiltonianenQ; TOL=1e-8) = ( consistency_check(bdg.HS,TOL=TOL) && consistency_check(bdg.ΔS,TOL=TOL) )


copy(BDGQ::BdGHamiltonianQ) = BdGHamiltonianQ{typeof(BDGQ.iσy)}( copy(BDGQ.H), copy(BDGQ.Δ), copy(BDGQ.iσy) )

copy(BDGQ::BdGHamiltonianenQ) = BdGHamiltonianenQ{typeof(BDGQ.iσy)}( copy(BDGQ.HS), copy(BDGQ.ΔS), copy(BDGQ.iσy) )


@inline convert(::Type{T}, x::T) where {T<:BdGHamiltonianQ} = x
@inline convert(::Type{T}, x::T) where {T<:BdGHamiltonianenQ} = x


function convert(::Type{T}, dq) where {T<:BdGHamiltonianQ}
    if is_HoppingHamiltonianQ(dq)
        iσ2 = iden_transform(dq)
        Tiσ2 = typeof(iσ2)
        return BdGHamiltonianQ{Tiσ2}( copy(dq), zero_gap(dq), iσ2 )
    elseif is_BdGHamiltonianQ(dq)
        dqΛ = sparse( dq.iσy )
        TDQΛ = typeof(dqΛ)
        return BdGHamiltonianQ{TDQΛ}( convert( HoppingHamiltonianQ, dq.H ),
                                      convert( SuperconductingGapFunctionQ, dq.Δ ),
                                      dqΛ )
    else
        throw("Type error of dq in convert(::Type{T}, dq) where {T<:BdGHamiltonianQ}")
    end
end

function convert(::Type{T}, dq) where {T<:BdGHamiltonianenQ}
    if is_HoppingHamiltonianQ(dq)
        iσ2 = iden_transform(dq)
        Tiσ2 = typeof(iσ2)
        dq1 = convert(HoppingHamiltonianenQ,dq)
        return BdGHamiltonianenQ{Tiσ2}( dq1, zero_gap(dq1), iσ2 )
    elseif is_HoppingHamiltonianenQ(dq)
        iσ2 = iden_transform(dq)
        Tiσ2 = typeof(iσ2)
        return BdGHamiltonianenQ{Tiσ2}( copy(dq), zero_gap(dq), iσ2 )
    elseif is_BdGHamiltonianQ(dq)
        iσ2 = sparse( dq.iσy )
        Tiσ2 = typeof(dqΛ)
        return BdGHamiltonianenQ{Tiσ2}( convert( HoppingHamiltonianenQ,         dq.H ),
                                        convert( SuperconductingGapFunctionenQ, dq.Δ ),
                                        dqΛ )
    elseif is_BdGHamiltonianenQ(dq)
        iσ2 = sparse( dq.iσy )
        Tiσ2 = typeof(dqΛ)
        return BdGHamiltonianenQ{Tiσ2}( convert( HoppingHamiltonianenQ,         dq.HS ),
                                        convert( SuperconductingGapFunctionenQ, dq.ΔS ),
                                        dqΛ )
    else
        throw("Type error of dq in convert(::Type{T}, dq) where {T<:BdGHamiltonianenQ}")
    end
end


function kspace_BdG_hamiltonian(
    hhq::HoppingHamiltonianQ,
    iσy::T1
    ) where {T1<:AbstractMatrix}
    bdg = convert(BdGHamiltonianQ, hhq)
    bdg.iσy = sparse(iσy)
    return bdg
end

function kspace_BdG_hamiltonian(
    hhqs::HoppingHamiltonianenQ,
    iσy::T1
    ) where {T1<:AbstractMatrix}
    bdg = convert(BdGHamiltonianenQ, hhqs)
    bdg.iσy = sparse(iσy)
    return bdg
end

function kspace_BdG_hamiltonian(
    dq::SuperconductingGapFunctionQ,
    iσy::T1
    ) where {T1<:AbstractMatrix}
    H0 = zero_HoppingHamiltonianQ(dq.LATT)
    SPiσy = sparse(iσy)
    return BdGHamiltonianQ{typeof(SPiσy)}(H0,copy(dq),SPiσy)
end

function kspace_BdG_hamiltonian(
    dqs::SuperconductingGapFunctionenQ,
    iσy::T1
    ) where {T1<:AbstractMatrix}
    H0 = zero_HoppingHamiltonianQ(dq.LATT)
    SPiσy = sparse(iσy)
    return BdGHamiltonianenQ{typeof(SPiσy)}(convert(HoppingHamiltonianenQ,H0),copy(dqs),SPiσy)
end

# --------------------------------------------------------------------------

function +(BdG::BdGHamiltonianQ, DQ::Union{SuperconductingGapFunctionQ,SuperconductingGapFunctionenQ})
    @assert BdG.H.LATT == DQ.LATT
    return BdGHamiltonianenQ(convert(HoppingHamiltonianenQ,BdG.H), BdG.Δ+DQ, BdG.iσy)
end

function +(DQ::Union{SuperconductingGapFunctionQ,SuperconductingGapFunctionenQ}, BdG::BdGHamiltonianQ)
    return +(BdG,DQ)
end

function +(BdG::BdGHamiltonianenQ, DQ::Union{SuperconductingGapFunctionQ,SuperconductingGapFunctionenQ})
    @assert BdG.HS.LATT == DQ.LATT
    return BdGHamiltonianenQ(BdG.HS, BdG.ΔS+DQ, BdG.iσy)
end

function +(DQ::Union{SuperconductingGapFunctionQ,SuperconductingGapFunctionenQ}, BdG::BdGHamiltonianenQ)
    return +(BdG,DQ)
end

function +(BdG1::BdGHamiltonianQ, BdG2::BdGHamiltonianQ)
    @assert BdG1.H.LATT == BdG2.H.LATT
    @assert BdG1.iσy == BdG2.iσy
    return BdGHamiltonianenQ(BdG1.H+BdG2.H, BdG1.Δ+BdG2.Δ, BdG1.iσy)
end

function +(BdG1::BdGHamiltonianQ, BdGS::BdGHamiltonianenQ)
    @assert BdG1.H.LATT == BdGS.HS.LATT
    @assert BdG1.iσy == BdGS.iσy
    return BdGHamiltonianenQ(BdG1.H+BdGS.HS, BdG1.Δ+BdGS.ΔS, BdG1.iσy)
end

function +(BdGS::BdGHamiltonianenQ,BdG1::BdGHamiltonianQ)
    return +(BdG1,BdGS)
end

function +(BdGS1::BdGHamiltonianenQ,BdGS2::BdGHamiltonianenQ)
    @assert BdGS1.HS.LATT == BdGS2.HS.LATT
    @assert BdGS1.iσy == BdGS2.iσy
    return BdGHamiltonianenQ(BdGS1.HS+BdGS2.HS, BdGS1.ΔS+BdGS2.ΔS, BdGS1.iσy)
end

function *(c::N,BdG::BdGHamiltonianQ) where {N<:Number}
    BdG1 = convert(BdGHamiltonianenQ,BdG)
    BdG1.HS.COEFFS[:] .= (c.*BdG1.HS.COEFFS[:])
    BdG1.ΔS.COEFFS[:] .= (c.*BdG1.ΔS.COEFFS[:])
    return BdG1
end

function *(c::N,BdG::BdGHamiltonianenQ) where {N<:Number}
    BdG1 = copy(BdG)
    BdG1.HS.COEFFS[:] .= (c.*BdG1.HS.COEFFS[:])
    BdG1.ΔS.COEFFS[:] .= (c.*BdG1.ΔS.COEFFS[:])
    return BdG1
end

# --------------------------------------------------------------------------

@inline blkσ(M,σ) = kron(σ,M)
σ11 = sparse([[1,0] [0,0]])
@inline blkσ11(M) = blkσ(M,σ11)
σ12 = sparse([[0,0] [1,0]])
@inline blkσ12(M) = blkσ(M,σ12)
σ21 = sparse([[0,1] [0,0]])
@inline blkσ21(M) = blkσ(M,σ21)
σ22 = sparse([[0,0] [0,1]])
@inline blkσ22(M) = blkσ(M,σ22)

# favourite function
@inline K(M) = conj.(M)


# HoppingHamiltonianQ
function BdGHamiltonian(
    q::Vector,
    HMAT::Dict,
    ΔMAT::Dict,
    iσy, a ) where {T<:Dict, N}
    H = qMAT(q,HMAT,a)
    Δ = qMAT(q,ΔMAT,a)
    return (   blkσ11(H)                   # the band Hamiltonian(en)Q
            .+ blkσ12(Δ)                   #
            .+ blkσ21(Δ')                  #
            .+ blkσ22(iσy*(-K(H))*(iσy'))  #
            )
end

# HoppingHamiltonianenQ
function BdGHamiltonian(
    q::Vector,
    HC::Tuple{Vector{T1},Vector{N1}},
    ΔC::Tuple{Vector{T2},Vector{N2}},
    iσy, 
    a 
    ) where {T1<:Dict,N1<:Number,T2<:Dict,N2<:Number}
    H = qMATS(q, a, HC...)
    Δ = qMATS(q, a, ΔC...)
    return (   blkσ11(H)                   # the band Hamiltonian(en)Q
            .+ blkσ12(Δ)                   #
            .+ blkσ21(Δ')                  #
            .+ blkσ22(iσy*(-K(H))*(iσy'))  #TODO ??? -H(-k)
            )
end

function BdGHamiltonian(
    q::Vector,
    BdGQ::BdGHamiltonianQ{T1}
    ) where {T1<:AbstractMatrix}
    return BdGHamiltonian( q, BdGQ.H.MAT, BdGQ.Δ.MAT, BdGQ.iσy, BdGQ.H.LATT.UC.a )
end

function BdGHamiltonian(
    q::Vector,
    BdGQ::BdGHamiltonianenQ{T1}
    ) where {T1<:AbstractMatrix}
    return BdGHamiltonian( q,
                           (BdGQ.HS.MATS,BdGQ.HS.COEFFS),
                           (BdGQ.ΔS.MATS,BdGQ.ΔS.COEFFS),
                           BdGQ.iσy, BdGQ.HS.LATT.UC.a )
end

# -----------------------------------------------------------
