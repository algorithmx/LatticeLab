##: ----------------------------------------------
##: projectors for DynamicalMatrixQ

## acoustic and optic

function norm2_sum_u_over_UC(v::Vector{T},UC::UnitCell,MASS::Dict) where {T<:Number}
    Minvsqrt = [1.0/sqrt(get(MASS,μ,1.0)) for μ ∈ UC.m]'
    u_sum_over_UC = vec( sum(reshape(v,(UC.dim,:)) .* Minvsqrt, dims=2) )
    return abs(dot(u_sum_over_UC,u_sum_over_UC))
end


function AcousticWeightDQ(v::Vector{T},UC::UnitCell,MASS::Dict) where {T<:Number}
    return sqrt(norm2_sum_u_over_UC(v,UC,MASS)/(num_sublattice(UC)))
end


function OpticalWeightDQ(v::Vector{T},UC::UnitCell,MASS::Dict) where {T<:Number}
    return sqrt(max(0.0, 1.0 - norm2_sum_u_over_UC(v,UC,MASS)/(num_sublattice(UC))))
end


## longitudinal and transerse

function LongitudinalProjectorDQ(kvec::Vector, UC::UnitCell)
    khat = normalize(kvec)
    BlockM = khat*khat'

    SublatticeM = speye(num_sublattice(UC))

    return kron( BlockM, SublatticeM )
end


function TransverseProjectorDQ(kvec::Vector, UC::UnitCell)
    khat = normalize(kvec)
    BlockM = eye(length(khat)) .- khat*khat'

    SublatticeM = speye(num_sublattice(UC))

    return kron( BlockM, SublatticeM )
end


function XYProjectorDQ(UC::UnitCell)
    khat = [0.0,0.0,1.0]
    BlockM = speye(3) .- khat*khat'

    SublatticeM = speye(num_sublattice(UC))

    return kron( BlockM, SublatticeM )
end


function ZProjectorDQ(UC::UnitCell)
    khat = [0.0,0.0,1.0]
    BlockM = khat*khat'

    SublatticeM = speye(num_sublattice(UC))

    return kron( BlockM, SublatticeM )
end


## subllattice-resolved

function SublatticeProjectorDQ(
    Projector::T,
    SublatticeIndexBasis::Vector{Int64},
    UC::UnitCell
    ) where { T<:AbstractMatrix }
    Nsubl = num_sublattice(UC)
    @assert all(SublatticeIndexBasis.<=Nsubl) && all(SublatticeIndexBasis.>0)
    SublatticeM = speye(UC.dim)

    BlockM = zero(ComplexF64).*speye(Nsubl)
    for (i,sli) ∈ enumerate(SublatticeIndexBasis)
        for (j,slj) ∈ enumerate(SublatticeIndexBasis)
            BlockM[sli,slj] = Projector[i,j]
        end
    end

    return kron( BlockM, SublatticeM )
end

@inline SublatticeProjectorDQ(Subl::UnitRange{Int64},UC::UnitCell) = SublatticeProjectorDQ(speye(length(Subl)), collect(Subl), UC)
@inline SublatticeProjectorDQ(Subl::Int64,UC::UnitCell) = SublatticeProjectorDQ(speye(1), [Subl,], UC)


## atom-type-resolved

function AtomProjectorDQ(Atms::Vector{Symbol},UC::UnitCell)
    AtomIndices = sort(unique( vcat([findall(UC.m.==A) for A ∈ unique(Atms) if A ∈ UC.m]...) ))
    return SublatticeProjectorDQ(speye(length(AtomIndices)), AtomIndices, UC)
end

@inline AtomProjectorDQ(A::Symbol,UC::UnitCell) = AtomProjectorDQ([A,],UC)


##: ----------------------------------------------
##: ----------------------------------------------
##: projectors for HoppingHamiltonianQ


function ProjectorSublatticeOrbitHQ(
    OrbitProjector::T1,
    SublatticeProjector::T2,
    SublatticeIndexBasis::Vector{Int64},
    UC::UnitCell
    ) where { T1<:AbstractMatrix, T2<:AbstractMatrix }
    BLK, dimH = generate_blocks_for_HQ(UC)
    P = zeros(ComplexF64,(dimH,dimH))
    Nsubl = num_sublattice(UC)
    @assert Nsubl == length(BLK)
    @assert all(SublatticeIndexBasis.>0) && all(SublatticeIndexBasis.<=Nsubl)

    for (i,sli) ∈ enumerate(SublatticeIndexBasis)
        for (j,slj) ∈ enumerate(SublatticeIndexBasis)
            α = SublatticeProjector[i,j]
            if abs(α)>1e-16
                P[BLK[sli],BLK[slj]] = α .* OrbitProjector
            end
        end
    end

    return sparse(P)
end


function ProjectorAtomHQ(
    OrbitProjector::T1,
    Atms::Vector{Symbol},
    UC::UnitCell
    ) where { T1<:AbstractMatrix }
    Subl = select_atoms_from_uc(UC,Atms)
    return ProjectorSublatticeOrbitHQ(OrbitProjector,speye(length(Subl)),Subl,UC)
end

@inline ProjectorAtomHQ(OrbitProjector,Atm::Symbol,UC::UnitCell) = ProjectorAtomHQ(OrbitProjector,[Atm,],UC)


function ProjectorAtomOrbitHQ(
    OrbitProjector::T1,
    AtomProjector::T2,
    AtomBasis::Vector{Symbol},
    UC::UnitCell
    ) where { T1<:AbstractMatrix, T2<:AbstractMatrix }
    BLK, dimH = generate_blocks_for_HQ(UC)
    P = zeros(ComplexF64,(dimH,dimH))
    Nsubl = num_sublattice(UC)
    @assert Nsubl == length(BLK)
    @assert all([a ∈ AtomBasis for a ∈ UC.m])

    for (i,atmi) ∈ enumerate(AtomBasis)
        SLi = findall(UC.m.==atmi)
        for (j,atmj) ∈ enumerate(AtomBasis)
            α = AtomProjector[i,j]
            if abs(α)>1e-16
                SLj = findall(UC.m.==atmj)
                for sli ∈ SLi
                    for slj ∈ SLj
                        P[BLK[sli],BLK[slj]] = α .* OrbitProjector
                    end
                end
            end
        end
    end

    return sparse(P)
end


## atom and orbit


function SublatticeOrbitHQ(
    Subl::Vector{Int64},
    Orbits::Vector{Symbol},
    UC::UnitCell
    )
    BLK, dimH = generate_blocks_for_HQ(UC)
    P = zeros(ComplexF64,(dimH,dimH))

    for s ∈ Subl
        P[BLK[s],BLK[s]] = diagm(0=>[((orbit ∈ Orbits) ? 1.0+0.0im : 0.0+0.0im) for orbit ∈ UC.ξ[s]])
    end

    return sparse(P)
end


@inline OrbitHQ(Orbits::Vector{Symbol}, UC::UnitCell) = SublatticeOrbitHQ(collect(1:UC.nsubl), Orbits, UC)

@inline OrbitHQ(Orbit::Symbol, UC::UnitCell) = OrbitHQ([Orbit,], UC)

@inline SublatticeHQ(Subl::Vector{Int64}, UC::UnitCell) = SublatticeOrbitHQ(Subl, all_orbits(UC), UC)

@inline SublatticeHQ(Subl::Int64, UC::UnitCell) = SublatticeHQ([Subl,], UC)

@inline AtomOrbitHQ(Atms::Vector{Symbol}, Orbits::Vector{Symbol}, UC::UnitCell) = SublatticeOrbitHQ(select_atoms_from_uc(UC,Atms), Orbits, UC)

@inline AtomHQ(Atms::Vector{Symbol}, UC::UnitCell) = AtomOrbitHQ(Atms, all_orbits(UC), UC)

@inline AtomHQ(Atm::Symbol, UC::UnitCell) = AtomHQ([Atm,], UC)
