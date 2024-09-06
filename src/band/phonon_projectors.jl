
##* ==================== for DynamicalMatrixQ ====================

# acoustic and optic

function CM_coords_portion2(v::Vector{T}, UC::UnitCell, MASS::Dict) where {T<:Number}
    M = [get(MASS,μ,1.0) for μ ∈ UC.m]'
    # ∑ᵢ mᵢ * uᵢ
    center_of_mass_coords = vec(sum(reshape(v,(UC.dim,:)).*M, dims=2))
    return abs((center_of_mass_coords'*center_of_mass_coords)/(M*M'))
end


function AcousticWeight(v::Vector{T}, UC::UnitCell, MASS::Dict) where {T<:Number}
    return CM_coords_portion2(v,UC,MASS)
end


function OpticalWeight(v::Vector{T}, UC::UnitCell, MASS::Dict) where {T<:Number}
    return max(0.0, 1-CM_coords_portion2(v,UC,MASS))
end


# longitudinal and transerse

function LongitudinalProjector(kvec::Vector, UC::UnitCell)
    khat = norm(kvec)<1e-10 ? kvec : normalize(kvec)
    BlockM = khat*khat'
    SublatticeM = speye(UC.nsubl)
    return SparseArrays.kron(BlockM, SublatticeM)
end


function TransverseProjector(kvec::Vector, UC::UnitCell)
    khat = norm(kvec)<1e-10 ? kvec : normalize(kvec)
    BlockM = Matrix(I,length(khat),length(khat)) .- khat*khat'
    SublatticeM = speye(UC.nsubl)
    return SparseArrays.kron(BlockM, SublatticeM)
end


function XYProjector(UC::UnitCell)
    khat = [0.0,0.0,1.0]
    BlockM = speye(3) .- khat*khat'
    SublatticeM = speye(UC.nsubl)
    return SparseArrays.kron(BlockM, SublatticeM)
end


function ZProjector(UC::UnitCell)
    khat = [0.0,0.0,1.0]
    BlockM = khat*khat'
    SublatticeM = speye(UC.nsubl)
    return SparseArrays.kron(BlockM, SublatticeM)
end


# subllattice-resolved
function SublatticeProjector(
    SublatticeIndexBasis::Vector{Int64},
    Weights::Vector{Float64},
    UC::UnitCell
    )
    @assert length(Weights)==length(SublatticeIndexBasis) 
    @assert all(SublatticeIndexBasis.<=UC.nsubl) && all(SublatticeIndexBasis.>0)
    BlockM = sparse(Int[],Int[],ComplexF64[],UC.nsubl,UC.nsubl)
    for (i,sli) ∈ enumerate(SublatticeIndexBasis)
        BlockM[sli,sli] = Weights[i]  #! TODO error !!!
    end
    return SparseArrays.kron(BlockM, speye(UC.dim))
end

function SublatticeProjectorI(
    SublatticeIndexBasis::Vector{Int64},
    UC::UnitCell
    )
    @assert all(SublatticeIndexBasis.<=UC.nsubl) && all(SublatticeIndexBasis.>0)
    dd = zeros(UC.nsubl)
    dd[SublatticeIndexBasis] .= 1
    return SparseArrays.kron(spdiagm(0=>dd), speye(UC.dim))
end

# subllattice-resolved
SublatticeProjector(SublatticeIndexBasis::Vector{Int64}, UC::UnitCell) = 
    SublatticeProjectorI(SublatticeIndexBasis, UC)


@inline SublatticeProjector(Subl::UnitRange{Int64}, UC::UnitCell) = SublatticeProjector(collect(Subl), UC)
@inline SublatticeProjector(Subl::Int64,            UC::UnitCell) = SublatticeProjector([Subl,],       UC)


# atom-type-resolved
function AtomProjector(Atms::Vector{Symbol},UC::UnitCell)
    # sublattice indices for atoms in list `Atms`
    AtomIndices = sort(unique(vcat([findall(UC.m.==A) for A ∈ unique(Atms) if A ∈ UC.m]...)))
    return SublatticeProjector(AtomIndices, UC)
end

@inline AtomProjector(A::Symbol,UC::UnitCell) = AtomProjector([A,],UC)

# atom-type-resolved
function AtomProjector(AtmMat::Dict{Tuple{Symbol,Symbol},T},UC::UnitCell) where {T<:Number}
    kAtmMat = keys(AtmMat)
    BlockM = sparse(Int[],Int[],ComplexF64[],UC.nsubl,UC.nsubl)
    for (i,si) ∈ enumerate(UC.m)
        for (j,sj) ∈ enumerate(UC.m)
            if (si,sj) ∈ kAtmMat
                BlockM[i,j] = AtmMat[(si,sj)]
            end
        end
    end
    return SparseArrays.kron(BlockM, speye(UC.dim))
end



function ChiralityProjector(kvec::Vector, UC::UnitCell)
    SublatticeM = speye(UC.nsubl)
    if norm(kvec)<1e-4
        return SparseArrays.kron(0 .*(kvec*(kvec')), SublatticeM)
    end
    (a,b,c) =  normalize(kvec)
    if abs(b)>1e-5 || abs(c)>1e-5
        α = normalize(cross([a,b,c], [a,-c,b]))
        β = normalize(cross([a,b,c], α))
    else
        α = normalize(cross([a,b,c], [-c,b,a]))
        β = normalize(cross([a,b,c], α))
    end 
    BlockM = im.*(α * β' .- β * α')
    return SparseArrays.kron(BlockM, SublatticeM)
end
