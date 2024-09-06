function CommonProjectors_AO_LT_ATOM(UC::LatticeLab.UnitCell, MASS::Dict)
    OPS= vcat([[(k,v)->AcousticWeight(v,UC,MASS)*abs(v'*(LongitudinalProjector(k,UC)*AtomProjector(a,UC))*v), 
                (k,v)->AcousticWeight(v,UC,MASS)*abs(v'*(TransverseProjector(k,UC)*AtomProjector(a,UC))*v), 
                (k,v)->OpticalWeight( v,UC,MASS)*abs(v'*(LongitudinalProjector(k,UC)*AtomProjector(a,UC))*v), 
                (k,v)->OpticalWeight( v,UC,MASS)*abs(v'*(TransverseProjector(k,UC)*AtomProjector(a,UC))*v),]
                for a ∈ unique(UC.m)]...)
    return OPS
end


function CommonProjectors_LT_ATOM(UC::LatticeLab.UnitCell)
    OPS= vcat([[(k,v)->abs(v'*(LongitudinalProjector(k,UC)*AtomProjector(a,UC))*v), 
                (k,v)->abs(v'*(TransverseProjector(k,UC)*AtomProjector(a,UC))*v), ]
                for a ∈ unique(UC.m)]...)
    return OPS
end


function CommonProjectors_CH_ATOM(UC::LatticeLab.UnitCell, kref)
    OPS= vcat([[(k,v)->max(+real(v'*(ChiralityProjector(kref,UC)*AtomProjector(a,UC))*v),0),
                (k,v)->max(-real(v'*(ChiralityProjector(kref,UC)*AtomProjector(a,UC))*v),0),] 
                for a ∈ unique(UC.m)]...)
    return OPS
end


function CommonProjectors_CH(UC::LatticeLab.UnitCell, kref)
    OPS= [(k,v)->max(+real(v'*(ChiralityProjector(kref,UC))*v),0),
          (k,v)->max(-real(v'*(ChiralityProjector(kref,UC))*v),0),]
    return OPS
end


@inline ChiralityProjector2D_T(UC::LatticeLab.UnitCell,T) = SparseArrays.kron(LatticeLab.speye(UC.nsubl),T)
CommonProjectors_σ_2D(UC::LatticeLab.UnitCell,σ) = [(k,v)->real(v'*ChiralityProjector2D_T(UC,σ)*v),]
CommonProjectors_σx_2D(UC::LatticeLab.UnitCell) = CommonProjectors_σ_2D(UC, [0   1 ; 1  0])
CommonProjectors_σy_2D(UC::LatticeLab.UnitCell) = CommonProjectors_σ_2D(UC, [0 -im ; im 0])
CommonProjectors_σz_2D(UC::LatticeLab.UnitCell) = CommonProjectors_σ_2D(UC, [1   0 ; 0 -1])

CommonProjectors_σ_ATOM_2D(UC::LatticeLab.UnitCell,σ) = vcat([ 
    [(k,v)->(v'*(ChiralityProjector2D_T(UC,σ)*AtomProjector(a,UC))*v),] 
    for a ∈ unique(UC.m)
]...)
CommonProjectors_σx_ATOM_2D(UC::LatticeLab.UnitCell) = CommonProjectors_σ_ATOM_2D(UC, [0   1 ; 1  0])
CommonProjectors_σy_ATOM_2D(UC::LatticeLab.UnitCell) = CommonProjectors_σ_ATOM_2D(UC, [0 -im ; im 0])
CommonProjectors_σz_ATOM_2D(UC::LatticeLab.UnitCell) = CommonProjectors_σ_ATOM_2D(UC, [1   0 ; 0 -1])

CommonProjectors_σ_ATOM_2D_PM(UC::LatticeLab.UnitCell,σ) = vcat([
    [(k,v)->max(+real(v'*(ChiralityProjector2D_T(UC,σ)*AtomProjector(a,UC))*v),0),
     (k,v)->max(-real(v'*(ChiralityProjector2D_T(UC,σ)*AtomProjector(a,UC))*v),0),] 
    for a ∈ unique(UC.m)
]...)
CommonProjectors_σx_ATOM_2D_PM(UC::LatticeLab.UnitCell) = CommonProjectors_σ_ATOM_2D_PM(UC, [0   1 ; 1  0])
CommonProjectors_σy_ATOM_2D_PM(UC::LatticeLab.UnitCell) = CommonProjectors_σ_ATOM_2D_PM(UC, [0 -im ; im 0])
CommonProjectors_σz_ATOM_2D_PM(UC::LatticeLab.UnitCell) = CommonProjectors_σ_ATOM_2D_PM(UC, [1   0 ; 0 -1])



function CommonProjectors_ATOM(UC::LatticeLab.UnitCell, weight::Dict{Symbol,Float64}=Dict{Symbol,Float64}())
    OPS= [(k,v)->abs(v'*(get(weight,a,1.0).*AtomProjector(a,UC)*v)) for a ∈ unique(UC.m)]
    return OPS
end

CommonProjectors_NX_4x4(UC::LatticeLab.UnitCell) = [(k,v)->real(v'*(AtomProjector(Dict((UC.m[1],UC.m[2])=>1.0,(UC.m[1],UC.m[2])=>1.0),UC)*v)),]
CommonProjectors_NY_4x4(UC::LatticeLab.UnitCell) = [(k,v)->real(v'*(AtomProjector(Dict((UC.m[1],UC.m[2])=>-1.0im,(UC.m[2],UC.m[1])=>1.0im),UC)*v)),]
CommonProjectors_NZ_4x4(UC::LatticeLab.UnitCell) = [(k,v)->real(v'*(AtomProjector(Dict((UC.m[1],UC.m[1])=>1.0,(UC.m[2],UC.m[2])=>-1.0),UC)*v)),]


function CommonProjectors_ATOM_PM(UC::LatticeLab.UnitCell, AtmMat::Dict{Tuple{Symbol,Symbol},T}) where {T<:Number}
    OPS= [ (k,v)-> max( 0.0, real(v'*(AtomProjector(AtmMat,UC)*v))),
           (k,v)-> min(-0.0,-real(v'*(AtomProjector(AtmMat,UC)*v))) ] 
    return OPS
end


arg(x) = abs(x)<1e-10 ? 0 : mod(angle(x)+4π,2π)


function CommonProjectors_ATOM_entanglement_4x4(UC::LatticeLab.UnitCell)
    @assert UC.nsubl==2 && UC.dim==2
    OPS= [ (k,v)->(mod((arg(v[1])+arg(v[2])-arg(v[3])-arg(v[4]))+4π,2π)) ]
    return OPS
end
