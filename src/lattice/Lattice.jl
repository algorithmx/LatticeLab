"""
    Lattice{TConnMat<:AbstractMatrix,TR<:Real}

`Lattice` is the basic data structure of `LatticeLab.jl`. It contains 
the geometry and connectivity information of the lattice.

Fields                   # Comments
--------------------     # -------------------------------------------------
R0::Coordinates{TR}      # equilibrium positions of lattice sites
BBOX::BoundingBox        # bounding box of the piece
LN::LinkInfo             # link between pairs of lattice sites [NOTE 1]
SL::Vector{Int64}        # list of sublattice number for each site in R0
f::TConnMat              # site connectivity matrix, element is the link symbol
UC::LatticeInfo          # unit cell data or random property
EqV::Index               # index of the equivalent vertex [NOTE 2]
TransPerm::Indices       # TransPerm[:,k] ==> perm. for trans. UC.a[ : , k ]
TransNegPerm::Indices    # TransNegPerm[:,k] ==> perm. for trans. -UC.a[ : , k ]

[NOTE 1]

[NOTE 2] meaning of EqV:
EqV[i] > 0 : site i at outer boundary, equivalent site EqV[i]
EqV[i] = 0 : site i inside the bounding box
EqV[i] < 0 : site i not used

"""
mutable struct Lattice{TConnMat<:AbstractMatrix,TR<:Real}
    R0::Coordinates{TR}                     # equilibrium positions of lattice sites
    #R0P::Vector{Int}                        # sortperm,  to accelerate the site index search
    BBOX::BoundingBox                       # bounding box of the piece
    LN::LinkInfo                            # link between pairs of lattice sites [NOTE 1]
    SL::Vector{Int64}                       # list of sublattice number for each site in R0
    f::TConnMat                             # site connectivity matrix, element is the link symbol
    UC::LatticeInfo                         # unit cell data or random property
    EqV::Index                              # index of the equivalent vertex
    TransPerm::Indices                      # TransPerm[:,k] is the index [NOTE 2]
    TransNegPerm::Indices                   # TransNegPerm[:,k] is the index [NOTE 3]
end



function isequal_upto_orbits(L1::Lattice, L2::Lattice)
    if ( typeof(L1.f)==typeof(L2.f) && eltype(L1.R0)==eltype(L2.R0) && size(L1.R0)==size(L2.R0) )
        if ( (L1.BBOX == L2.BBOX) && (L1.SL == L2.SL) && (L1.EqV == L2.EqV)
            && (L1.LN == L2.LN) && isequal_upto_orbits(L1.UC,L2.UC) )
            return (L1.f==L2.f) && norm(L1.R0.-L2.R0)<1e-10
        else
            return false
        end
    else
        return false
    end
end

function isequal(L1::Lattice, L2::Lattice)
    return (isequal_upto_orbits(L1,L2) && L1.UC.ξ==L2.UC.ξ)
end

function ==(L1::Lattice, L2::Lattice)
    return isequal(L1, L2)
end


@inline is_Lattice(latt) = (   Set(fieldnames(Lattice))==Set(fieldnames(typeof(latt)))
                            && is_Coordinates(latt.R0)
                            #&& isa(latt.R0P,  Vector{Int})
                            && is_BoundingBox(latt.BBOX)
                            && is_LinkInfo(latt.LN)
                            && is_Index(latt.SL)
                            && is_Index(latt.EqV)
                            && is_Indices(latt.TransPerm)
                            && is_Indices(latt.TransNegPerm)  )

copy(latt::Lattice) = Lattice( copy( latt.R0            ),
                               #copy( latt.R0P           ),
                               copy( latt.BBOX          ),
                               copy( latt.LN            ),
                               copy( latt.SL            ),
                               copy( latt.f             ),
                               copy( latt.UC            ),
                               copy( latt.EqV           ),
                               copy( latt.TransPerm     ),
                               copy( latt.TransNegPerm  ) )


@inline lattice_dimensions(latt::Lattice) = size(latt.R0)

@inline dimensions(latt::Lattice) = lattice_dimensions(latt)[1]

@inline num_sites(latt::Lattice) = lattice_dimensions(latt)[2]

@inline num_sublattice(latt::Lattice) = latt.UC.nsubl

@inline num_inner_sites(latt::Lattice) = length(findz(latt.EqV))

@inline inner_site_id_from_eqv(eqv) = findz(eqv)

@inline inner_site_id(latt::Lattice) = inner_site_id_from_eqv(latt.EqV)

@inline outer_border_site_id(latt::Lattice) = findall(latt.EqV.>0)

@inline inner_equivalent_of_outer_border_site_id(latt::Lattice) = latt.EqV[outer_border_site_id(latt)]

@inline is_inner_site(latt::Lattice, site_i::Int) = ( latt.EqV[site_i]==0  )

@inline function outer_border_and_inner_site_id(latt::Lattice)
    I,J,F = findnz(latt.f)
    IJ = collect(zip(I,J))
    return sort(unique(vcat([ [i,j]
                              for (i,j) ∈ IJ
                              if (is_inner_site(latt,i) || is_inner_site(latt,j)) ]...)))
end

@inline inner_site_coords(latt::Lattice) = latt.R0[:,inner_site_id(latt)]

@inline inner_site_sublattices(latt::Lattice) = latt.SL[inner_site_id(latt)]

@inline inner_site_atom_labels(latt::Lattice,atom_labels::Vector) = map(x->findfirst(y->y==x,atom_labels), atom_labels[inner_site_sublattices(latt)])

@inline inner_site_atom_labels(latt::Lattice) = inner_site_atom_labels(latt, latt.UC.m)

@inline inner_site_orbit_labels(latt::Lattice) = latt.UC.ξ[inner_site_sublattices(latt)]

@inline orbit_numbers(latt::Lattice) = orbit_numbers(latt.UC)

@inline orbit_number(latt::Lattice, subl_i::Int) = orbit_number(latt.UC,subl_i)

@inline total_num_orbits(latt::Lattice) = total_num_orbits(latt.UC)

@inline mass_site(latt::Lattice, site_i::Int) = mass_sublattice(latt.UC, latt.SL[site_i])

@inline is_outer_site_at_boundary(latt::Lattice, site_i::Int) = ( latt.EqV[site_i]>0   )

@inline is_outer_site_irrelavent(latt::Lattice, site_i::Int)  = ( latt.EqV[site_i]==-1 )

@inline equiv_site(latt::Lattice, site_i::Int) = (is_inner_site(latt,site_i) ? site_i : latt.EqV[site_i])

@inline reciprocal_basis(latt::Lattice) = reciprocal_basis(latt.UC)

@inline relative_reciprocal_vector(b_abs::Vector,latt::Lattice) = inv(reciprocal_basis(latt.UC))*b_abs

@inline function num_unitcells(latt::Lattice)
    @assert ( num_inner_sites(latt) % num_sublattice(latt) ) == 0
    return num_inner_sites(latt) ÷ num_sublattice(latt)
end

@inline origin_shift(latt::Lattice) = origin_shift(latt.BBOX,latt.UC.a)

@inline unique_link_types(latt::Lattice) = sort(unique(last(findnz(latt.f))))

check_compat(latt::Lattice) = ((isa(latt.UC,UnitCell) ? latt.UC.dim==dimensions(latt) : true)
                            # && (isa(latt.UC,RandomProperty) ? latt.UC.Nsites==num_sites(latt) : true)
                            && size(latt.f,1)==size(latt.f,2)==num_sites(latt)
                            && check_valid_EqV(latt.EqV)
                            && length(latt.EqV)==num_sites(latt)
                            && size(latt.TransPerm,1)==num_sites(latt)
                            && size(latt.TransPerm,2)==dimensions(latt)
                            && size(latt.TransNegPerm,1)==num_sites(latt)
                            && size(latt.TransNegPerm,2)==dimensions(latt)
                            && length(latt.SL)==num_sites(latt)
                            && check_compat(latt.UC)
                            && ( num_inner_sites(latt) % latt.UC.nsubl ) == 0
                            && check_compat(latt.BBOX)
                            && check_compat(latt.LN)  )


typeparam1(::Type{Lattice{A,B}}) where {A, B} = A
typeparam1(::Type{T}) where {T} = nothing
typeparam2(::Type{Lattice{A,B}}) where {A, B} = B
typeparam2(::Type{T}) where {T} = nothing


@inline convert(::Type{T}, latt::T) where {T<:Lattice} = latt

function convert(::Type{T}, latt) where {T<:Lattice}
    @assert is_Lattice(latt)
    T1 = typeof(latt.f)
    @assert T1<:AbstractMatrix
    # lazy, no type-promotion of elements
    A = (typeparam1(T)===nothing ? T1 : typeparam1(T))

    TR = eltype(latt.R0)
    @assert TR<:Real
    # lazy, no type-promotion of elements
    B = (typeparam2(T)===nothing ? TR : typeparam2(T))

    #XXX other kinds of LatticeInfo ?
    TLI = is_UnitCell(latt.UC) ? UnitCell : LatticeInfo
    LI = latt.UC
    if TLI==UnitCell
        LI = UnitCell( Int64(latt.UC.dim),
                       Int64(latt.UC.nsubl),
                       Coordinates{B}(latt.UC.a),
                       Coordinates{B}(latt.UC.δ),
                       Masses(latt.UC.m),
                       Vector{Orbits}(latt.UC.ξ) )
    end

    Lattice{A,B}( Coordinates{B}( latt.R0              ),
                  #Vector{Int}(    latt.R0P             ),
                  convert( BoundingBox,      latt.BBOX ),
                  convert( LinkInfo,    latt.LN   ),
                  Index(          latt.SL              ),
                  convert( A,                latt.f    ),
                  convert( TLI,                   LI   ),
                  Index(          latt.EqV             ),
                  Indices(        latt.TransPerm       ),
                  Indices(        latt.TransNegPerm    )   )
end


comp_vec(v1,v2,EPS) = (abs(v1[1]-v2[1])<EPS && (abs(v1[2]-v2[2])<EPS && (length(v1)==3 ? abs(v1[3]-v2[3])<EPS : true)))


function find_index_nocheck(
    v,
    latt::Lattice,
    EPS::Float64
    )::Union{Int64,Nothing}
    return findfirst(i->comp_vec(latt.R0[:,i],vec(v),EPS),1:size(latt.R0,2))
end


function find_index(
    v,
    latt::Lattice,
    EPS::Float64
    )::Union{Int64,Nothing}
    @assert length(vec(v)) == dimensions(latt)
    return find_index_nocheck(v, latt, EPS)
end



"""

    sortperm_inner_id(
        id0EqV::Vector{Int64},
        latt::Lattice,
        EPS::Float64
    )

#NOTE

===

Sort `id0EqV` according to sublattice.

"""
function sortperm_inner_id(
    id0EqV::Vector{Int64},
    latt::Lattice,
    EPS::Float64
    )
    δ1 = latt.UC.δ[:,1]
    sl(i) = latt.SL[i]
    uc_δ1_pos(i) = find_index_nocheck(latt.R0[:,i].+δ1.-latt.UC.δ[:,sl(i)], latt, EPS)
    info(id) = (uc_δ1_pos(id),sl(id),id)
    return sortperm(map(info,id0EqV))
end


function verify_monolayer(
    latt::Lattice,
    EPS::Float64        # used by find_index_nocheck()
    )::Bool
    (origin,shift) = origin_shift(latt)
    #% must be : 3D, [m,n,z] = [.., .., 1], pbc = [.., .., false]
    if ! ( dimensions(latt)==3 && latt.BBOX[3][3]==1 && (!latt.BBOX[4][3]) )
        return false
    end
    #% verify the monolayer settings
    R0_plus  = inner_site_coords(latt) .+ shift[:,3]
    R0_minus = inner_site_coords(latt) .- shift[:,3]
    find1 = [find_index_nocheck(R0_plus[:,i],  latt, EPS) for i ∈ axes(R0_plus, 2)]
    find2 = [find_index_nocheck(R0_minus[:,i], latt, EPS) for i ∈ axes(R0_minus,2)]
    test1 = Bool[ latt.EqV[f1]==-1 for f1 ∈ find1 if f1 !== nothing ]
    test2 = Bool[ latt.EqV[f2]==-1 for f2 ∈ find2 if f2 !== nothing ]
    return all(test1) && all(test2)
end


function verify_cylinder_2d(
    latt::Lattice,
    boundary::Symbol
    )::Bool
    t1 =  (boundary==:x || boundary==:y)
    if !t1 @error "Boundary not x nor y !" end
    t2 =  dimensions(latt) == 2
    if !t2 @error "Dimension not 2 !" end
    (_,_,units,pbc) = latt.BBOX
    bdn = boundary==:x ? 1 : (boundary==:y ? 2 : 3)
    t3 =  units[bdn]==1
    if !t3 @error "units[$(bdn)] = $(units[bdn]) !" end
    t4 =  (!pbc[bdn]) && (pbc[3-bdn])
    if !t4 @error "pbc is inconsistent with boundary !" end
    return (t1 && t2 && t3 && t4)
end


function verify_cylinder(
    latt::Lattice,
    boundary::Symbol,
    EPS::Float64
    )::Bool
    t0 =  verify_monolayer(latt, EPS)  #!FIXME
    if !t0 @error "Not a monolayer !" end
    t1 =  (boundary==:x || boundary==:y)
    if !t1 @error "Boundary not x nor y !" end
    t2 =  dimensions(latt) == 3
    if !t2 @error "Dimension not 3 !" end
    (origin,trans,units,pbc) = latt.BBOX
    bdn = boundary==:x ? 1 : (boundary==:y ? 2 : 3)
    t3 =  units[bdn]==1
    if !t3 @error "units[$(bdn)] = $(units[bdn]) !" end
    t4 =  (!pbc[bdn]) && (pbc[3-bdn])
    if !t4 @error "pbc is inconsistent with boundary !" end
    return (t0 && t1 && t2 && t3 && t4)
end


function check_zeroth_unitcell(
    latt::Lattice,
    inbd_test::Function
    )::Bool
    (origin,shift) = origin_shift(latt)
    invS = inv(shift)
    # the following line verifies that
    # the zeroth unit cell are totally contained
    all([ inbd_test(latt.UC.δ[:,i].-origin, invS) for i ∈ axes(latt.UC.δ,2) ])
end
