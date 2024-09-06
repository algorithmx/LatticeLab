mutable struct UnitCell <: LatticeInfo
    dim::Int64              # dimension of the lattice, 2 or 3
    nsubl::Int64            # number of sublattices in a UC
    a::Coordinates          # Basis
    δ::Coordinates          # shift vectors for sublattices
    m::Masses               # mass
    ξ::Vector{Orbits}       # atomic orbital basis
end


function isequal_upto_orbits(UC1::UnitCell,UC2::UnitCell)
    return ( UC1.dim==UC2.dim && UC1.nsubl==UC2.nsubl
          && norm(UC1.a.-UC2.a)<1e-10 && norm(UC1.δ.-UC2.δ)<1e-10
          && UC1.m==UC2.m )
end

function isequal(UC1::UnitCell,UC2::UnitCell)
    return ( isequal_upto_orbits(UC1,UC2) && UC1.ξ==UC2.ξ )
end


function ==(UC1::UnitCell,UC2::UnitCell)
    return isequal(UC1,UC2)
end

copy(uc::UnitCell) = UnitCell(  uc.dim,
                                uc.nsubl,
                                copy(uc.a),
                                copy(uc.δ),
                                copy(uc.m),
                                copy(uc.ξ), )


# ensure that all sublattices uc.δ are inside the parallelgon formed by the basis uc.a
δ_in_a(uc::UnitCell) = inbbox(uc.δ.+(uc.a*((1e-8).*rand(uc.dim))), inv(uc.a))


function check_compat(uc::UnitCell)
    c1 = size(uc.a,1)==size(uc.a,2)==size(uc.δ,1)==uc.dim
    if !c1 @error "UC dimension mismatch." end
    c2 = abs(det(uc.a))>1e-10
    if !c2 @error "UC basis orientation det<0." end
    c3 = size(uc.δ,2)==uc.nsubl
    if !c3 @error "UC sublattice list and nsubl mismatch." end
    c4 = length(uc.m)==uc.nsubl
    if !c4 @error "UC mass label numbers and nsubl mismatch." end
    c5 = length(uc.ξ)==uc.nsubl
    if !c5 @error "UC orbits numbers and nsubl mismatch." end
    c6 = δ_in_a(uc)
    if !c6 @error "UC sublattice not inside paralegon." end
    return c1 && c2 && c3 && c4 && c5 && c6
end


@inline is_UnitCell(UC) = (  Set(fieldnames(UnitCell))==Set(fieldnames(typeof(UC)))
                          && is_Coordinates(UC.a)
                          && is_Coordinates(UC.δ)
                          && is_Masses(UC.m)
                          && all(is_Orbits.(UC.ξ)) )


@inline convert(::Type{T}, uc::T) where {T<:UnitCell} = uc


function convert(::Type{T}, uc) where {T<:UnitCell}
    UnitCell( Int64(uc.dim),
              Int64(uc.nsubl),
              Coordinates(uc.a), # lazy, no type-promotion of elements
              Coordinates(uc.δ), # lazy, no type-promotion of elements
              Masses(uc.m),
              Vector{Orbits}(uc.ξ) )
end


@inline orbit_numbers(UC::UnitCell) = length.(UC.ξ)

@inline orbit_number(UC::UnitCell, subl_i::Int) = length(UC.ξ[subl_i])

@inline total_num_orbits(UC::UnitCell) = sum( orbit_numbers(UC) )

@inline num_sublattice(uc::UnitCell) = uc.nsubl

@inline mass_sublattice(uc::UnitCell, subl_i::Int) = uc.m[subl_i]

@inline reciprocal_basis(uc::UnitCell) = 2π.*inv(uc.a)'

@inline all_orbits(uc::UnitCell) = unique(vcat(uc.ξ...))

@inline select_atoms_from_uc(uc::UnitCell,atoms::Vector{Symbol}) = Int64[i for i=1:uc.nsubl if uc.m[i] ∈ atoms]

# ---------------------------------------


# accumulate a list of sizes of consecutive blocks and convert into ranges for blocks
@inline block_ranges(X) = Vector{UnitRange{Int64}}(
            map( x->(x[2]-x[1]+1):x[2],
                 zip(X,accumulate(+,X)) ) )


function generate_blocks_for_HQ(uc::UnitCell)
    blocks = block_ranges( orbit_numbers(uc) )
    dimA = sum(length.(blocks))
    return blocks, dimA
end


function generate_blocks_for_DQ(uc::UnitCell)
    blocks = block_ranges(uc.dim ⨰ uc.nsubl)
    dimA = sum(length.(blocks))
    return blocks, dimA
end
