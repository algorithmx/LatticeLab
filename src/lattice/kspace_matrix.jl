# useful functions
function get_mn(dist,basis,inv_basis)
    mn = Int64.(round.(inv_basis*dist))
    # the following means that dist must be
    # a Bravais lattice vector
    diff = norm( basis*mn .- dist )
    if diff > 1e-6
        @warn "norm(basis*mn.-dist) = $(diff) too large."
    end
    @assert diff < 1e-4
    return mn
end


"""
    kspace_matrix(
        latt::Lattice,
        BLK::Vector{UnitRange{Int64}},
        block_constructor::Function,
        block_constructor_diagonal::Function,
        inbd_test::Function
    )


###Note:

============

Variables and their meanings:

`MN` : 

`(origin,shift)` : 

`SP` : 

`NB` : 

`all_dR_Φ` :  List of (`i`,`sl`,`MN`,`ϕ`) where 

`i` is the index in `SP` (from the keys of `latt.LN.SPNB`)

`sl = (jj,ii)` is a pair of sublattice id with the i'th type of link attached

`MN` is the shift of the 2nd site in Bravais lattice to the 1st site, with respect to direction `d ∈ NB[iln][(jj,ii)]`

`ϕ` is the block matrix

"""
function kspace_matrix(
    latt::Lattice,
    BLK::Vector{UnitRange{Int64}},
    block_constructor::Function,
    block_constructor_diagonal::Function,
    inbd_test::Function
    )
    dim   = dimensions(latt)
    Nsubl = num_sublattice(latt)
    @assert length(BLK) == Nsubl
    a     = latt.UC.a
    ainv  = inv(a)
    δ     = latt.UC.δ
    nlnk  = length(latt.LN.SPNB)
    SP    = sck(latt.LN.SPNB)
    NB    = [ latt.LN.SPNB[sp] for sp ∈ SP ]
    Norbits = sum(length.(BLK))
    (origin,shift) = origin_shift(latt)
    invS  = inv(shift)

    #NOTE CONVENTION 1<--2
    # H(1,2) C†(1) C(2)
    # the following list connsists of pairs (i, sl, mn, mat)
    # where i is the index in the keys of latt.LN.SPNB
    # sl is a pair of sublattice id with the i'th type of link attached
    # mn is the shift of the 2nd site in Bravais lattice to the 1st site, with respect to direction d
    # and the last element mat is the corresponding block matrix
    all_dR_Φ =[( iln, (jj,ii), get_mn(d.-(δ[:,ii].-δ[:,jj]), a, ainv),
                 BIJ(BLK[ii], BLK[jj], block_constructor(ii,jj,d,SP[iln]), Norbits) )
                 for iln ∈ 1:nlnk
                     for (jj,ii) ∈ keys(NB[iln])
                         for d ∈ NB[iln][(jj,ii)]    #XXX CONVENTION jj-->ii
                             if ( inbd_test(δ[:,jj].+d.-origin, invS) )] # see next line #!
    # the function inbd_test(r, A) test 0.0 <= A*r < 1.0
    # it is used for open boundary conditions (at certain direction)
    # it is always true for one direction if pbc==true for that direction
    # for example, a cylinder has open bd. cond. in the x direction
    # and closed bd. cond. in the other direction y

    for sl ∈ 1:Nsubl
        push!( all_dR_Φ,
               (1, (sl,sl), 0 ⨰ dim,
                    BIJ(BLK[sl], BLK[sl], block_constructor_diagonal(sl), Norbits) ) )
    end

    MN_unique = unique([ MN for (i,p,MN,ϕ) ∈ all_dR_Φ ])

    mat = Dict( Vector{Int64}(mnuniq) => sum( ϕ for (i,p,MN,ϕ) ∈ all_dR_Φ if MN==mnuniq )
                for mnuniq ∈ MN_unique 
              ) |> delete_empty_sparse

    if length(mat)==0
        mat = Dict(Vector{Int64}(0 ⨰ dim)=>sparse(Int[],Int[],ComplexF64[],Norbits,Norbits))
    end

    return mat
end
