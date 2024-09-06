#NOTE permutation[i,d] is the index of v = ( R0[:,i] + basis[:,d] ) in R0
#NOTE permutation[i,d] is meaningfuul only when EqV[i]==0 (inner sites)
#NOTE permutation[permutation] = translate twice
# ------------------------------------------

function translate_perm(
    R0::Coordinates,
    R0P,
    trans::Coordinates, # unitcell / supercell basis
    EqV::Vector{Int64},
    pbc::Vector{Bool},
    EPS::Float64
    )::Indices

    # this function is caclled by build_lattice() and apply_periodic_boundary_condition()
    (dim,Nsites) = (size(trans,2),size(R0,2))

    # permutation[:,1] --> Tx ;  permutation[:,2] --> Ty ; permutation[:,3] --> Tz
    # println("[INFO --- build_lattice()] generate permutations for translation ...")
    permutation = fill( -one(Int64), Nsites, dim )
    for d ∈ 1:dim
        if pbc[d]
            A = trans[:,d]
            for k ∈ 1:Nsites
                if EqV[k]==0 # for all inner sites
                    p = index_R0( R0[:,k].+A, R0, R0P, EPS )

                    # +A translate the inner site out of the lattice
                    if p <= 0
                        println((k, R0[:,k], p, R0[:,k].+A))
                        @assert p > 0  "+A translate the inner site out of the lattice!"
                    end

                    # +A moves site i too far from the bounding box
                    if EqV[p]<0
                        println((k, R0[:,k], p, R0[:,k].+A, EqV[p]))
                        @assert EqV[p] >= 0   "+A moves site i too far from the bounding box" #TODO
                    end

                    # EqV[p]==0 : +A still inside
                    # EqV[p] >0 : +A move the site just across the boundary
                    permutation[k,d] = (EqV[p]==0) ? p : EqV[p]
                else
                    # if site k is not inside bounding box, don't care
                    permutation[k,d] = k
                end
            end
        end
    end
    return permutation
end


# private
function translate_perm_0(
    R0::Coordinates,
    trans::Coordinates, # unitcell / supercell basis
    )::Indices

    # this function is caclled by build_lattice() and apply_periodic_boundary_condition()
    (dim,Nsites) = (size(trans,2),size(R0,2))

    # permutation[:,1] --> Tx ;  permutation[:,2] --> Ty ; permutation[:,3] --> Tz
    # println("[INFO --- build_lattice()] generate permutations for translation ...")
    permutation = fill( -one(Int64), Nsites, dim )
    return permutation
end

# private
@inline make_sense( PERM::Index,
                    EqV::Index ) = all( PERM[EqV.!=0].==collect(1:length(PERM))[EqV.!=0] )


function translate_dN_units(
    dN::Vector{Int64},
    latt::Lattice
    )::Vector{Int64}
    (dim,Nsites) = lattice_dimensions(latt)
    @assert check_compat(latt)
    @assert dim==length(dN)
    @assert all([( make_sense( latt.TransPerm[:,k], latt.EqV )
                && make_sense( latt.TransNegPerm[:,k], latt.EqV ))
                for k ∈ 1:dim if (dN[k]!=0 && latt.BBOX[4][k])
                ]) "[ERROR] translate_dN_units() doesn't make sense !!! Check the PBC !!!!!"
    permutation = collect(1:Nsites)
    for d ∈ 1:dim
        if latt.BBOX[4][d]  # pbc
            if dN[d]>0
                for r ∈ 1:dN[d]
                    permutation[:] = copy(permutation[latt.TransPerm[:,d]])
                end
            elseif dN[d]<0
                for r ∈ 1:(-dN[d])
                    permutation[:] = copy(permutation[latt.TransNegPerm[:,d]])
                end
            end
        end
    end
    return permutation
end


# matrix transformation acting on D(r1, s1; r2, s2)
function Ta(latt::Lattice; d=1)
    tp = nothing
    if d==1
        tp = translate_dN_units( (latt.UC.dim==2 ? [1,0] : [1,0,0]), latt )
    elseif d==2
        tp = translate_dN_units( (latt.UC.dim==2 ? [0,1] : [0,1,0]), latt )
    elseif d==3 && latt.UC.dim==3
        tp = translate_dN_units( [0,0,1], latt )
    else
        nothing
    end
    @assert make_sense( tp, latt.EqV )
    id0EqV, C = position_index(latt.EqV) # index of inner sites in R0
    ND = length(id0EqV)
    # this line is difficult
    # C[tp[id0EqV[k]]] means the k'th index in EqV with EqV[id]==0
    # translated by a lattice vector, to index i <-- tp[id0EqV[k]]
    # this index i is the C[i]'th index in EqV with EqV[i]==0 (inner site)
    id0EqV_translated = [ C[tp[id0EqV[k]]] for k ∈ 1:ND ]
    @assert all(id0EqV_translated.>0) && ispermute(id0EqV_translated)
    # C[i] == -1 when eqv[i] < 0
    P = sum( Eij(id0EqV_translated[k],k,ND) for k ∈ 1:ND )
    return kron( P, spdiagm(0=>(latt.UC.dim==2 ? [1,1] : [1,1,1])) )
end

TX(latt::Lattice) = Ta( latt, d=1 )
TY(latt::Lattice) = Ta( latt, d=2 )
TZ(latt::Lattice) = Ta( latt, d=3 )
