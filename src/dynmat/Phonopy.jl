#TODO 
#+ understand precisions and decimals in Phonopy
#+ symprec
#+ degeneracy_tolerance
#+ decimals
#+ dm_decimals
#+ dynamical_matrix_decimals
#+ fc_decimals
#+ force_constants_decimals


mutable struct PhonopyCopy
    UC::LatticeInfo
    SUC::LatticeInfo
    MASS::Dict
    FC::Array{Float64,4}
    P2S
    S2P
    P2P
    SMVEC::Array{Float64,4}
    MULT::Array{Int,2}
end


@inline fake_orbits(s)  = [:X,]

@inline make_UC(x) = UnitCell(  3,  x.get_number_of_atoms(), 
                                    x.get_cell()', x.get_positions()', 
                                    Symbol.(x.get_chemical_symbols()), 
                                    fake_orbits.(x.get_chemical_symbols())  )

function to_mat_3x3(s)
    L = Float64.(eval.(Meta.parse.(split(s," ",keepempty=false))))
    return [L[1:3], L[4:6], L[7:9]]
end


function phonony_from_poscar(
    phonopy,
    poscar_filename::String,
    conf_filename::String,
    force_constants_filename::String
    )
    I3 = to_mat_3x3("1 0 0 0 1 0 0 0 1")
    conf = phonopy.cui.settings.PhonopyConfParser(filename=conf_filename)
    primitive_matrix = ("primitive_axis" ∈ keys(conf.confs)) ? (conf.confs["primitive_axis"] |> to_mat_3x3) : I3
    supercell = ("dim" ∈ keys(conf.confs)) ? (conf.confs["dim"] |> to_mat_3x3) : I3
    if force_constants_filename==""
        unitcell = phonopy.interface.vasp.read_vasp(poscar_filename)
        return  phonopy.Phonopy(
                            unitcell, 
                            supercell, 
                            primitive_matrix=primitive_matrix
                            ) 
    else
        return  phonopy.load(
            unitcell_filename=poscar_filename, 
            supercell_matrix=supercell, 
            primitive_matrix=primitive_matrix,
            force_constants_filename=force_constants_filename)
    end
end


function PhonopyCopy_from_vasp(
    phonopy,
    poscar_filename::String,
    force_filename::String,
    conf_filename::String
    )::PhonopyCopy
    phonon = phonony_from_poscar(
                    phonopy, 
                    poscar_filename, 
                    conf_filename, 
                    (endswith(force_filename,"FORCE_CONSTANTS") ? force_filename : "") )
    prim = phonon.get_primitive()
    (smallest_vectors,multiplicity) = prim.get_smallest_vectors()
    (sa,sb,sc,sd) = size(smallest_vectors)
    smallest_vectors_cartesian = permutedims(
        reshape(reshape(smallest_vectors,(sa*sb*sc,sd)) * prim.get_cell(), (sa,sb,sc,sd)), 
        [4,3,2,1])

    UC  = make_UC(prim)
    SUC = make_UC(phonon.get_supercell())

    if !endswith(force_filename,"FORCE_CONSTANTS") && endswith(force_filename,"FORCE_SETS")
        # FORCE_SETS should be present after VASP calculation
        # produce force constants
        #* unit:
        #* 
        force_sets = phonopy.file_IO.parse_FORCE_SETS(filename=force_filename)
        phonon.set_displacement_dataset(force_sets)
        phonon.produce_force_constants()
    end

    mass = prim.get_masses()
    symb = prim.get_chemical_symbols()
    MASS = Dict(k=>v for (k,v) ∈ unique(zip(Symbol.(symb),mass)))
    fc = permutedims(phonon._force_constants, [4,3,1,2])  #!!!!!
    return PhonopyCopy( UC, SUC, MASS, fc, 
                        prim.p2s_map, prim.s2p_map, prim.p2p_map, 
                        smallest_vectors_cartesian, multiplicity  )
end


function extract_force_constants(
    PH::PhonopyCopy, 
    LN::Dict; 
    SVEC_NORM_EPS=1e-8,
    SVEC_PREC_CTRL=8, 
    FC_PREC_CTRL=8, 
    FC_EPS=1e-6,
    )

    #: ----------------------------------------------------
    # from phonopy :
    # 
    # check class ShortestPairs in cell.py
    # in particular _run_sparse()
    # and finally phpy_set_smallest_vectors_sparse() in phonopy.c
    #
    #: ----------------------------------------------------
    # get_smallest_vectors seeks the shortest distance of a pair of sites
    # from j in primitive cell to i in supercell
    # by adding appropriate Bravais lattice vectors lattice_points[k]
    # to the position difference R_ij = pos_to[i] - pos_from[j]
    # Bravais lattice vectors introduces a multiplicity
    # which is recorded in multiplicity[to_i][from_j] 
    #: ----------------------------------------------------
    #
    # meaning of vec, len
    #    R_ij = pos_to[i] - pos_from[j]
    #    vec[k] = R_ij + lattice_points[k]
    #    len[k] = norm(reduced_basis * vec[k])
    # meaning of smallest_vectors
    #    c = 0
    #    if  (len[k] - min_len < symprec)
    #       smallest_vectors[to_i * num_pos_from + from_j][c] = trans_mat * vec[k]
    #       c += 1
    #    end
    # or 
    #    c = 0
    #    if  (len[k] - min_len < symprec)
    #       smallest_vectors[i][j][c] = trans_mat * vec[k]
    #       c += 1
    #    end
    # meaning of multiplicity (C-trick on array indices are used)
    #    multiplicity[to_i * num_pos_from + from_j] = count
    # or 
    #    multiplicity[to_i][from_j] = count
    #: ----------------------------------------------------

    maxabs3(x) = maximum(abs.(x))
    SP = []
    NB = []
    for (i, s_i) in enumerate(PH.P2S)
        for (j, s_j) in enumerate(PH.P2S)
            @inline find_v(svec_phonopy) = findall(x->maxabs3(svec_phonopy-x[2])<SVEC_NORM_EPS, LN[(i,j)])
            SP_local = []
            NB_local = []
            for k = 1:length(PH.SUC.m)
                if maxabs3(PH.FC[:,:,s_i+1,k])>FC_EPS
                    if s_j == PH.S2P[k]
                        multi = PH.MULT[k,i]  #+ in phonopy.c ?
                        for v = 1:multi
                            svec = PH.SMVEC[:,v,i,k]  #+ in phonopy.c ?
                            #+ is it possible that svec = [0,0,0] ?
                            if maxabs3(svec)>SVEC_NORM_EPS   #*** PRECISION CONTROL ***
                                links = [LN[(i,j)][p] for p in find_v(svec)]  #*** PRECISION CONTROL ***
                                if length(links)!=1
                                    msg = [ "extract_force_constants() error ",
                                            "(i, s_i) = ($i, $s_i) (j, s_j) = ($j, $s_j)",
                                            "k = $k, v = $v,",
                                            "links = $(links)." ]                                    
                                    throw(SmallVecMismatchError(join(msg,"\n\t\t")))
                                end
                                push!(  SP_local, 
                                        links[1][1] => 
                                        round.( (1.0/multi).*PH.FC[:,:,s_i+1,k], digits=FC_PREC_CTRL ))  #*** PRECISION CONTROL ***
                                push!(  NB_local, 
                                        links[1][1] => 
                                        round.( svec, digits=SVEC_PREC_CTRL ))  #*** PRECISION CONTROL ***
                            end
                        end
                    end
                end
            end
            if length(SP_local)>0  push!(SP, (i,j)=>Dict(SP_local))  end
            if length(NB_local)>0  push!(NB, (i,j)=>Dict(NB_local))  end
        end
    end
    if  !allunique(vcat([collect(keys(v)) for (k,v) in SP]...))
        throw(KeysNotUniqueError("SP"))
    end
    if  !allunique(vcat([collect(keys(v)) for (k,v) in NB]...))
        throw(KeysNotUniqueError("NB"))
    end
    neighbors = [link_symbol=>Spring(ij=>[dir,])  for (ij,L) in Dict(NB) for (link_symbol,dir) in L]
    force3x3  = [link_symbol=>(:FC3X3, [fc,])     for (ij,L) in Dict(SP) for (link_symbol,fc ) in L]
    return Dict(force3x3),  LinkInfo(copy(PH.UC), Dict(neighbors))
end


##------------------------------------------------------

##------------------------------------------------------