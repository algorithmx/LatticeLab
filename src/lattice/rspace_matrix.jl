@inline in01(x) = (x>=0&&x<1.0)
@inline grt1(x) = (x>=1.0)


function across_upper_border(fracRi,fracRj)
    if in01(fracRj[1]) && grt1(fracRi[1])
        return true
    elseif (in01(fracRj[1])&&in01(fracRi[1])) && (in01(fracRj[2]) && grt1(fracRi[2]))
        return true
    else
        return  ( (in01(fracRj[1])&&in01(fracRi[1]))
               && (in01(fracRj[2])&&in01(fracRi[2]))
               && (in01(fracRj[3])&&grt1(fracRi[3])) )
    end
end


function make_dic(l)
    @inline F(i) = findall(x->first(x)==i, l)
    KEYS = unique(first.(l))
    Dict(k=>Dict(v[1]=>v[2] for (p,v) in l[F(k)]) for k in KEYS)
end

#inner_j_outer_i = [(1,(2,:s)),(1,(3,:s)),(4,(2,:t)),(5,(3,:t))]
#inner_i_outer_j = [(2,(1,:s)),(3,(1,:s))]
#make_dic(inner_j_outer_i)

##


function boundary_hopping_phase_correction!(
    M,
    inner_j_outer_i,
    inner_i_outer_j,
    latt,
    BLK::Vector{UnitRange{Int64}},
    block_constructor::Function,
    block_constructor_diagonal::Function
    )
    @inline find_min_x(l) = last(findmin(map(x->first(frac(latt.R0[:,x])),l)))
    inJ_outI = make_dic(inner_j_outer_i)
    inI_outJ = make_dic(inner_i_outer_j)
    # find reference edge: j_ref -> i_ref
    j_ref = find_min_x(first.(inner_j_outer_i))
    i_ref = find_min_x(collect(keys(inJ_outI[j_ref])))
    v_ref = inJ_outI[j_ref][i_ref]
    
end

##
##
##

function rspace_matrix(
    latt::Lattice,
    BLK::Vector{UnitRange{Int64}},
    block_constructor::Function,
    block_constructor_diagonal::Function;
    correct_boundary_hopping_phase = false
    )
    (origin,shift)  = origin_shift(latt)
    invS       = inv(shift)
    @inline frac(x) = invS*(x.-origin)
    id0EqV, C  = position_index( latt.EqV )
    dimM       = sum( length.(BLK) ) 
    M          = spzeros( ComplexF64, dimM, dimM )
    X,Y,V      = findnz( latt.f )

    inner_j_outer_i = []
    inner_i_outer_j = []
    for (i0,j0,v0) ∈ zip(X,Y,V)
        if is_inner_site(latt,j0)
            if is_inner_site(latt,i0)
                M[BLK[C[i0]], BLK[C[j0]]] += block_constructor(i0,j0,v0)  #% fold
            elseif across_upper_border(frac(latt.R0[:,i0]), frac(latt.R0[:,j0]))
                push!(inner_j_outer_i, (j0,(i0,v0)))
                ieqv = equiv_site(latt,i0)
                M[BLK[C[ieqv]], BLK[C[j0]]] += block_constructor(i0,j0,v0)  #% fold
            end
        elseif is_inner_site(latt,i0) && across_upper_border(frac(latt.R0[:,j0]), frac(latt.R0[:,i0]))
            push!(inner_i_outer_j, (i0,(j0,v0)))
            jeqv = equiv_site(latt,j0)
            M[BLK[C[i0]], BLK[C[jeqv]]] += block_constructor(i0,j0,v0)  #% fold
        end
    end
    
    if correct_boundary_hopping_phase
        boundary_hopping_phase_correction!(
            M, 
            inner_j_outer_i,
            inner_i_outer_j,
            latt,
            BLK,
            block_constructor,
            block_constructor_diagonal
        )
    end

    # the diagonal
    M1 = copy(M)
    for i0 ∈ inner_site_id(latt)
        M[BLK[C[i0]], BLK[C[i0]]] += block_constructor_diagonal(i0,M1)
    end
    return M
end
