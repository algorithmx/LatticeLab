index_t(t, N_min::Vector{Int}, N_max::Vector{Int}, dim::Int)::Int = 1+sum((t[d]-N_min[d])*prod(Int[(N_max[l]-N_min[l]+1) for l=1:d-1]) for d=1:dim)

index_t0(  N_min::Vector{Int}, N_max::Vector{Int}, dim::Int)::Int = 1+sum((-N_min[d])*prod(Int[(N_max[l]-N_min[l]+1) for l=1:d-1]) for d=1:dim)

iter_t(N_min::Vector{Int}, N_max::Vector{Int}, dim::Int) = Iterators.product([N_min[k]:N_max[k] for k=1:dim]...)


#: =====================================================
#: compute_Bravais_cutoff

@inline max3(x) = maximum(Int.(ceil.(abs.(x))))

@inline collect_vec_v1(SPNB) = vcat([vcat(collect(values(nbd))...) for nbd ∈ values(SPNB)]...)

## EXTRACT ALL DIRECTION VECTORS IN  fc.SPNB
## COMPUTE THE MAXIMUM OF THEIR FRACTIONAL COORDINATES
#TODO tighten the margin for speed !!!
@inline compute_Bravais_cutoff(fc::LinkInfo, inv_basis) = 
                ( length(fc.SPNB)>0
                    ? maximum(map(v->max3(inv_basis*v), unique(collect_vec_v1(fc.SPNB))))
                    : 3 )  #* shouldn't be 0; margin is to reserve space for lattice translations

#: =====================================================
#: generate_R0

function generate_R0_compute_EqV(
    geom::UnitCell,
    bbox::BoundingBox,
    margin::Int64
    )::Tuple{Coordinates,Vector{Int},Vector{Int64}}

    # prepare
    (a, δ, Nsubl, dim) = (geom.a, geom.δ, geom.nsubl, geom.dim)
    pbc = bbox[4]
    (origin, shift)    = origin_shift(bbox, a)
    invS         = inv(shift)
    N_min, N_max = bounding_box_Nmin_Nmax(a, bbox, margin)
    Trans        = bbox[2]*diagm(0=>bbox[3])  #! supercell trans Bravais vects

    @info "generate_R0_compute_EqV() : \n\t\t bbox = $(bbox)\n\t\t N_min, N_max = $((N_min,N_max))"
    @assert all(N_min .<= 0) && all(N_max .>= 0)
    iter_N_min_max = iter_t(N_min, N_max, dim)  # regardless of bd conds
    all_t = collect(iter_N_min_max)
    # generate all sites regardless of the boundary conditions
    all_sites = hcat([  # r = 
                        δ[:,s].+(a*[t...,]) 
                        for t ∈ iter_N_min_max for s ∈ 1:Nsubl ]...)
    sublattice_index = [ s     for t ∈ iter_N_min_max for s ∈ 1:Nsubl ]
    translation_vecs = [ t     for t ∈ iter_N_min_max for s ∈ 1:Nsubl ]
    subl_trans       = [ (s,t) for t ∈ iter_N_min_max for s ∈ 1:Nsubl ]
    Nsites = length(sublattice_index)
    ids    = Dict((t,s)=>i for (i,(t,s)) in enumerate(zip(translation_vecs,sublattice_index)))

    @info "generate_R0_compute_EqV() : compute EqV ..."
    print(""); flush(stdout);
    #TODO bottleneck
    EqV = fill(-1,Nsites)
    @time begin
        @inline inside_bbox(s,t) = (t ∈ all_t && s ∈ 1:Nsubl) && inbbox(all_sites[:,ids[(t,s)]] .- origin,invS)
        #inside = Dict((t,s)=>i for (i,(t,s)) ∈ enumerate(zip(translation_vecs, sublattice_index))  if inside_bbox(s,t))
        inside = Dict((t,s)=>i for (i,(s,t)) ∈ enumerate(subl_trans)  if inside_bbox(s,t))
        SHIFTS = collect(Iterators.product([-(7*Int(pbc[d])):(7*Int(pbc[d])) for d=1:dim]...)) ;
        SHIFTS_ABS = map(v->Trans*[v...,], SHIFTS);
        p0  = findfirst(v->all(v.==0), SHIFTS)
        ii  = 1
        for t ∈ iter_N_min_max
            for s ∈ 1:Nsubl
                #p = findfirst(v->inside_bbox(s,((t.+Trans*[v...,])...,)), SHIFTS)
                p = findfirst(v->inside_bbox(s,((t.+v)...,)), SHIFTS_ABS)
                if p!==nothing
                    #EqV[ii] = (p==p0 ? 0 : inside[(t.+(Trans*[SHIFTS[p]...,]...,),s)])
                    EqV[ii] = (p==p0 ? 0 : inside[(t.+(SHIFTS_ABS[p]...,),s)])
                end
                ii += 1
            end
        end
    end
    print(""); flush(stdout);
    return all_sites, sublattice_index, EqV
end


#: =====================================================
#: generate_f_compute_EqV


function reorganize_SPNB(SPNB, Nsubl::Int)::Dict
    dic = Dict(ij=>[] for ij ∈ Iterators.product(1:Nsubl,1:Nsubl))
    for (label, dir_dict) ∈ SPNB
        for (ij,dirs) ∈ dir_dict
            push!(dic[(ij[1],ij[2])], (label,dirs))
        end
    end
    return copy(dic)
end


function check_conflicts(I,J,V)
    conflicts = []
    IJVsort = sort(collect(zip(I,J,V)))
    for l=1:length(IJVsort)-1
        if IJVsort[l][1:2]==IJVsort[l+1][1:2] && IJVsort[l][3]!=IJVsort[l+1][3]
            push!(conflicts, IJVsort[l][1:2], (IJVsort[l][3],IJVsort[l+1][3]))
        end
    end
    return conflicts
end



"""
    function generate_f_compute_EqV(
        LnInfo::LinkInfo, 
        R0::Array{T,2},
        eqv0::Vector{Int}, 
        bbox::BoundingBox,
        margin::Int, 
        Nsites::Int
        ) where {T<:Real}


### PSEUDOCODE

```
for each sublattice pair (iδ,jδ)
    for each label
        for each direction
            # case I
            for each site jj in bounding box on sublattice jδ
                compute all ii such that  direction .+ R0[:,jj] == R0[:,ii]
                check
                record (ii,jj,label)
            end
            # case II
            for each site ii in bounding box on sublattice iδ
                compute all jj such that   R0[:,jj] == R0[:,ii] .+ direction
                check
                record (ii,jj,label)
            end
        end
    end
end
```

"""
function generate_f(
    LnInfo::LinkInfo, 
    R0::Array{T,2},
    bbox::BoundingBox,
    margin::Int, 
    Nsites::Int
    ) where {T<:Real}

    (Nsubl, dim, a, δ) = (LnInfo.UC.nsubl, LnInfo.UC.dim, LnInfo.UC.a, LnInfo.UC.δ)
    inva = inv(a)
    N_min, N_max = bounding_box_Nmin_Nmax(a, bbox, margin)
    #ig0 = index_t0(N_min, N_max, dim)
    (origin, shift) = origin_shift(bbox, a)
    invS = inv(shift)

    # helper function
    @inline trans(x) = Int.(round.(inva*x))
    @inline check(a,b,x) = norm((R0[:,a].-R0[:,b]).-x) < 1e-5

    # boundary condition check
    boundaries = findall(bbox[4].==false)
    @inline in_bbox_check(r) = inbboxt(r.-origin,invS,boundaries)

    # sparse matrix to return
    I=Int64[]; J=Int64[]; V=keytype(LnInfo.SPNB)[];

    # inner id
    ind0 = [iv for iv=1:size(R0,2) if inbbox(R0[:,iv].-origin,invS)]
    innerid = Dict(jδ=>ind0[findall(x->mod(x-1,Nsubl)==jδ-1,ind0)] for jδ in 1:Nsubl) #! 1-based

    # re-organize SPNB dictionary : (iδ,jδ)=>(label,directions)
    IJFD = reorganize_SPNB(LnInfo.SPNB, Nsubl)
    @inline outside(t) = (any(t.-N_min.<0) || any(N_max.-t.<0))

    #! main loop
    @info "generate_f() : main loop ..."
    @time begin
        #% for each sublattice pair (iδ,jδ)
        for (iδ,jδ) ∈ Iterators.product(1:Nsubl,1:Nsubl)
            # for (i,j) ∈ (I,J), either i or j in bounding box
            #% for each label
            for (label, directions) ∈ IJFD[(jδ,iδ)]
                # (jδ,iδ) NOT (iδ,jδ)
                # because  CONVENTION j-->i for Lattice.f[i,j]
                # and      CONVENTION s-->t for spring[(s,t)]
                #% for each direction
                for d ∈ directions
                    d_m_δi = d.-δ[:,iδ]
                    d_p_δj = d.+δ[:,jδ]
                    #% case I : jδ = j%t in bounding box
                    # R[j] + d = R[i] = δ[i%t] + t[i]
                    # t[i] = R[j] + d - δ[i%t]
                    #% for each site jj in bounding box on sublattice jδ
                    for j0 ∈ innerid[jδ]
                        #% compute all ii such that   direction .+ R0[:,jj] == R0[:,ii]
                        Bravais = trans(d_m_δi.+R0[:,j0])
                        if outside(Bravais)  continue  end
                        i1 = iδ+(index_t(Bravais,N_min,N_max,dim)-1)*Nsubl
                        #% check
                        if !check(i1,j0,d)
                            throw(error("(iδ,jδ,i1,j0)=($iδ,$jδ,$i1,$j0) incompatible ! dR=$(R0[:,i1].-R0[:,j0]), dir=$(d)"))
                        end
                        #% record
                        if in_bbox_check(R0[:,i1])
                            push!(I,i1); push!(J,j0); push!(V,label);
                        end
                    end
                    #% case II : iδ = i%t in bounding box
                    # R[j] + d = δ[j%t] + t[j] + d = R[i]
                    # t[j] = R[i] - d - δ[j%t]
                    #% for each site ii in bounding box on sublattice iδ
                    for i0 ∈ innerid[iδ]
                        #% compute all jj such that   R0[:,jj] == R0[:,ii] .+ direction
                        Bravais = trans(R0[:,i0].-d_p_δj)
                        if outside(Bravais)  continue  end
                        j1 = jδ+(index_t(Bravais,N_min,N_max,dim)-1)*Nsubl
                        #% check
                        if !check(i0,j1,d)
                            throw(error("(iδ,jδ,i0,j1)=($iδ,$jδ,$i0,$j1) incompatible ! dR=$(R0[:,i0].-R0[:,j1]), dir=$(d)"))
                        end
                        #% record
                        if in_bbox_check(R0[:,j1])
                            push!(I,i0); push!(J,j1); push!(V,label);
                        end
                    end
                end
            end
        end
    end
    print(""); flush(stdout);
    conflicts = check_conflicts(I,J,V)
    ⊛(a,b) = b
    return sparse(I, J, V, Nsites, Nsites, ⊛), conflicts
end


#: =====================================================
#: build_lattice


function build_lattice(
    LnInfo::LinkInfo,
    bbox::BoundingBox
    )::Lattice

    @assert check_compat(bbox)   "[ERROR] build_lattice() bbox incompatible."
    @assert check_compat(LnInfo) "[ERROR] build_lattice() LnInfo incompatible."

    #: margin = maximum of fractional coordinates of 
    #: the direction vectors in LnInfo.SPNB
    #TODO tighten the margin for speed !!!
    margin  = compute_Bravais_cutoff(LnInfo, inv(LnInfo.UC.a))
    @info "build_lattice() : \n\t\t margin = $(margin)"
    print(""); flush(stdout);

    R0, Subl, EqV = generate_R0_compute_EqV(LnInfo.UC, bbox, margin)
    Nsites  = size(R0,2)
    print(""); flush(stdout);

    f, conflicts  = generate_f(LnInfo, R0, bbox, margin, Nsites)
    @assert  length(conflicts)==0   "conflicts = $conflicts"
    print(""); flush(stdout);

    _DEBUG_MODE_ ? (length(conflicts)>0  &&  @warn("\n"*join(string.(conflicts), "\n"))) : @assert(length(conflicts)==0)

    # generate index array of permutation for lattice translation
    perm  = translate_perm_0( R0, LnInfo.UC.a*bbox[2] )
    Nperm = translate_perm_0( R0, LnInfo.UC.a*bbox[2] )

    return Lattice( R0,
                    bbox,
                    copy(LnInfo),
                    Subl,
                    f,
                    copy(LnInfo.UC),
                    EqV,
                    perm, 
                    Nperm )
end

