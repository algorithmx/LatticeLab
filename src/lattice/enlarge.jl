function lattice_as_FC(
    latt0::Lattice;
    Translate = nothing,
    EPS = 1e-8
    )
    dim    = dimensions(latt0)
    Ninner = num_inner_sites(latt0)
    (origin,trans,units,pbc) = latt0.BBOX
    # extension should preserve ordering !!!
    id0EqV0, C0 = position_index(latt0.EqV)
    # sort id0EqV according to sublattice
    ASS    = sortperm_inner_id(id0EqV0, latt0, EPS)
    ASSINV = sortperm(ASS)
    id0EqV = id0EqV0[ASS]
    # map from indices of latt0.R0 to inner indices 
    C = map(x->((x>0) ? ASSINV[x] : -1), C0)

    # UC
    # the new basis is NOT determined by bbox_new
    # the new basis is determined by latt0.BBOX,
    # which records how the piece latt0 was generated
    # since latt0 is considered as a unit cell
    # the new basis is the larger bounding box
    T = (Translate===nothing) ? zeros(eltype(latt0.R0),latt0.UC.dim) : Translate
    nsubl = Ninner
    δ = latt0.R0[:,id0EqV]
    sl_inner = latt0.SL[id0EqV]
    m = latt0.UC.m[sl_inner]
    ξ = latt0.UC.ξ[sl_inner]
    a = latt0.UC.a * (trans*diagm(0=>units))
    UC = UnitCell( dim, nsubl, a, δ.+T, m, ξ )
    @assert check_compat(UC)

    # LN
    A,B,F = findnz(latt0.f)
    SP = keys(latt0.LN.SPNB)
    SPNB = Dict(sp=>Spring() for sp ∈ SP)
    for sp ∈ SP
        directions = Spring((i,j)=>[] for i ∈ 1:Ninner for j ∈ 1:Ninner)
        F_eq_sp = findall(F.==sp)
        for i ∈ F_eq_sp
            (a0,b0) = (A[i],B[i])
            if is_inner_site( latt0, b0 )
                push!(directions[(C[b0],C[a0])], latt0.R0[:,a0].-latt0.R0[:,b0])
            end
        end
        SPNB[sp] = delete_empty(directions)
    end
    LN = LinkInfo( UC, SPNB )
    @assert check_compat(LN)
    return LN
end


function enlarge(
    latt0::Lattice,
    bbox_new::BoundingBox;
    Translate = nothing,
    EPS = 1e-8
    )::Lattice
    @assert check_compat(latt0)
    #NOTE extension should preserve sublattice ordering !!!
    FC_new = lattice_as_FC(latt0, EPS=EPS, Translate=Translate)
    build_lattice(FC_new, bbox_new)
end


function enlarge(
    uc0::UnitCell,
    bbox_new::BoundingBox;
    EPS = 1e-8
    )::UnitCell
    @assert check_compat( uc0 )
    latt0 = build_lattice(LinkInfo(uc0,Dict{Symbol,Spring}()), bbox_new)
    #NOTE extension should preserve sublattice ordering !!!
    ln = lattice_as_FC(latt0, EPS=EPS, Translate=nothing)
    return ln.UC
end
