#TODO 
#+ margin

@inline rotate_2d_by_degrees(v,deg) = ( Float64( cos(π*(deg//180))*v[1] - sin(π*(deg//180))*v[2] ),
                                        Float64( cos(π*(deg//180))*v[2] + sin(π*(deg//180))*v[1] ) )


function find_bounding_box(V)
    (vx_min, vx_max) = (100000.0, -100000.0)
    (vy_min, vy_max) = (100000.0, -100000.0)
    for v in V
        vx_max  = max(v[1],vx_max)
        vx_min  = min(v[1],vx_min)
        vy_max  = max(v[2],vy_max)
        vy_min  = min(v[2],vy_min)
    end
    return (vx_min, vx_max, vy_min, vy_max)
end


function lattice_to_graph(
    latt0::Lattice{T1,TR},
    view_axes::Vector{Int64},
    bounds6
    ) where { T1<:AbstractMatrix,TR<:Real }

    latt = convert(Lattice{T1,Float64}, latt0)

    (dim,Nsites) = lattice_dimensions(latt)
    @assert ( ((view_axes==[1,2] || view_axes==[2,3] || view_axes==[1,3]) && dim==3) ||
              ((view_axes==[1,2]) && dim==2) )

    #
    latt = convert(Lattice{T1,Float64}, latt0)
    # find inner and boundary sites / trim the R0 site coordinate array
    (xmin,xmax,ymin,ymax,zmin,zmax) = bounds6
    @inline in_bound(v) = ((v[1]>=xmin && v[1]<=xmax)
                        && (v[min(2,dim)]>=ymin && v[min(2,dim)]<=ymax)
                        && (v[min(3,dim)]>=zmin && v[min(3,dim)]<=zmax))
    A0 = outer_border_and_inner_site_id(latt)
    A = [ i for i ∈ A0 if in_bound(latt.R0[:,i]) ]
    C = Dict(k=>i for (i,k) ∈ enumerate(A)) # position index
    Ain = inner_site_id(latt)
    Cin = Dict(s=>i for (i,s) ∈ enumerate(Ain))
    Ninner = length(A)

    # edges
    I,J,F = findnz(latt.f)
    IJ = collect(zip(I,J))
    # select edges between two sites
    # with at least one inner site
    E = [(C[i],C[j]) for (i,j)   ∈ zip(I,J)   if (i ∈ A && j ∈ A)]
    F = [f           for (i,j,f) ∈ zip(I,J,F) if (i ∈ A && j ∈ A)]

    # vertices
    V = Array(latt.R0[view_axes,A])
    M = latt.UC.m[latt.SL[A]]

    return V,M, E,F, A,C, Ain,Cin
end


function lattice_to_graph_V(
    latt0::Lattice{T1,TR},
    view_axes::Vector{Int64},
    bounds6
    ) where { T1<:AbstractMatrix,TR<:Real }
    (dim,Nsites) = lattice_dimensions(latt0)
    @assert ( ((view_axes==[1,2] || view_axes==[2,3] || view_axes==[1,3]) && dim==3) ||
              ((view_axes==[1,2]) && dim==2) )

    #
    latt = convert(Lattice{T1,Float64}, latt0)
    # find inner and boundary sites / trim the R0 site coordinate array
    (xmin,xmax,ymin,ymax,zmin,zmax) = bounds6
    @inline in_bound(v) = ((v[1]>=xmin && v[1]<=xmax)
                        && (v[min(2,dim)]>=ymin && v[min(2,dim)]<=ymax)
                        && (v[min(3,dim)]>=zmin && v[min(3,dim)]<=zmax))
    A0 = outer_border_and_inner_site_id(latt)
    A = [ i for i ∈ A0 if in_bound(latt.R0[:,i]) ]
    V = Array(latt.R0[view_axes,A])
    return V
end


function lattice_to_graph(
    latt0::Lattice{T1,TR},
    M::MillerIndex,
    bounds6
    ) where { T1<:AbstractMatrix,TR<:Real }
    (dim,Nsites) = lattice_dimensions(latt0)
    @assert is_MillerIndex(M)

    # find inner and boundary sites / trim the R0 site coordinate array
    (xmin,xmax,ymin,ymax,zmin,zmax) = bounds6
    @inline in_bound(v) = ((v[1]>=xmin && v[1]<=xmax)
                        && (v[min(2,dim)]>=ymin && v[min(2,dim)]<=ymax)
                        && (v[min(3,dim)]>=zmin && v[min(3,dim)]<=zmax))
    #
    latt = convert(Lattice{T1,Float64}, latt0)
    A0 = outer_border_and_inner_site_id(latt)
    A = [ i for i ∈ A0 if in_bound(latt.R0[:,i]) ]
    C = Dict(k=>i for (i,k) ∈ enumerate(A)) # position index
    Ain = inner_site_id(latt)
    Cin = Dict(s=>i for (i,s) ∈ enumerate(Ain))
    Ninner = length(A)

    # edges
    I,J,F = findnz(latt.f)
    IJ = collect(zip(I,J))
    # select edges between two sites
    # with at least one inner site
    E = [(C[i],C[j]) for (i,j)   ∈ zip(I,J)   if (i ∈ A && j ∈ A)]
    F = [f           for (i,j,f) ∈ zip(I,J,F) if (i ∈ A && j ∈ A)]

    # vertices
    V = Array(to_projector(M,latt.UC.a) * latt.R0[:,A])
    M = latt.UC.m[latt.SL[A]]

    return V,M, E,F, A,C, Ain,Cin
end


function lattice_to_graph_V(
    latt0::Lattice{T1,TR},
    M::MillerIndex,
    bounds6
    ) where { T1<:AbstractMatrix,TR<:Real }
    (dim,Nsites) = lattice_dimensions(latt0)
    @assert is_MillerIndex(M)

    # find inner and boundary sites / trim the R0 site coordinate array
    (xmin,xmax,ymin,ymax,zmin,zmax) = bounds6
    @inline in_bound(v) = ((v[1]>=xmin && v[1]<=xmax)
                        && (v[min(2,dim)]>=ymin && v[min(2,dim)]<=ymax)
                        && (v[min(3,dim)]>=zmin && v[min(3,dim)]<=zmax))
    #
    latt = convert(Lattice{T1,Float64}, latt0)
    A0 = outer_border_and_inner_site_id(latt)
    A = [ i for i ∈ A0 if in_bound(latt.R0[:,i]) ]
    # vertices
    V = Array(to_projector(M,latt.UC.a) * latt.R0[:,A])
    return V
end


# the order of transformation is important to align images with margin
# shift S at final stage
function transform_vertiecs(V::Array{Float64,2}, Λ, S, θ, H, MAR)
    V1 = [ Tuple(rotate_2d_by_degrees(Λ.*V[:,i],θ)) for i ∈ 1:size(V,2) ]
    (vx_min, vx_max, vy_min, vy_max) = find_bounding_box(V1)
    ORIGIN = (S == nothing) ? (-vx_min, -vy_min) : S
    YMAX = H > 0 ? H : (vy_max-vy_min)
    @inline shiftOflipY(t) = ( t[1]+ORIGIN[1]+MAR, YMAX-(t[2]+ORIGIN[2])+MAR )
    return shiftOflipY.(V1), ORIGIN, YMAX
end


function transform_vertiecs(V::Vector{Tuple{Float64,Float64}}, Λ, S, θ, H, MAR)
    transform_vertiecs(hcat([[t[1],t[2]] for t ∈ V]...), Λ, S, θ, H, MAR)
end


function transform_vertiecs(V::Vector{Vector{Float64}}, Λ, S, θ, H, MAR)
    transform_vertiecs(hcat(V...), Λ, S, θ, H, MAR)
end


function generate_vertex_style(VSTYLE, VS, VPH, default_VS, default_R, A, Ain, M, Cin, Ninner)
    grey(s) = ("gray",default_R*s[2],s[3])
    defV = ( "gray", default_VS, :circle )
    get_vstyle(a,m) = ( (a∈Ain) ? get(VSTYLE,m,defV) : grey(get(VSTYLE,m,defV)) )
    vertextuples = map( v->get_vstyle(v[1],v[2]), zip(A,M) )
    vertexstyle  = map( m->m[3], vertextuples )

    if VS===nothing
        vertexrad   = map( m->m[2], vertextuples )
    else
        #NOTE this branch is used by show_ham_eigenmode_svg()
        #NOTE which sets vertex_sizes to be amplitudes of eigenmode
        @assert length(VS)==Ninner #num_inner_sites(latt)
        vertexrad = map( a->((a∈Ain) ? VS[Cin[a]] : 0.0), A )
    end

    if VPH===nothing
        vertexcolor = map( m->m[1], vertextuples )
    else
        #NOTE this branch is used by show_ham_eigenmode_svg()
        #NOTE which sets vertex_sizes to be amplitudes of eigenmode
        @assert length(VPH)==Ninner #num_inner_sites(latt)
        @inline HSV(r) = colorschemes[:hsv][1+Int(floor(100*(r%1)))]
        col(φ) = "#"*hex(HSV((mod(φ+8π,2π)/2π)))
        vertexcolor = map( a->((a∈Ain) ? col(VPH[Cin[a]]) : "gray"), A )
    end
    return vertexcolor, vertexrad, vertexstyle
end

function generate_edge_style(ESTYLE, default_ES, default_R, A, Ain, E, F)
    grey(s) = ("gray",default_R*s[2],s[3])
    defE = ( "gray", default_ES, :solid  )
    get_estyle(e,f) = ( (A[e[1]]∈Ain && A[e[2]]∈Ain) ? get(ESTYLE,f,defE) : grey(get(ESTYLE,f,defE)) )
    edgetuples = map( v->get_estyle(v[1],v[2]), zip(E,F) )
    edgecolor = map( v->v[1], edgetuples )
    edgewidth = map( v->v[2], edgetuples )
    edgestyle = map( v->v[3], edgetuples )
    return edgecolor, edgewidth, edgestyle
end

# ------------------------------------------------
#NOTE
# bounds6 : to restrict the lattice site
# miller_index overwrites view_axes
# ------------------------------------------------
function show_lattice_svg(
    latt0::Lattice{T1,TR},
    VSTYLE::Dict,
    ESTYLE::Dict;
    ORIGIN       = nothing,
    YMAX         = -1,
    bounds6      = (-1e5,1e5,-1e5,1e5,-1e5,1e5),
    rotdeg       = 0,
    view_axes    = [1,2],
    miller_index = nothing,
    upscale      = 20.0,
    margin       = 20.0,
    default_R    = 0.618,
    default_ES   = 1.0,
    default_VS   = 2.0,
    VS           = nothing, #NOTE used by show_ham_eigenmode_svg()
    VPH          = nothing, #NOTE used by show_ham_eigenmode_svg()
    highlight    = []
    ) where { T1<:AbstractMatrix,TR<:Real }
    
    #
    latt = convert(Lattice{T1,Float64}, latt0)

    #
    (V0,M,E,F,A,C,Ain,Cin) = ( (miller_index===nothing)
                                ? lattice_to_graph(latt,view_axes,bounds6)
                                : lattice_to_graph(latt,MillerIndex(latt.UC.dim,miller_index),bounds6) )

    # add margin and invert y-axis ( the canvas of svg has downward y-axis )
    (V, S, Height) = transform_vertiecs(V0, upscale, ORIGIN, rotdeg, YMAX, margin)

    # highlights
    HI = [ (transform_vertiecs(point_list,upscale,S,rotdeg,Height,margin)|>first,col)
            for (point_list,col) ∈ highlight ]
    HI = vcat( HI, [(transform_vertiecs([(0 .* latt.R0[view_axes,1]),],upscale,S,rotdeg,Height,margin)|>first,"yellow"),] )

    # edge styles
    (edgecolor, edgewidth, edgestyle) = generate_edge_style(
                                                ESTYLE,
                                                default_ES, default_R,
                                                A, Ain, E, F )

    # vertices styles
    (vertexcolor, vertexrad, vertexstyle) = generate_vertex_style(
                                                VSTYLE, VS, VPH,
                                                default_VS, default_R,
                                                A, Ain, M, Cin,
                                                num_inner_sites(latt) )

    return plot_graph(  V,  E,
                        MARGIN = margin,
                        VERTEXCOLOR=vertexcolor, VERTEXRAD=vertexrad, VERTEXSTYLE=vertexstyle,
                        EDGECOLOR=edgecolor, EDGEWIDTH=edgewidth, EDGESTYLE=edgestyle,
                        HIGHLIGHT= HI  )
end

function generate_edge_weight(A,Ain,Cin,E,D)
    # C = Dict(k=>i for (i,k) ∈ enumerate(A)) # position index
    # E = [(C[i],C[j]) for (i,j)   ∈ zip(I,J)   if (i ∈ A && j ∈ A)]
    # F = [f           for (i,j,f) ∈ zip(I,J,F) if (i ∈ A && j ∈ A)]
    return [ ((A[i]∈Ain && A[j]∈Ain) ? abs(D[Cin[A[i]],Cin[A[j]]]) : 0.0) for (i,j) ∈ E ]
end

function generate_edge_style_by_weight(ESTYLE, D, default_ES, default_R, A, Ain, E, F, Cin)
    W = generate_edge_weight(A,Ain,Cin,E,D)
    grey(s) = ("gray",default_R*s[2],s[3])
    defE    = ("gray", default_ES, :solid)
    get_estyle(e,f) = ( (A[e[1]]∈Ain && A[e[2]]∈Ain) ? get(ESTYLE,f,defE) : grey(get(ESTYLE,f,defE)) )
    edgetuples = map( v->(get_estyle(v[1],v[2]),v[3]), zip(E,F,W) )
    edgecolor  = map( v->v[1][1],      edgetuples )
    edgewidth  = map( v->v[1][2]*v[2], edgetuples )
    edgestyle  = map( v->v[1][3],      edgetuples )
    return edgecolor, edgewidth, edgestyle
end

function show_lattice_edge_weight_svg(
    latt0::Lattice{T1,TR},
    VSTYLE::Dict,
    ESTYLE::Dict,
    D::T;
    ORIGIN       = nothing,
    YMAX         = -1,
    bounds6      = (-1e5,1e5,-1e5,1e5,-1e5,1e5),
    rotdeg       = 0,
    view_axes    = [1,2],
    miller_index = nothing,
    upscale      = 20.0,
    margin       = 20.0,
    default_R    = 0.618,
    default_ES   = 1.0,
    default_VS   = 2.0,
    VS           = nothing, #NOTE used by show_ham_eigenmode_svg()
    VPH          = nothing, #NOTE used by show_ham_eigenmode_svg()
    highlight    = []
    ) where { T1<:AbstractMatrix, TR<:Real, T<:AbstractMatrix }

    #
    latt = convert(Lattice{T1,Float64}, latt0)

    # E =
    (V0,M,E,F,A,C,Ain,Cin) = ( (miller_index===nothing)
                                ? lattice_to_graph(latt,view_axes,bounds6)
                                : lattice_to_graph(latt,MillerIndex(latt.UC.dim,miller_index),bounds6) )

    # add margin and invert y-axis ( the canvas of svg has downward y-axis )
    (V, S, Height) = transform_vertiecs(V0, upscale, ORIGIN, rotdeg, YMAX, margin)

    # highlights
    HI = [ (transform_vertiecs(point_list,upscale,S,rotdeg,Height,margin)|>first,col)
            for (point_list,col) ∈ highlight ]
    HI = vcat( HI, [(transform_vertiecs([(0 .* latt.R0[view_axes,1]),],upscale,S,rotdeg,Height,margin)|>first,"yellow"),] )

    # edge styles
    (edgecolor, edgewidth, edgestyle) = generate_edge_style_by_weight(
                                                ESTYLE, D,
                                                default_ES, default_R,
                                                A, Ain, E, F, Cin)

    # vertices styles
    (vertexcolor, vertexrad, vertexstyle) = generate_vertex_style(
                                                VSTYLE, VS, VPH,
                                                default_VS, default_R,
                                                A, Ain, M, Cin,
                                                num_inner_sites(latt) )

    return plot_graph(  V,  E,
                        MARGIN = margin,
                        VERTEXCOLOR=vertexcolor, VERTEXRAD=vertexrad, VERTEXSTYLE=vertexstyle,
                        EDGECOLOR=edgecolor, EDGEWIDTH=edgewidth, EDGESTYLE=edgestyle,
                        HIGHLIGHT= HI  )
end


#=
function show_unitcell_svg(
    UC::UnitCell,
    neighbor_num_max;
    view_axes = [1,2],
    upscale=20.0,
    vertex_styles = [("black",3.0,:dot),],
    edge_styles   = [("blue",4.0,:solid),("red",3.2,:solid),("green",2.4,:solid),("gray",1.6,:solid)],
    highlight = [],
    use_origin=false, ORIGIN=nothing, VY_MINMAX=0, rotdeg=0, #XXX hot-fix VY_MINMAX
    )
    @inline symext(s,i,j) = Symbol(string(s)*string(i)*string(j))
    Nsubl = UC.nsubl
    LD1 = [ (i,j)=>(symext(:f,i,j),neighbor_num_max)
                for i=1:Nsubl for j=1:Nsubl if i>j ]
    LD2 = [ (j,i)=>(symext(:f,i,j),neighbor_num_max)
                for i=1:Nsubl for j=1:Nsubl if i>j ]
    LD3 = [ (i,i)=>(symext(:g,i,i),neighbor_num_max)
                for i=1:Nsubl ]
    LD = Dict(vcat(LD1,LD2,LD3)...)
    BBOX = BoundingBox( (-0.01.+0.0.*UC.a[:,1],
                         [[1,0,0] [0,1,0] [0,0,1]],
                         [1,1,1],
                         [true,true,true]) )
    LN = link_info_sublattice_pairs_by_nth(LD, UC)
    LATT = build_lattice(UC, LN, BBOX)
    show_lattice_svg( LATT,
                      use_origin=use_origin, ORIGIN=ORIGIN, VY_MINMAX=VY_MINMAX,
                      rotation_degrees=rotation_degrees, #XXX hot-fix VY_MINMAX
                      view_axes=view_axes, upscale=upscale,
                      vertex_styles=vertex_styles, edge_styles=edge_styles,
                      highlight=highlight )
end
=#

function DR_mode_to_ΔR(
    DR_mode::Vector{Float64},
    latt::Lattice,
    MASS::Dict
    )::Array{Float64,2}

    dim = dimensions(latt)
    @assert length(vec(DR_mode)) == dim*num_inner_sites(latt)
    inner_site_ids = inner_site_id(latt)
    outer_border_site_ids = outer_border_site_id(latt)
    C = Dict( n=>i for (i,n) ∈ enumerate(inner_site_ids) )
    inner_equivalents = [C[k] for k ∈ inner_equivalent_of_outer_border_site_id(latt)]

    Msqrtinv = map( m->1.0/sqrt(get(MASS,m,1.0)), latt.UC.m[latt.SL[inner_site_ids]] )
    DR_mode_R0 = reshape( real.(vec(DR_mode)), dim, : ) .* (Msqrtinv')

    DR_mode_ext = zeros(Float64, size(latt.R0))
    DR_mode_ext[:, inner_site_ids       ] = DR_mode_R0[:,:]
    DR_mode_ext[:, outer_border_site_ids] = DR_mode_R0[:,inner_equivalents]

    return DR_mode_ext
end


function DQ_mode_to_ΔR(
    kvec::Vector{Float64},
    DQ_mode::Vector{ComplexF64},
    M::Dict,
    latt::Lattice;
    φ0 = 0.0
    )::Matrix{Float64}
    #%   latt1.R0[:,:] += ΔR
    R0 = zeros(Float64, size(latt.R0))
    inner_site_ids = inner_site_id(latt)
    DQ_R0   = reshape( vec(DQ_mode), dimensions(latt), : )
    ΔR = DQ_R0[:,latt.SL[inner_site_ids]] .* cis.((kvec'*latt.R0[:,inner_site_ids]).+φ0)
    mass = [1.0/sqrt(M[mi]) for mi ∈ latt.UC.m[latt.SL[inner_site_ids]]]
    R0[:,inner_site_ids] .= mass.*real.(ΔR)
    return R0
end


function show_diff_lattice_svg(
    ΔR_list::Vector{Array{Float64,2}},
    latt0::Lattice{T1,TR},
    VSTYLE::Dict,
    ESTYLE::Dict;
    view_axes    = [1,2],
    miller_index = nothing,
    bounds6      = (-1e5,1e5,-1e5,1e5,-1e5,1e5),
    upscale      = 20.0,
    margin       = 20.0,
    default_R    = 0.618,
    default_ES   = 1.0,
    default_VS   = 2.0,
    rotation_degrees = 0
    ) where { T1<:AbstractMatrix, TR<:Real }

    # plot lattice twice
    latt = convert(Lattice{T1,Float64}, latt0) # latt.R0 could be Array{Int64,2}
    latt1 = copy(latt)
    R0 = latt1.R0
    
    V0 = ( (miller_index===nothing)
            ? lattice_to_graph_V(latt,view_axes,bounds6)
            : lattice_to_graph_V(latt,MillerIndex(latt.UC.dim,miller_index),bounds6) )
    # add margin and invert y-axis ( the canvas of svg has downward y-axis )
    (V, ORIGIN, YMAX) = transform_vertiecs(V0, upscale, nothing, rotation_degrees, -1, margin)

    # gray layer for original lattice
    VSTYLE_gray = Dict( k=>("gray",     v[2], v[3]) for (k,v) ∈ VSTYLE )
    ESTYLE_gray = Dict( k=>("gray", 0.3*v[2], v[3]) for (k,v) ∈ ESTYLE )

    # latt0
    G0 = show_lattice_svg(  latt, VSTYLE_gray, ESTYLE_gray;
                            ORIGIN    = ORIGIN,
                            YMAX      = YMAX,
                            bounds6   = bounds6,
                            rotdeg    = rotation_degrees,
                            view_axes = view_axes,
                            miller_index = miller_index,
                            upscale   = upscale,
                            margin    = margin,
                            default_R = default_R,
                            default_ES = default_ES,
                            default_VS = default_VS,
                            highlight = []  )

    G = G0
    for ΔR ∈ ΔR_list
        # colored layer for lattice with mode on it
        latt1.R0 = R0 .+ ΔR
        G1 = show_lattice_svg(  latt1, VSTYLE, ESTYLE;
                                ORIGIN    = ORIGIN,
                                YMAX      = YMAX,
                                bounds6   = bounds6,
                                rotdeg    = rotation_degrees,
                                view_axes = view_axes,
                                miller_index = miller_index,
                                upscale   = upscale,
                                margin    = margin,
                                default_R = default_R,
                                default_ES = default_ES,
                                default_VS = default_VS,
                                highlight = []  )
        G = join_svg(G,G1)
    end
    return G
end


function show_diff_lattice_svg(
    ΔR::Array{Float64,2},
    latt0::Lattice{T1,TR},
    VSTYLE::Dict,
    ESTYLE::Dict;
    view_axes    = [1,2],
    miller_index = nothing,
    bounds6      = (-1e5,1e5,-1e5,1e5,-1e5,1e5),
    upscale      = 20.0,
    margin       = 20.0,
    default_R    = 0.618,
    default_ES   = 1.0,
    default_VS   = 2.0,
    rotation_degrees = 0
    ) where { T1<:AbstractMatrix, TR<:Real }

    show_diff_lattice_svg(
        [ΔR,],
        latt0,
        VSTYLE,
        ESTYLE;
        view_axes    = view_axes,
        miller_index = miller_index,
        bounds6      = bounds6,
        upscale      = upscale,
        margin       = margin,
        default_R    = default_R,
        default_ES   = default_ES,
        default_VS   = default_VS,
        rotation_degrees = rotation_degrees
    )

end



function show_dynmat_eigenmode_svg(
    m0::Vector{Float64},
    latt0::Lattice{T1,TR},
    MASS::Dict,
    VSTYLE::Dict,
    ESTYLE::Dict;
    view_axes    = [1,2],
    miller_index = nothing,
    bounds6      = (-1e5,1e5,-1e5,1e5,-1e5,1e5),
    upscale      = 20.0,
    margin       = 20.0,
    default_R    = 0.618,
    default_ES   = 1.0,
    default_VS   = 2.0,
    mode_scale   = 0.1,
    rotation_degrees = 0
    ) where { T1<:AbstractMatrix, TR<:Real }

    ΔR = mode_scale .* DR_mode_to_ΔR(vec(m0),latt,MASS)
    show_diff_lattice_svg(
        ΔR,
        latt0,
        MASS,
        VSTYLE,
        ESTYLE;
        view_axes    = view_axes,
        miller_index = miller_index,
        bounds6      = bounds6,
        upscale      = upscale,
        margin       = margin,
        default_R    = default_R,
        default_ES   = default_ES,
        default_VS   = default_VS,
        rotation_degrees = rotation_degrees
        )
end
