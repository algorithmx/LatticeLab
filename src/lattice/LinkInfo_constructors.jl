#: =====================================================

# used in link_info_by_***()
function generate_T(
    uc::UnitCell,
    bbox::BoundingBox,
    margin::Int,
    max_distance=1e10
    )::Array
    #* Q: how can I use it for 2d case ? 
    #* set bbox pbc = [true, true, false]
    N_min, N_max = bounding_box_Nmin_Nmax(uc.a, bbox, margin)
    Tvecs = [(uc.a*collect(t)) for t ∈ iter_t(N_min, N_max, uc.dim)]
    return hcat([v for v in Tvecs if norm(v)<max_distance]...)
end

#: =====================================================


#: =====================================================
#:  LinkInfo constructors

#>  Dict(:link_name => (d₁,[]),) 
#>  Dict(:link_name => (d₁,[n₁, n₂, ...]),) 
#>  =>  link_info_by_distance_direction

#>  Dict(:link_name => (distance,subl),) 
#>  (uc.m[j],uc.m[i])==subl
#>  =>  link_info_sublattice_pairs_by_distance

#>  ...
#>  =>  link_info_all_different

#>  Dict((s₁,s₂)=>(:link_name_prefix, max_nth_neighbor),)
#>  sᵢ sublattice id 
#>  =>  link_info_sublattice_pairs_by_nth

#>  Dict((el₁,el₂)=>(:link_name_prefix, max_nth_neighbor),)
#>  elᵢ element label ( ∈ uc.m ) 
#>  =>  link_info_element_pairs_by_nth

#+  UC, SUC, FC[subl_1, subl_2, i, j]
#+  Phonopy: p2s_map, p2p_map, s2p_map
#+  Dict()
#+  =>  link_info_phonopy_fc

#: =====================================================

global const _DEFAULT_BBOX_ = BoundingBox((
        [-0.001,-0.0004,-0.0006],
        [1 0 0; 0 1 0; 0 0 1],
        [1,1,1],
        [true,true,true]))

global const _DEFAULT_BBOX_2D_Z_ = BoundingBox((
        [-0.001,-0.0004,-0.0006],
        [1 0 0; 0 1 0; 0 0 1],
        [1,1,1],
        [true,true,false]))

#: =====================================================


function link_info_by_distance_direction(
    link_list::Dict,    #> Dict(:link_name => (d₁,[]),)
                        #> Dict(:link_name => (d₁,[n₁, n₂, ...]),)
    uc::UnitCell;
    bounding_box=_DEFAULT_BBOX_,
    rounding_digits=15
    )::LinkInfo

    EPS = 1e5*f1e(rounding_digits)
    @assert check_compat(uc)
    δ = uc.δ

    #% using a large `margin` doesn't cause problems
    #% the distance d₁ controls the span of the result
    #% the result does not contain R0-indices     #### c'est quoi ?!
    if uc.nsubl>  1
        shortest_δ_norm = minimum([norm(δ[:,i]-δ[:,j]) for i=1:uc.dim for j=i+1:uc.dim])
        large_enough_margin = Int(ceil(maximum(first.(values(link_list)))/shortest_δ_norm))+1
    else
        large_enough_margin = 5
    end
    dR = generate_T(uc, bounding_box, large_enough_margin)

    # helper functions
    @inline δδΔ(i,j,x) = (δ[:,j].+x).-δ[:,i] # r_i - r_j == x
    test1(i,j,Δ,R) = abs(norm(δδΔ(i,j,Δ))-R)<EPS
    @inline cosv(u,v)  = dot(u,v)/(norm(u)*norm(v))
    test12(i,j,Δ,R,d) = test1(i,j,Δ,R) && abs(cosv(δδΔ(i,j,Δ),d)-1)<EPS

    neighbors = Dict(k=>Spring() for k ∈ keys(link_list))
    for (k,sp) ∈ link_list
        (r0, directions) = sp
        if length(directions)==0
            #: case I example,  :link_name => [(d₁,[]), ...] with [] ~ ALL directions
            neighbors[k] = Spring(  (i,j) => [ δδΔ(i,j,dR[:,dr_i])
                                                    for dr_i ∈ 1:size(dR,2) 
                                                        if test1(i,j, dR[:,dr_i], r0) ]
                                    for i ∈ 1:uc.nsubl for j ∈ 1:uc.nsubl  ) |> delete_empty
        else
            #: case II example, :link_name => [(d₁,[n₁, n₂, ...]), ...]
            #:    with n₁ = [1,0,1], n₁ = [0,1,-1]
            neighbors[k] = Spring(  (i,j) => [ δδΔ(i,j,dR[:,dr_i])  
                                                    for dr_i ∈ 1:size(dR,2) 
                                                        for dir ∈ directions 
                                                            if test12(i,j, dR[:,dr_i], r0, dir) ]
                                    for i ∈ 1:uc.nsubl for j ∈ 1:uc.nsubl  ) |> delete_empty
        end
    end
    return LinkInfo( copy(uc), neighbors )
end


# modified from link_info_by_distance_direction()
# example for link_list:  Dict(:f1=>(BF(1),(:C,:P)),) ;
#* The link is directional
function link_info_sublattice_pairs_by_distance(
    link_list::Dict,    #> Dict(:link_name => (d₁,(s₁,s₂)),)
    uc::UnitCell;
    bounding_box = _DEFAULT_BBOX_,
    rounding_digits = 10
    )::LinkInfo

    EPS = f1e(rounding_digits)
    @assert check_compat(uc)
    δ = uc.δ

    # here, using a large `margin` doesn't cause problems
    # because it is the distance d₁ which controls the span of the result
    # the result does not contain R0-indices
    shortest_δ_norm = minimum([norm(δ[:,i]) for i=1:uc.dim])
    large_enough_margin = Int(ceil(maximum(first.(values(link_list)))/shortest_δ_norm))+1
    dR = generate_T(uc,bounding_box,large_enough_margin)
    cosv(u,v)  = dot(u,v)/(norm(u)*norm(v))
    δδΔ(i,j,x) = (δ[:,j].+x).-δ[:,i]
    test1(i,j,Δ,R) = abs(norm(δδΔ(i,j,Δ))-R)<EPS

    SP = keys(link_list) # symbols for links
    neighbors = Dict(k=>Spring((i,j)=>[] for i ∈ 1:uc.nsubl for j ∈ 1:uc.nsubl) for k ∈ SP)
    for k ∈ SP # match the sublattice labels and update neighbors[k]
        link_list_k = ( (link_list[k] isa Tuple) && length(link_list[k])==2 
                        ? ((link_list[k][2] isa Tuple) && length(link_list[k][2])==2
                            ? [link_list[k],]                                         # (BF(1),(:C,:P))
                            : ( eltype(link_list[k][2])<:Tuple                        # (BF(1),[(:C,:P), (:Cs,:P)])
                                ? Tuple[(link_list[k][1],t) for t in link_list[k][2]] # convert to [(BF(1),(:C,:P)), (BF(1),(:Cs,:P))]
                                : Tuple[ ] ))                                         # unknwon format  
                        : link_list[k])                                               # [(BF(1),(:C,:P)), (BF(1),(:Cs,:P))]
        # example link_list[:f1] is (BF(1),(:C,:P))
        for (distance,subl) in link_list_k 
            if length(subl) == 2 # avoid wrong formats
                # match the sublattice labels and update neighbors[k]
                for i ∈ 1:uc.nsubl
                    for j ∈ 1:uc.nsubl
                        if (uc.m[j],uc.m[i])==subl  
                            # this means that the pair (:C,:P) instructs  a link :C-->:P 
                            # since δδΔ(i,j,dr) dictates j --> i
                            neighbors[k][(i,j)] = vcat( neighbors[k][(i,j)],
                                                        [δδΔ(i,j,dR[:,dr_i]) for dr_i ∈ 1:size(dR,2) 
                                                         if test1(i,j,dR[:,dr_i],distance)] )
                        end
                    end
                end
            else
                @error "link_info_sublattice_pairs_by_distance() got dictionary of wrong format ($distance,$subl)"
            end
        end # (distance,subl) in link_list_k 
        neighbors[k] = delete_empty(neighbors[k])
    end
    return LinkInfo( copy(uc), neighbors )
end


#: =====================================================

function link_info_sublattice_pairs_by_nth(
    neighbor_num_for_each_nb_pair::Dict,    #> # Dict((1,2)=>(:f,2),(1,1)=>(:g,1),(2,2)=>(:g,1))
    uc::UnitCell;
    bounding_box = _DEFAULT_BBOX_,
    rounding_digits = 14
    )::LinkInfo

    EPS = 1e4*f1e(rounding_digits)
    P = neighbor_num_for_each_nb_pair
    δ = uc.δ
    large_enough_margin = maximum(last.(values(neighbor_num_for_each_nb_pair)))+1
    dR = generate_T(uc, bounding_box, large_enough_margin)

    # complete the neighbor_num_for_each_nb_pair::Dict
    # (i,j) implies (j,i)
    neighbor_keys = sck(P)
    for k ∈ neighbor_keys
        if (k[2],k[1]) ∉ neighbor_keys
            #NOTE CONVENTION in this way, (1,2) === (2,1)
            P[(k[2],k[1])] = neighbor_num_for_each_nb_pair[k]
        end
    end

    #
    neighbors = Dict{Symbol,Spring}( symbol_suffix(symb_nth[1],i)=>Spring()
                                     for (subl_pair,symb_nth) ∈ P for i ∈ 1:symb_nth[2] )
    for (subl_pair,symbol_nth) ∈ P
        (sl1,sl2) = subl_pair
        (symbol,nth) = symbol_nth
        # we don't distinguish directions here
        directions = sort([(δ[:,sl2].+dR[:,dr_i]).-δ[:,sl1] for dr_i ∈ 1:size(dR,2)], by=norm)
        distances = unique(x->round(x,digits=rounding_digits) , map(x->norm(x), directions))
        if sl1==sl2
            # remove distances[1] because it is zero
            @assert distances[1]<EPS
            distances = copy(distances[2:end])
        end
        # put numbers at the end of the spring constant symbol
        # to distinguish distances / n'th nearest neighbor
        symbols = [symbol_suffix(symbol,i) for i ∈ 1:nth]
        for i ∈ 1:nth
            symbol_i = symbols[i]
            dist_i = distances[i]
            l = sort( [dir for dir ∈ directions if abs(norm(dir)-dist_i)<EPS],
                      by=(x->mod(atan(x[2],x[1])+2π+1e-16,2π)) )
            if subl_pair ∈ keys(neighbors[symbol_i])
                throw("link_info_sublattice_pairs_by_nth()")
                push!(neighbors[symbol_i][subl_pair], l)
            else
                neighbors[symbol_i][subl_pair] = l
            end
        end
    end

    return LinkInfo( copy(uc), neighbors )
end


#NOTE example for neighbor_num_for_each_nb_pair
# Dict((1,2)=>(:f,2),(1,1)=>(:g,1),(2,2)=>(:g,1))
function link_info_sublattice_pairs_by_nth(
    Nth::Int64,
    uc::UnitCell;
    bounding_box = _DEFAULT_BBOX_,
    rounding_digits = 14
    )::LinkInfo

    @inline SJ(s1,s2) = Symbol(s1)>Symbol(s2) ? symbol_join(s1,s2) : symbol_join(s2,s1)
    
    LD = Dict( ((i==j) ? (i,i)=>(SJ(uc.m[i],uc.m[i]), Nth) : (i,j)=>(SJ(uc.m[i],uc.m[j]), Nth))
                for i ∈ 1:uc.nsubl for j ∈ 1:i )
    
    return link_info_sublattice_pairs_by_nth( LD, uc; 
                                              bounding_box=bounding_box, 
                                              rounding_digits=rounding_digits )
end

 
#: =====================================================

# one sublattice pair, one direction (DISTINGUISH_IJ), one link
#TODO rewrite
function link_info_all_different(
    DISTINGUISH_IJ::Bool,
    subl_pair_cutoffs::Dict,
    uc::UnitCell;
    bounding_box = _DEFAULT_BBOX_,
    rounding_digits = 14
    )::LinkInfo
    EPS = 1e4*f1e(rounding_digits)
    @assert check_compat(uc)
    δ = uc.δ
    dR = generate_T(uc,bounding_box,margin)
    δδΔ(i,j,Δ)        = (δ[:,j].+Δ).-(δ[:,i])
    test1(i,j,Δ,R)    = abs(norm(δδΔ(i,j,Δ))-R)<EPS
    dists(i,j)        = sort(uniquen([ norm(δδΔ(i,j,dR[:,dr_i])) for dr_i ∈ 1:size(dR,2) ], rounding_digits))
    directions(i,j,R) = [ δδΔ(i,j,dR[:,dr_i]) for dr_i ∈ 1:size(dR,2) if test1(i,j,dR[:,dr_i],R) ]
    neighbors = []
    #! CONVENTION FOR LINK DIRECTIONS : sli --> slj
    for sli ∈ 1:num_sublattice(uc)
        for slj ∈ 1:num_sublattice(uc)
            if (!DISTINGUISH_IJ) && slj>sli
                #* in this case, (sli,slj) is equivalent to (slj,sli)
                #* i.e. (sli,slj) in the lattice implies (slj,sli)
                #* and they have identical labels
                continue
            else
                #* for two cases:
                #* (1) (sli,slj) <=!!=> (slj,sli) 
                #* (2) (sli,slj) <==> (slj,sli) and sli<=slj
                link_symbol0 = symbol_join( uc.m[sli], uc.m[slj] )
                distij = dists(sli, slj)
                cutoff = get( subl_pair_cutoffs, (sli,slj), 1 )
                if sli == slj
                    @assert distij[1]<EPS
                    @assert length(distij) >= cutoff+1
                    distij = distij[2:end]
                else
                    @assert length(distij) >= cutoff
                end
                for l ∈ 1:cutoff
                    # append the l'th neighbor label l to the link symbol
                    link_symbol = symbol_join(link_symbol0, l)
                    # find distinguished directions from subl i to subl j
                    dirij = directions(sli, slj, distij[l])
                    dirs = sort( uniquen(dirij, rounding_digits), by=(x->mod(atan(x[2],x[1])+2π+1e-16,2π)) )
                    # for each direction, construct link info
                    dir_for_ii = []
                    for (k,dir) ∈ enumerate(dirs)
                        link_symbol_k = symbol_join(link_symbol,k)
                        if DISTINGUISH_IJ
                            #* one sublattice pair, one direction, one link
                            push!( neighbors, link_symbol_k=>Spring((sli,slj)=>[dir,]) )
                        else
                            if sli==slj
                                #* one sublattice pair, one direction, one link, BUT
                                #* for sublattice (i,i), dir and -dir is indistinguishable 
                                #* (consider 2nd nearest nb hopping in honeycomb lattice)
                                #* so we need to avoid double-counting
                                if all([norm(-d.-dir)>EPS for d ∈ dir_for_ii])
                                    push!( neighbors, link_symbol_k=>Spring((sli,slj)=>[dir,-dir]) )
                                    push!( dir_for_ii, dir )
                                end
                            else
                                push!( neighbors, link_symbol_k=>Spring((sli,slj)=>[dir,], (slj,sli)=>[-dir,]) )
                            end
                        end
                    end
                end # for l ∈ 1:cutoff
            end # (!DISTINGUISH_IJ) && slj>sli
        end # slj ∈ 1:num_sublattice(uc)
    end # sli ∈ 1:num_sublattice(uc)
    return LinkInfo( copy(uc), Dict(neighbors) )
end


# one sublattice pair, one direction (DISTINGUISH_IJ), one link
#TODO rewrite
function link_info_all_different(
    DISTINGUISH_IJ::Bool,
    all_cutoffs::Int64,
    uc::UnitCell;
    bounding_box = _DEFAULT_BBOX_,
    rounding_digits = 14
    )::LinkInfo
    subl_pair_cutoffs = Dict( (i,j)=>all_cutoffs for i ∈ 1:uc.nsubl for j ∈ 1:uc.nsubl )
    return link_info_all_different( DISTINGUISH_IJ, subl_pair_cutoffs, uc;
                                    bounding_box=bounding_box, 
                                    rounding_digits=rounding_digits )
end


## --------------------------------------------------------------


search_interval(sorted_list,L,R) = searchsortedfirst(sorted_list,L):searchsortedlast(sorted_list,R)

is_continuous(V) = all(last.(V[1:end-1]).+1 .== first.(V[2:end]))


function compute_distance_group(
    dRij::Array, 
    maxlen::Int, 
    rdg::Int, 
    EPS::Float64, 
    cond::Bool
    )
    N  = [norm(dRij[:,s]) for s ∈ 1:size(dRij,2)]
    Np = sortperm(N)
    Ns = N[Np]
    Na = sort(uniquen(N,rdg))
    if  cond
        if !(Na[1]<EPS)  throw(DistanceListHeadError(""))  end
        if !(length(Na)>=maxlen+1)  throw(DistanceListTooShortError("length(Na)=$(length(Na))"))  end
        Na = Na[2:end]
    else
        if !(length(Na)>=maxlen)  throw(DistanceListTooShortError("length(Na)=$(length(Na))"))  end
    end
    _ran = [search_interval(Ns, R-(1-1e-15)*EPS, R+(1-1e-15)*EPS) for R ∈ Na[1:maxlen]]
    if !is_continuous(_ran)
        throw(DiscontinuousDistanceGroupError("ranges = $_ran"))
    end
    return [Np[p] for p in _ran]
end



"""

    function all_links_between_sublattices(
        uc::UnitCell;             # unit cell
        dR_cutoff = 8,            # determined by supercell size
        distance_maxlen = 30,     # should be smaller than the length of generate_dR result
        rounding_digits = 10
    )::Dict

# Note :

CONVENTION FOR LINK DIRECTIONS : sli --> slj

It does not contradict to the convention that  f[b,k] = sp from k to b .

It returns a dictionary : 

Dict( (i,j)=>[(:si_i_sj_j_distlevel_samedistk,[x,y,z]),], )

TODO : test 2D case

"""
function all_links_between_sublattices(
    uc::UnitCell;
    bounding_box = _DEFAULT_BBOX_,
    maxlen=30,
    max_distance=15.0,
    large_enough_margin=8,
    rounding_digits=4
    )::Dict{Tuple{Int64,Int64}, Vector{Tuple{Symbol,Vector}}}

    ⊛(a,b) = symbol_join(a,b)
    if !check_compat(uc)  throw(IncorrectUnitCellError("$uc"))  end
    dR = generate_T(uc, bounding_box, large_enough_margin, max_distance)
    Nsubl = num_sublattice(uc)

    #*** PRECISION CONTROL ***
    neighbors = []
    for sli ∈ 1:Nsubl
        for slj ∈ 1:Nsubl
            dRij = (uc.δ[:,slj].-uc.δ[:,sli]).+dR  #: CONVENTION FOR LINK DIRECTIONS : sli --> slj
            dist_group = compute_distance_group( dRij, maxlen, 
                                                 rounding_digits, (1//2)*f1e(rounding_digits), 
                                                 sli==slj )
            link_symbol0 = (uc.m[sli] ⊛ sli) ⊛ (uc.m[slj] ⊛ slj) 
            #: LINK DIRECTIONS : sli --> slj
            push!(  neighbors, 
                    (sli,slj) => [ (link_symbol0 ⊛ l ⊛ k, dRij[:,idx])
                                        for l = 1:maxlen
                                            for (k,idx) ∈ enumerate(dist_group[l]) ] )
        end # slj
    end # sli
    return Dict(neighbors) |> delete_empty
end

#% format of V
#* [(:C3_3_C3_3_1_1, [-3.59591717, 0.0, 0.0]), (:C3_3_C3_3_1_2, [3.59591717, 0.0, 0.0]), 
#*  (:C3_3_C3_3_2_1, [0.0, -4.68684147, 0.0]), ... ]

pick_nearest(v, nth::Int, d=6) =
    v[sortperm(v, by=x->round(norm(x[2]),digits=d))[1:nth]]

pick_nearest_by_ij(all_ln::Dict, nn_dict::Dict{Tuple{Int,Int},Int}; digits=6) = 
    Dict( k=>pick_nearest(v,get(nn_dict,k,0),digits) for (k,v) ∈ all_ln )

pick_nearest_by_ij(all_ln::Dict, nn::Int; digits=6) = 
    Dict( k=>pick_nearest(v,nn,digits) for (k,v) ∈ all_ln )

pick_dist(v, dist::Real, d=10) =
    [x for x in v if round(norm(x[2]),digits=d)<dist]

pick_dist_by_ij(all_ln::Dict, dd_dict::Dict{Tuple{Int,Int},R}; digits=10) where {R<:Real} = 
    Dict( k => pick_dist(v, get(dd_dict,k,0.0), digits) for (k,v) ∈ all_ln ) |> delete_empty

pick_dist_by_ij(all_ln::Dict, dd::Real; digits=10) = 
    Dict( k=>pick_dist(v,dd,digits) for (k,v) ∈ all_ln ) |> delete_empty


function all_links_dict_to_sp_no_directions(
    all_ln_dict::Dict
    )
    @inline SPX(sb) = sort(split(string(sb),"_",keepempty=false)[[1,3,5]])
    @inline JN(l) = Symbol(join(l[[2,3,1]],"_"))
    @inline elem01(x) = length(x)==0 ? [] : first.(x)
    labels = vcat(elem01.(values(all_ln_dict))...) .|> SPX |> unique .|> JN
    ret = Dict(k=>Pair[] for k in labels)
    for (k,v) ∈ all_ln_dict
        for (s,d) ∈ v
            ss = JN(SPX(s))
            push!(ret[ss], k=>[d,])
        end
    end
    return Dict(k=>list_of_pairs_to_dict_merge_v(v) for (k,v) ∈ ret)
end


function all_links_dict_to_sp_no_directions_new_BUG(
    all_ln_dict::Dict
    )
    @inline take_first_value(x) = length(values(x))==0 ? [] : first.(values(x))
    @inline splt_ndrscr_135(sb) = sort(split(string(sb),"_",keepempty=false)[[1,3,5]])
    @inline make_new_symbol(l) = Symbol(join(l[[2,3,1]],"_"))
    @inline transform_symbol(sb) = make_new_symbol(splt_ndrscr_135(sb))
    labels = vcat(take_first_value(all_ln_dict)...) .|> transform_symbol |> unique
    ret = Dict(k=>Pair[] for k in labels)
    for (k,v) ∈ all_ln_dict
        for (s,d) ∈ v
            ss = transform_symbol(s)
            push!(ret[ss], k=>[d,])
        end
    end
    return Dict(k=>list_of_pairs_to_dict_merge_v(v) for (k,v) ∈ ret)
end
