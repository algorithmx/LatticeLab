##: ============================================================

"""
    MinimalLattice{NL,TN,TE}

A `MinimalLattice` consists of the following minimal data of a lattice:

    + `V::Vector{Pair{VL,VF}}` : a list of vertex_label => vertex_feature pairs

    + `E::Vector{Pair{Tuple{VL,VL},EF}}` : a list of edge_label => edge_feature pairs

"""
mutable struct MinimalLattice{VL,VF,EF}
    V::Vector{Pair{VL,VF}}
    E::Vector{Pair{Tuple{VL,VL},EF}}
end


function MinimalLattice(vl::Vector, el::Vector)
    NL = typeof(vl[1][1])
    TN = typeof(vl[1][2])
    TE = typeof(el[1][2])
    @assert all([(e[1][1] isa NL) && (e[1][2] isa NL) for e in el])
    MinimalLattice{NL,TN,TE}(vl, el)
end


function make_featured_vertex_list(latt::Lattice)::Vector{Pair{Int,Tuple}}
    IN = inner_site_id(latt)
    r0 = latt.R0[:,IN[1]]
    pos(i) = latt.R0[:,i].-r0
    sl(i)  = latt.SL[i]
    ms(i)  = latt.UC.m[sl(i)]
    orb(i) = latt.UC.ξ[sl(i)]
    Pair{Int,Tuple}[(i=>(pos(p), ms(p), orb(p))) for (i,p) ∈ enumerate(IN)]
end


function make_featured_edge_list(latt::Lattice)::Vector{Pair{Tuple,Tuple}}
    IN, C = position_index(latt.EqV)
    Nsites = length(C)
    P,Q,V = findnz(latt.f)
    PQ = collect(zip(P,Q))
    sPQ = Set(PQ)
    sIN = Set(IN)
    linked(i,j) = ((i,j) ∈ sPQ)
    inside(i,j) = (i ∈ sIN) || (j ∈ sIN)
    direc(i,j)  = latt.R0[:,i].-latt.R0[:,j]
    lab(i,j)    = V[findfirst(x->x==(i,j),PQ)]
    Pair{Tuple,Tuple}[
        ((C[i], C[j]) => (direc(i, j), lab(i, j)))
        for i=1:Nsites for j=1:Nsites 
            if linked(i,j) && inside(i,j)] |> unique
end


@inline MinimalLattice(latt::Lattice) =
            MinimalLattice(make_featured_vertex_list(latt), make_featured_edge_list(latt))


@inline function is_valid(minilatt::MinimalLattice) 
    v   = first.(minilatt.V)
    v1s = first.(first.(minilatt.E))
    v2s =  last.(first.(minilatt.E))
    all([(p ∈ v) for p ∈ v1s]) && all([(p ∈ v) for p ∈ v2s])
end


##: ============================================================

#vl = [1=>:A, 2=>:B, 3=>:C]
#el = [(1,2)=>(:AB, 1.0)]
#latt = MinimalLattice(vl,el)
#is_valid(latt)