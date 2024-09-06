
check_valid_EqV(eqv::Vector{Int64}) = all(eqv[unique(eqv[eqv.> 0])].==0)

"""
    position_index(eqv::Vector{Int64})

#NOTE

===

returns id, cc

`id`  is the inner site id 

`cc`  is the (equivalent) position of site `i` in `id`

"""
function position_index(eqv::Vector{Int64})
    @assert check_valid_EqV(eqv)
    # index of inner sites in R0
    id = inner_site_id_from_eqv(eqv) # inner_site_id(latt::Lattice)
    # position index for inner sites
    dic = Dict(id[i]=>i for i ∈ 1:length(id))
    # auxiliary function for index i of R0
    # for inner site return site position (in array id)
    # for outer (irrelevant) site return -1
    # for boundary sites return site position (in array id)
    # of the euivalent inner site
    f(i) = ((eqv[i]==0) ? dic[i] : (eqv[i]>0 ? dic[eqv[i]] : -1))
    # (index of R0) --> (index of inner sites)
    cc = map(f,1:length(eqv)) # "pull back" to inner
    @assert all(cc.<=length(id)) && all((cc.>0).|(cc.==-1)) # all sites in R0 are mapped
    @assert all([cc[id[k]]==k for k ∈ 1:length(id)])  # it is really the inverse of id
    return id, cc
end
