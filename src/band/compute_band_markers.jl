@inline transform_mat(kpoint,mat,OPS) = Float64[op(kpoint,mat[bandi,:]) for bandi=1:size(mat,1), op∈OPS]

compute_band_markers(BS0::BandStructure, OPS::Vector; TM=transform_mat) = BandStructure{Float64}(
    BS0.Kpath,
    BS0.Bands,
    # Vector{Pair{String, VectorKpointMatMTuple{T}}}
    [ s => KpointMatMTuple{Float64}[ (kpoint,TM(kpoint,mat,OPS)) 
                                        for (kpoint,mat) ∈ KMT ]
        for  (s,KMT) ∈ BS0.Markers ]
)

@inline flatten(A) = A[:]
compute_band_markers_and_collapse(BS0::BandStructure, OPS::Vector; TM=transform_mat) = hcat([ 
    [hcat(first.(KMT)...);  hcat(map(x->flatten(TM(x[1],x[2],OPS)), KMT)...)] 
    for (_,KMT) in BS0.Markers ]...)

CBMC(BS0::BandStructure, OPS::Vector; TM=transform_mat) = compute_band_markers_and_collapse(BS0, OPS; TM=TM)

compute_band_energy_markers_and_collapse(BS0::BandStructure, OPS::Vector; TM=transform_mat) = hcat([ 
    [hcat(first.(KMT)...);  hcat(last.(BND)...);  hcat(map(x->flatten(TM(x[1],x[2],OPS)), KMT)...)] 
    for ((_,BND),(_,KMT)) in zip(BS0.Bands,BS0.Markers) ]...)

CBEMC(BS0::BandStructure, OPS::Vector; TM=transform_mat) = compute_band_energy_markers_and_collapse(BS0, OPS; TM=TM)
