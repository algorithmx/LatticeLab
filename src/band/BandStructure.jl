# type of v = [(K1, [a1,a2,a3,...]), (K2, [b1,b2,b3,...]), ...]
KpointVecMTuple{M} = Tuple{Vector{Float64},Vector{M}}
VectorKpointVecMTuple{M} = Vector{KpointVecMTuple{M}}

# type of v = [(K1, A1[:,:]), (K2, A2[:,:]), ...]
KpointMatMTuple{M} = Tuple{Vector{Float64},Matrix{M}}
VectorKpointMatMTuple{M} = Vector{KpointMatMTuple{M}}

KpointVec = Vector{Vector{Float64}}

# design: independent of Hamiltonian / Dynamical matrix
mutable struct BandStructure{T}

    Kpath::Vector{Pair{String,KpointVec}}

    Bands::Vector{Pair{String, VectorKpointVecMTuple{Float64}}}

    Markers::Vector{Pair{String, VectorKpointMatMTuple{T}}}

end


@inline check_compat(BS::BandStructure) = ( BS.Markers==[]
                                            || (length(BS.Bands)==length(BS.Markers)
                                                && all([ first(b)==first(m) && length(last(b))==length(last(m))
                                                         for (b,m) ∈ zip(BS.Bands,BS.Markers) ])) )


@inline copy_k_v_pair(kv::Tuple) = (kv[1],copy(kv[2]))

@inline copy_bands(b::Pair) = string(b[1]) => copy_k_v_pair.(b[2])

@inline copy_markers(b::Pair) = string(b[1]) => copy_k_v_pair.(b[2])

@inline copy_kpath(b::Pair) = string(b[1]) => copy(b[2])

copy(bs::BandStructure{T}) where T = BandStructure{T}( copy_kpath.(bs.Kpath), copy_bands.(bs.Bands), copy_markers.(bs.Markers) )


# =========================================================


function extract_kpath_to_line(
    bands::Vector;
    scale=1.0
    )
    k_prev = last(bands[1])[1][1]
    accum = 0.0
    bands_k = Dict()
    for band ∈ bands
        x_axis = Float64[]
        for (kpoint,eig) ∈ last(band)
            accum += scale*norm(kpoint.-k_prev)
            push!(x_axis, accum)
            k_prev = kpoint
        end
        bands_k[first(band)] = x_axis
    end
    return bands_k
end


function extract_kpath_to_line(
    BS::BandStructure;
    scale=1.0
    )
    extract_kpath_to_line(BS.Bands; scale=scale)
end
