# design: independent of Hamiltonian / Dynamical matrix
mutable struct BandStructureGrid{T}

    Kgrid::Array{Float64,2}

    Bands::Array{T,2}

    Markers::Array{ComplexF64,2}

end

@inline check_compat(BS::BandStructureGrid) = true

copy(bs::BandStructureGrid{T}) where T = BandStructureGrid{T}( copy(bs.Kgrid), copy(bs.Bands), copy(bs.Markers) )

# =========================================================
