function diff_eig(
    en1::Vector, markers1::Vector,
    en2::Vector, markers2::Vector
    )
    @assert length(en1)==length(markers1)==length(en2)==length(markers2)
    deig = []
    for (e1,m1,e2,m2) ∈ zip(en1,markers1,en2,markers2)
        #@assert e1[1] ≈ m1[1] ≈ e2[1] ≈ m2[1]
        N = length(e1[2])
        @assert length(e2[2]) == N
        @assert size(m1[2]) == size(m2[2]) == (N,N)

        gram = abs.( m1[2] * m2[2]' )
        de = zeros(Float64,(N,2))
        for i=1:N
            (vmax,imax) = findmax(vec(gram[i,:]'))
            @assert abs(vmax) > (1.0/N)
            diff_en = real(e2[2][imax] - e1[2][i])
            de[i,1] = (diff_en>=0) ?  diff_en : 0
            de[i,2] = (diff_en <0) ? -diff_en : 0
        end
        push!(deig,(e1[1],de))
    end
    return deig
end


function diff_bands(B1::BandStructure, B2::BandStructure)
    @assert first.(B1.Kpath) == first.(B2.Kpath)
    @assert length(B1.Bands) == length(B2.Bands) == length(B1.Markers) == length(B2.Markers)
    Markers = [kv[1] => [] for kv ∈ B1.Markers]
    for i ∈ 1:length(B1.Bands)
        bands1 = B1.Bands[i]
        markers1 = B1.Markers[i]
        bands2 = B2.Bands[i]
        markers2 = B2.Markers[i]
        @assert bands1[1] == bands2[1] == markers1[1] == markers2[1]
        Markers[i] = (bands1[1] => diff_eig(bands1[2],markers1[2],bands2[2],markers2[2]))
    end
    db = copy(B1)
    db.Markers = Markers
    return db
end


@inline function diff_bands(B0::BandStructure, Bs...)
    join_markers(((diff_bands(B0,B) for B ∈ Bs))...)
end
