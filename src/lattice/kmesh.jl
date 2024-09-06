function b_krel(latt::Lattice)
    a = latt.UC.a
    (origin,trans,units,pbc) = latt.BBOX
    #NOTE the meaming of trans in BBOX is
    #(2) translation unit in terms of Bavais basis
    #    example t1=1*a1+2*a2, t2=2*a1+1*a2, trans = [[1,2] [2,1]]
    b1 = 2π.*inv(a*trans)'
    klists = [ ( t ? (collect(0:(units[i]-1)).//units[i]) : ([0//1,]) )
                 for (i,t) ∈ enumerate(pbc) ]
    # kx = 0 if pbc[1]==false
    return b1, klists
end


# TODO all dimensions
function kmesh(latt::Lattice)
    dim = dimensions(latt)
    (b1, klists) = b_krel(latt)
    mesh = zeros(eltype(b1),tuplejoin(dim,length.(klists)))
    if dim==2
        for ik ∈ 1:length(klists[1])
            for jk ∈ 1:length(klists[2])
                mesh[:,ik,jk] = b1 * [ klists[1][ik], klists[2][jk] ]
            end
        end
    elseif dim==3
        for ik ∈ 1:length(klists[1])
            for jk ∈ 1:length(klists[2])
                for lk ∈ 1:length(klists[3])
                    mesh[:,ik,jk,lk] = b1 * [ klists[1][ik], klists[2][jk], klists[3][lk] ]
                end
            end
        end
    end
    return mesh
end

@inline is_in_FBZ(v,b) = all([ abs(dot(normalize(b[:,i]),v))<=0.5*norm(b[:,i])
                               for i ∈ 1:size(b,2) ])


function b_krel_FBZ(latt::Lattice,N)
    @assert length(N)==dimensions(latt)
    b = reciprocal_basis(latt)
    klists = [ (t ? (collect(-(N[i]-1):(N[i]-1)).//N[i])
                 : ([0//1,]))
                   for (i,t) ∈ enumerate(latt.BBOX[4]) ]
    # kx = 0 if pbc[1]==false
    return b, klists
end

function kmeshFBZ(latt::Lattice; N=[])
    dim = dimensions(latt)
    (b, klists) = b_krel_FBZ(latt,N)
    mesh = []
    if dim==2
        for ik ∈ 1:length(klists[1])
            for jk ∈ 1:length(klists[2])
                v = b * [klists[1][ik], klists[2][jk]]
                if is_in_FBZ(v,b)
                    push!(mesh,v)
                end
            end
        end
    elseif dim==3
        for ik ∈ 1:length(klists[1])
            for jk ∈ 1:length(klists[2])
                for lk ∈ 1:length(klists[3])
                    v = b * [ klists[1][ik], klists[2][jk], klists[3][lk] ]
                    if is_in_FBZ(v,b)
                        push!(mesh,v)
                    end
                end
            end
        end
    end
    return mesh
end
