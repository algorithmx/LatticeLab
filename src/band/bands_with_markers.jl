
# requirement for OP
# op(k,vec) -> value
# dot(vec,ModeQ(k,UC)*vec)
#NOTE structure of the result: V[:,1] = eigenvalues,
# V[:,2] = op1(eigenvecs)
# V[:,3] = op2(eigenvecs)
function _eigen_with_markers_(kvec,M,OP)
    sol = eigen(M)
    #vectors = hcat([normalize(sol.vectors[:,i]) for i ∈ 1:length(sol.values)]...)
    ASS = sortperm(sol.values,by=real)
    V = []
    push!(V, sol.values[ASS])
    for op in OP
        #push!(V, [op(kvec,vectors[:,i]) for i in ASS])
        push!(V, [op(kvec,sol.vectors[:,i]) for i in ASS])
    end
    return hcat(V...)
end


function _eigen_with_markers_(kvec,M,S,OP)
    sol = eigen(M,S)
    vectors = hcat([normalize(sol.vectors[:,i]) for i ∈ 1:length(sol.values)]...)
    ASS = sortperm(sol.values,by=real)
    V = []
    push!(V, sol.values[ASS])
    for op in OP
        push!(V, [op(kvec,vectors[:,i]) for i in ASS])
    end
    return hcat(V...)
end

# ------------------------------------------------------

function band_structure_with_markers(
    kpath::Vector{Pair{String,Vector{T}}},
    MATQ,
    COMP::Function,
    TEST::Function,
    POST::Function,
    MARK_OPS::Vector;
    Δk = 0.02,
    test_eigen = true,
    eps = 1e-8
    )::BandStructure  where {T<:Number}
    comp(λ,a,b) = COMP(λ,a,b,MATQ,MARK_OPS)
    test(eig) = TEST(eig, eps)
    latt = (MATQ isa Tuple) ? MATQ[1].LATT : MATQ.LATT
    bands_markers = band_structure_base( latt, kpath, comp, test;
                                         Δk=Δk, test_eigen=test_eigen )
    # split the markers from the join result
    KPATH   = [ q => first.(EIG) for (q,EIG) ∈ bands_markers ]
    bands   = [ q => map(x->POST((x[1],vec(x[2][:,1]),)), EIG) for (q,EIG) ∈ bands_markers ]
    markers = [ q => map(x->(x[1],Matrix(x[2][:,2:end])), EIG) for (q,EIG) ∈ bands_markers ]
    # bost-process : taking real part
    return BandStructure( KPATH, bands, markers )
end


function band_structure_with_markers(
    kpath::Vector{Pair{String,Vector{T}}},
    DQ::Union{DynamicalMatrixQ,DynamicalMatricesQ,DynamicalMatrixKp},
    OP::Vector;
    Δk = 0.02,
    test_eigen = true,
    eps = 1e-8
    )::BandStructure  where {T<:Number}
    return band_structure_with_markers( kpath, DQ, EIGEN2M, all_non_negative_real, SQRT_EIGEN2, OP;
                                        Δk = Δk, test_eigen = test_eigen, eps = eps )
end


##: ---------------------------------------------------------------------


function band_structure_with_eigenvectors0(
    n::Int64,
    kpath::Vector{Pair{String,Vector{T}}},
    MATQ;
    Δk         = 1e-3,
    test_eigen = true,
    eps        = 1e-8
    )::BandStructure  where {T<:Number}
    OP = [(k,v)->vec(v)[l] for l=1:n]  #! k not used
    return band_structure_with_markers(kpath, MATQ, OP, Δk=Δk, test_eigen=test_eigen, eps=eps)
end


function band_structure_with_eigenvectors(
    kpath::Vector{Pair{String,Vector{T}}},
    DQ::Union{DynamicalMatrixQ,DynamicalMatricesQ};
    Δk         = 1e-3,
    test_eigen = true,
    eps        = 1e-8
    )::BandStructure  where {T<:Number}
    n = DQ.LATT.UC.dim * DQ.LATT.UC.nsubl
    return band_structure_with_eigenvectors0(n, kpath, DQ, Δk=Δk, test_eigen=test_eigen, eps=eps)
end


##: ------------------------------------------------------------
##: ------------------------------------------------------------


function band_structure_with_markers(
    kpath::Vector{Pair{String,Vector{T}}},
    HQ::HoppingHamUnion,
    OP::Vector;
    Δk = 0.02,
    test_eigen = true,
    eps = 1e-8
    )::BandStructure  where {T<:Number}
    return band_structure_with_markers( kpath, HQ, EIGENM, all_real, REAL_EIGEN, OP;
                                        Δk = Δk, test_eigen = test_eigen, eps = eps )
end


function band_structure_with_markers(
    kpath::Vector{Pair{String,Vector{T}}},
    HSQ::Tuple{T1,T2},
    OP::Vector;
    Δk = 0.02,
    test_eigen = true,
    eps = 1e-8
    )::BandStructure  where { T1<:HoppingHamUnion, T2<:HoppingHamUnion, T<:Number }
    return band_structure_with_markers( kpath, HSQ, EIGENM, all_real, REAL_EIGEN, OP;
                                        Δk = Δk, test_eigen = test_eigen, eps = eps )
end


function band_structure_with_markers(
    kpath::Vector{Pair{String,Vector{T}}},
    BdGQ::Union{BdGHamiltonianQ,BdGHamiltonianenQ},
    OP::Vector;
    Δk = 0.02,
    test_eigen = true,
    eps = 1e-8
    )::BandStructure  where {T<:Number}
    return band_structure_with_markers( kpath, BdGQ, EIGENMBdG, all_real, REAL_EIGEN, OP;
                                        Δk = Δk, test_eigen = test_eigen, eps = eps )
end


function band_structure_with_eigenvectors(
    kpath::Vector{Pair{String,Vector{T}}},
    HQ::Union{HoppingHamiltonianQ,HoppingHamiltonianenQ};
    Δk         = 1e-3,
    test_eigen = true,
    eps        = 1e-8
    )::BandStructure  where {T<:Number}
    n = sum(length.(HQ.LATT.UC.ξ))
    return band_structure_with_eigenvectors0(n, kpath, HQ, Δk=Δk, test_eigen=test_eigen, eps=eps)
end


function band_structure_with_eigenvectors(
    kpath::Vector{Pair{String,Vector{T}}},
    HSQ::Tuple{T1,T2};
    Δk         = 1e-3,
    test_eigen = true,
    eps        = 1e-8
    )::BandStructure  where { T1<:HoppingHamUnion, T2<:HoppingHamUnion, T<:Number }
    n = sum(length.(HSQ[1].LATT.UC.ξ))
    return band_structure_with_eigenvectors0(n, kpath, HSQ, Δk=Δk, test_eigen=test_eigen, eps=eps)
end


##: ------------------------------------------------------------
##: ------------------------------------------------------------


@inline zip_marker(M1,M2) = (M1[1]=>[(m1[1],hcat(m1[2],m2[2])) for (m1,m2) ∈ zip(M1[2],M2[2])])

function join_markers(Bs...)
    @assert length(Bs) >= 1
    if length(Bs)==1
        return Bs[1]
    else
        @assert all([length(Bs[i].Markers)==length(Bs[i].Bands)==length(Bs[1].Markers) for i=1:length(Bs)])
        Markers = [kv[1] => [] for kv ∈ Bs[1].Markers]
        for i ∈ 1:length(Bs[1].Bands)
            bands1 = Bs[1].Bands[i]
            markers1 = Bs[1].Markers[i]
            # accumulate markers in 2:n to markers in 1
            for j ∈ 2:length(Bs)
                bands2   = Bs[j].Bands[i]
                markers2 = Bs[j].Markers[i]
                @assert bands1[1] == bands2[1] == markers1[1] == markers2[1]
                markers1 = zip_marker(markers1,markers2)
            end
            Markers[i] = (bands1[1] => markers1[2])
        end
        db = copy(Bs[1])
        db.Markers = Markers
        return db
    end
end
