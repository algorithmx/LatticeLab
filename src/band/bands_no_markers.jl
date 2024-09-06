function _eigen_no_markers_(M)
    return sort(eigvals(M),by=real)
end


function _eigen_no_markers_(M,S)
    return sort(eigvals(M,S),by=real)
end


##% ---------------------------------------------------------------


# specially for dynamical matrix
function EIGEN2(
    λ,A,B,
    DQ::Union{DynamicalMatrixQ,DynamicalMatricesQ,DynamicalMatrixKp}
    )
    K = (1-λ).*A.+(λ).*B
    M = convert(Matrix{ComplexF64}, DynamicalMatrix(K,DQ))
    return ( real.(K) , _eigen_no_markers_(M) )
end


# specially for dynamical matrix
function EIGEN2M(
    λ,A,B,
    DQ::Union{DynamicalMatrixQ,DynamicalMatricesQ,DynamicalMatrixKp},
    OP::Vector
    )
    K = (1-λ).*A.+(λ).*B
    M = convert(Matrix{ComplexF64}, DynamicalMatrix(K,DQ))
    return ( real.(K) , _eigen_with_markers_(K,M,OP) )
end


# specially for dynamical matrix
@inline SQRT_EIGEN2(EIGEN2_result) = ( EIGEN2_result[1], map(x2->(2Int(real(x2)>1e-20)-1)*sqrt(abs(real(x2))),EIGEN2_result[2]) )


# specially for dynamical matrix
@inline function all_non_negative_real(list_EIGEN2_results, ϵ0)
    t1(t) = all(real.(t).>(-ϵ0)) && all(abs.(imag.(t)).<ϵ0)
    t2(t) = all(real.(t[:,1]).>(-ϵ0)) && all(abs.(imag.(t[:,1])).<ϵ0)
    test(t) = (ndims(t)==1) ? t1(t) : t2(t)
    return all( test(last(eig)) for eig ∈ list_EIGEN2_results )
end


##% ---------------------------------------------------------------


function band_structure(
    kpath::Vector{Pair{String,Vector{T}}},
    MATQ,
    COMP::Function,
    TEST::Function,
    POST::Function;
    Δk = 0.02,
    test_eigen = true,
    eps = 1e-8
    ) where {T<:Number}
    comp(λ,a,b) = COMP(λ,a,b,MATQ)
    test(eig)   = TEST(eig, eps)
    latt  = (MATQ isa Tuple) ? MATQ[1].LATT : MATQ.LATT
    bands = band_structure_base( latt, kpath, comp, test;
                                 Δk=Δk, test_eigen=test_eigen )
    KPATH = [q=>first.(EIG) for (q,EIG) in bands]
    # post-process : taking real part
    return  BandStructure{Float64}(
                KPATH, 
                map(x->(x[1]=>POST.(x[2])),bands), 
                Pair{String, VectorKpointMatMTuple{Float64}}[] )
end


function band_structure(
    kgrid::Array{Float64,2},
    MATQ,
    COMP::Function,
    TEST::Function,
    POST::Function;
    test_eigen = true,
    eps = 1e-8
    )
    comp(q)   = COMP(0,q,q,MATQ)
    test(eig) = TEST(eig, eps)
    bands = band_structure_base(nothing, kgrid, comp, test; test_eigen=test_eigen)
    # post-process : taking real part
    return  BandStructureGrid{Float64}(kgrid, POST(bands), 0im.*bands)
end


function band_structure(
    kpath::Vector{Pair{String,Vector{T}}},
    DQ::Union{DynamicalMatrixQ,DynamicalMatricesQ,DynamicalMatrixKp};
    Δk = 0.02,
    test_eigen = true,
    eps = 1e-8
    )::BandStructure  where {T<:Number}
    return band_structure( kpath, DQ, EIGEN2, all_non_negative_real, SQRT_EIGEN2;
                           Δk = Δk, test_eigen = test_eigen, eps = eps )
end


function band_structure(
    kgrid::Array{Float64,2},
    DQ::Union{DynamicalMatrixQ,DynamicalMatricesQ,DynamicalMatrixKp};
    test_eigen = true,
    eps = 1e-8
    )::BandStructureGrid
    return band_structure( kgrid, DQ, EIGEN2, (x,e)->all(x.>-e), x->sqrt.(real.(x));
                           test_eigen = test_eigen, eps = eps )
end



##: ------------------------------------------------------
##: ------------------------------------------------------


# specially for hopping hamiltonian
function EIGEN(
    λ,A,B,
    HQ::HoppingHamUnion
    )
    K = (1-λ).*A.+(λ).*B
    M = convert(Matrix{ComplexF64}, HoppingHamiltonian(K,HQ))
    return ( real.(K) , _eigen_no_markers_(M) )
end


# specially for LCAO hamiltonian
function EIGEN_LCAO(
    λ,A,B,
    HSQ::Tuple{T1,T2}  #* for LCAO Hamiltonian, (H,S)
    ) where { T1<:HoppingHamUnion, T2<:HoppingHamUnion }
    K = (1-λ).*A.+(λ).*B
    M = convert(Matrix{ComplexF64}, HoppingHamiltonian(K,HSQ[1]))
    S = convert(Matrix{ComplexF64}, HoppingHamiltonian(K,HSQ[2]))
    return ( real.(K) , _eigen_no_markers_(M,S) )
end


# specially for hopping hamiltonian
function EIGENM(
    λ,A,B,
    HQ::HoppingHamUnion,
    OP::Vector
    )
    K = (1-λ).*A.+(λ).*B
    M = convert(Matrix{ComplexF64}, HoppingHamiltonian(K,HQ))
    return ( real.(K) , _eigen_with_markers_(K,M,OP) )
end


# specially for hopping hamiltonian
function EIGENM(
    λ,A,B,
    HSQ::Tuple{T1,T2},
    OP::Vector
    ) where { T1<:HoppingHamUnion, T2<:HoppingHamUnion }
    K = (1-λ).*A.+(λ).*B
    M = convert(Matrix{ComplexF64}, HoppingHamiltonian(K,HSQ[1]))
    S = convert(Matrix{ComplexF64}, HoppingHamiltonian(K,HSQ[2]))
    return ( real.(K) , _eigen_with_markers_(K,M,S,OP) )
end


# specially for hopping hamiltonian and BdG hamiltonian
@inline REAL_EIGEN(EIGEN_result) = ( EIGEN_result[1], real.(EIGEN_result[2]) )


# specially for hopping hamiltonian and BdG hamiltonian
@inline function all_real(list_EIGEN_results, ϵ0)
    t1(t) = all(abs.(imag.(t)).<ϵ0)
    t2(t) = all(abs.(imag.(t[:,1])).<ϵ0)
    test(t) = (ndims(t)==1) ? t1(t) : t2(t)
    return all( test(last(eig)) for eig ∈ list_EIGEN_results )
end


function band_structure(
    kpath::Vector{Pair{String,Vector{T}}},
    HQ::HoppingHamUnion;
    Δk = 0.02,
    test_eigen = true,
    eps = 1e-8
    )::BandStructure  where {T<:Number}
    return band_structure( kpath, HQ, EIGEN, all_real, REAL_EIGEN;
                           Δk = Δk, test_eigen = test_eigen, eps = eps )
end


function band_structure(
    kgrid::Array{Float64,2},
    HQ::HoppingHamUnion;
    test_eigen = true,
    eps = 1e-8
    )::BandStructureGrid
    return band_structure( kgrid, HQ, EIGEN, (x,e)->all(abs.(imag.(x)).<e), x->real.(x);
                           test_eigen = test_eigen, eps = eps )
end



function band_structure(
    kpath::Vector{Pair{String,Vector{T}}},
    HSQ_tuple::Tuple{T1,T2};   #* for LCAO Hamiltonian, (H,S)
    Δk = 0.02,
    test_eigen = true,
    eps = 1e-8
    )::BandStructure where {T1<:HoppingHamUnion, T2<:HoppingHamUnion, T<:Number}
    return band_structure( kpath, HSQ_tuple, EIGEN_LCAO, all_real, REAL_EIGEN;
                           Δk = Δk, test_eigen = test_eigen, eps = eps )
end



function band_structure(
    kgrid::Array{Float64,2},
    HSQ_tuple::Tuple{T1,T2};   #* for LCAO Hamiltonian, (H,S)
    test_eigen = true,
    eps = 1e-8
    )::BandStructureGrid where {T1<:HoppingHamUnion, T2<:HoppingHamUnion}
    return band_structure( kgrid, HSQ_tuple, EIGEN_LCAO, (x,e)->all(abs.(imag.(x)).<e), x->real.(x);
                           test_eigen = test_eigen, eps = eps )
end

