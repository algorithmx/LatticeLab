#: -------------------------------------------------------------
#: definition of struct
#: -------------------------------------------------------------

"""
    mutable struct DynamicalMatrixQ
        LATT::Lattice
        MAT::Dict       # R_mnl => Φ_mnl
        SP::Dict        # spring constants Dict(:f1=>1.0, :f2=>1.2)
        M::Dict         # mass Dict(:mA=>1.0, :mB=>1.2)
    end

======

# Note :

Data struxture `DynamicalMatrixQ` is used for fast construction of 
the momentum space Dynamical matrix 

D(q) = ∑_mnl exp(i q.R_mnl) Φ_mnl

where R_mnl and Φ_mnl are the key and value of Dict `DynamicalMatrixQ.MAT`.

"""
mutable struct DynamicalMatrixQ
    LATT::Lattice
    MAT::Dict       # R_mnl => Φ_mnl
    SP::Dict        # spring constants Dict(:f1=>1.0, :f2=>1.2)
    M::Dict         # mass Dict(:mA=>1.0, :mB=>1.2)
end

#: -------------------------------------------------------------
#: cpoy, convert and chekc_compat and consistency
#: -------------------------------------------------------------


#TODO addition of DQ with identical LATT and M
#TODO merge SP field

@inline is_DynamicalMatrixQ(dmq) = (   Set(fieldnames(DynamicalMatrixQ))==Set(fieldnames(typeof(dmq)))
                                    && is_Lattice(dmq.LATT)
                                    && isa(dmq.MAT, Dict)
                                    && isa(dmq.SP,  Dict)
                                    && isa(dmq.M,   Dict) )


check_compat(DQ::DynamicalMatrixQ) = ( check_compat(DQ.LATT)
                                    && all(size(DQ.MAT[r],1)==size(DQ.MAT[r],2)==dimensions(DQ.LATT)
                                            for r ∈ keys(DQ.MAT))
                                    && all(length(r)==dimensions(DQ.LATT) for r ∈ keys(DQ.MAT))
                                    && Set(collect(keys(DQ.SP)))==Set(keys(DQ.LATT.LN.SPNB))
                                    && all([ ( (typeof(DQ.SP[k])<:Number)
                                            || (typeof(DQ.SP[k])<:Tuple && isa(first(DQ.SP[k]),Number)) )
                                             for k ∈ keys(DQ.SP) ])
                                    && Set(collect(keys(DQ.M)) )==Set(DQ.LATT.UC.m) )


copy(DQ::DynamicalMatrixQ) = DynamicalMatrixQ( copy(DQ.LATT), copy(DQ.MAT),
                                               copy(DQ.SP), copy(DQ.M) )


@inline convert(::Type{T}, x::T) where {T<:DynamicalMatrixQ} = x


function convert(::Type{T}, dmq) where {T<:DynamicalMatrixQ}
    @assert is_DynamicalMatrixQ(dmq)
    DynamicalMatrixQ( convert( Lattice, dmq.LATT ),
                      convert( Dict,    dmq.MAT  ),
                      convert( Dict,    dmq.SP   ),
                      convert( Dict,    dmq.M    )  )
end


function consistency_check(
    MAT::Dict{Vector{Int64},V};
    TOL = 1e-8
    ) where { V }
    Bravais = sck(MAT)
    t1 = [(-v ∈ Bravais) && norm(MAT[-v]'.-MAT[v])<TOL for v ∈ Bravais if v[1]>=0 ]
    t2 = [(-v ∈ Bravais) && norm(MAT[-v]'.-MAT[v])<TOL for v ∈ Bravais if v[2]>=0 ]
    t3 = [(-v ∈ Bravais) && norm(MAT[-v]'.-MAT[v])<TOL for v ∈ Bravais if v[3]>=0 ]
    return all(t1) && all(t2) && all(t3)
end

@inline consistency_check(dq::DynamicalMatrixQ; TOL=1e-8) = consistency_check(dq.MAT,TOL=TOL)


#: -------------------------------------------------------------
#: D(q)
#: -------------------------------------------------------------


function DynamicalMatrix(
    q::Vector,
    DQ::DynamicalMatrixQ
    )
    # the Fourier sum
    # qMAT(q,M,a) = ∑_mnl exp[i q . (m*a[1]+n*a[2]+l*a[3])] * M[l,m,n]
    return qMAT(q, DQ.MAT, DQ.LATT.UC.a)
end


#: -------------------------------------------------------------
#: ∂D∂qα(q)
#: -------------------------------------------------------------


function ∂DynamicalMatrix∂qα(
    q::Vector,
    α::Int,
    DQ::DynamicalMatrixQ
    )
    # the Fourier sum
    # ∂qMAT∂q(q,M,a) = ∑_mnl (i) (m*a[1]+n*a[2]+l*a[3])_α * exp[i q . (m*a[1]+n*a[2]+l*a[3])] * M[l,m,n]
    return ∂qMAT∂q(q, α, DQ.MAT, DQ.LATT.UC.a)
end


#: -------------------------------------------------------------

#: -------------------------------------------------------------

@inline generate_blocks_for_DQ(latt::Lattice) = generate_blocks_for_DQ(latt.UC)

function enforce_acoustic_sum_rule(
    M::Dict,
    latt::Lattice,
    Dq0
    )
    BLK, _ = generate_blocks_for_DQ(latt)
    Nsubl = num_sublattice(latt)
    Msqrt = [ sqrt(get(M,latt.UC.m[i],1.0)) for i ∈ 1:Nsubl ]
    D = zeros(Float64, size(Dq0))
    for s ∈ 1:Nsubl
        D[BLK[s],BLK[s]] .= LinearAlgebra.Matrix((-1.0/Msqrt[s]).*sum(Msqrt[sp].*Dq0[BLK[s],BLK[sp]] for sp ∈ 1:Nsubl))
        #NOTE the above line doesn't guarantee flexural mode
    end
    return sparse(D)
end


"""
    function kspace_dynamical_matrix(
        K::Dict,             # spring constants Dict(:f1=>1.0, :f2=>1.2) 
        M::Dict,             # mass Dict(:mA=>1.0, :mB=>1.2)
        latt::Lattice;       # 
        boundary = 0,        # 
        inbboxbd = ((x,A)->true)
    )::DynamicalMatrixQ


# Note :

Return : `DynamicalMatrixQ(latt, mat, K, M)`

where `mat = kspace_matrix(latt, BLK, Λ, Λd, inbboxbd)`


"""
function kspace_dynamical_matrix( #NOTE test Kagome fail
    K::Dict,
    M::Dict,
    latt::Lattice;
    boundary = 0,             #+  ???
    inbboxbd = ((x,A)->true)  #+  ???
    )::DynamicalMatrixQ

    # verify
    @assert check_zeroth_unitcell(latt, inbboxbd) "Zeroth unit cell outside bounding box."

    # block structure for hopping Hamiltonian
    BLK, dimsD = generate_blocks_for_DQ(latt)
    dim = dimensions(latt)

    # ---------------------------------------------------------------------
    # block constructors
    onsite_label = (:ONSITE ∈ keys(K) ? :ONSITE : :DEFAULT)
    onsite_param = get(K,:ONSITE,(0.0,))
    @inline msqrt(subl)  = sqrt(get(M,latt.UC.m[subl],1.0))
    Λ0(ii,jj,Δ,hk) = ( (1.0/(msqrt(ii)*msqrt(jj))) .* interatomic_force_matrix(hk[1],(Δ,ii,jj),hk[2]) )
    Λ(ii,jj,Δ,sp) = Λ0( ii, jj, Δ, get(K,sp,(:DEFAULT,(0.0,))) )
    Λd(ii) = interatomic_force_matrix(onsite_label,(latt.UC.a[:,1],ii,ii),onsite_param)
    # ---------------------------------------------------------------------

    # compute matrices for all reciprocal vectors
    mat = kspace_matrix(latt, BLK, Λ, Λd, inbboxbd)

    # enforce the acoustic sum rule
    Dq0 = sum(values(mat))
    D = enforce_acoustic_sum_rule(M, latt, Dq0)
    OOO = zeros(Int64,dim)
    if OOO ∈ keys(mat)
        mat[OOO] = sparse(mat[OOO].+D)
    else
        mat[OOO] = sparse(D)
    end

    return DynamicalMatrixQ( copy(latt), mat, copy(K), copy(M) )
end


function kspace_dynamical_matrix_monolayer(
    K::Dict,
    M::Dict,
    latt::Lattice;
    EPS=1e-6
    )::DynamicalMatrixQ
    @assert verify_monolayer(latt, EPS)
    #  pbc[1] and pbc[2] are taken as 'true' regardless of their actual settings
    kspace_dynamical_matrix( K, M, latt, boundary=0, inbboxbd=inbbox3 )
end


function kspace_dynamical_matrix_cylinder(
    K::Dict,
    M::Dict,
    latt::Lattice;
    boundary=:x,
    EPS=1e-6
    )::DynamicalMatrixQ
    @assert verify_cylinder(latt, boundary, EPS)
    inbboxbd_cylinder(v,A) = ( boundary==1 ? inbboxt(v,A,[1,3]) : inbboxt(v,A,[2,3]) )
    kspace_dynamical_matrix( K, M, latt, boundary=boundary, inbboxbd=inbboxbd_cylinder )
end


function kspace_dynamical_matrix_ribbon_2D(
    K::Dict,
    M::Dict,
    latt::Lattice;
    boundary=:x
    )::DynamicalMatrixQ
    @assert verify_cylinder_2d(latt, boundary)
    inbboxbd_cylinder_2d(v,A) = ( boundary==1 ? inbboxt(v,A,[1]) : inbboxt(v,A,[2]) )
    kspace_dynamical_matrix(K, M, latt, boundary=boundary, inbboxbd=inbboxbd_cylinder_2d)
end


# ----------------------------------------------------
