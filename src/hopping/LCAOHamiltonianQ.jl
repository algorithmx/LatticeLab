# for fast construction of
# the momentum space Hopping Hamiltonian
# H(q) = ∑_l exp(i q.r[l]) h[l]

mutable struct LCAOHamiltonianQ
    LATT::Lattice
    MAT::Dict               # mn => Φ
    OVERLAP::Dict           # mn => Φ
    T::HoppingParameter     # spring constants Dict(:f1=>1.0, :f2=>1.2)
    S::HoppingParameter     # spring constants Dict(:f1=>1.0, :f2=>1.2)
    V::Dict                 # on-site potential
end

@inline is_LCAOHamiltonianQ(hhq) = ( Set(fieldnames(LCAOHamiltonianQ))==Set(fieldnames(typeof(hhq)))
                                     && is_Lattice(          hhq.LATT )
                                     && isa(  hhq.MAT,       Dict     )
                                     && isa(  hhq.OVERLAP,   Dict     )
                                     && is_HoppingParameter( hhq.T    )
                                     && is_HoppingParameter( hhq.S    )
                                     && isa(  hhq.V,         Dict     ) 
                                    )


check_compat(HQ::LCAOHamiltonianQ) =  ( check_compat(HQ.LATT)
                                       && all( size(HQ.MAT[r],1)==size(HQ.MAT[r],2)==total_num_orbits(HQ.LATT)
                                               for r ∈ keys(HQ.MAT) )
                                       && all( size(HQ.OVERLAP[r],1)==size(HQ.OVERLAP[r],2)==total_num_orbits(HQ.LATT)
                                               for r ∈ keys(HQ.OVERLAP) )
                                       && all(length(r)==dimensions(HQ.LATT) for r ∈ keys(HQ.MAT))
                                       && all(length(r)==dimensions(HQ.LATT) for r ∈ keys(HQ.OVERLAP))
                                       && Set(sck(HQ.T.HPLN))==Set(sck(HQ.LATT.LN.SPNB))
                                       && Set(sck(HQ.S.HPLN))==Set(sck(HQ.LATT.LN.SPNB))
                                       && HQ.T.UC == HQ.LATT.UC && check_compat(HQ.T)
                                       && HQ.S.UC == HQ.LATT.UC && check_compat(HQ.S)
                                      )


copy(HQ::LCAOHamiltonianQ) = LCAOHamiltonianQ( copy(HQ.LATT), copy(HQ.MAT), copy(HQ.OVERLAP),
                                               copy(HQ.T),    copy(HQ.S),   copy(HQ.V) )


@inline convert(::Type{T}, x::T) where {T<:LCAOHamiltonianQ} = x


function convert(::Type{T}, hhq) where {T<:LCAOHamiltonianQ}
    @assert is_LCAOHamiltonianQ(hhq)
    LCAOHamiltonianQ(    convert( Lattice,          hhq.LATT ),
                         convert( Dict,             hhq.MAT  ),
                         convert( Dict,             hhq.OVERLAP  ),
                         convert( HoppingParameter, hhq.T    ),
                         convert( HoppingParameter, hhq.S    ),
                         convert( Dict,             hhq.V    )  )
end


@inline consistency_check(hq::LCAOHamiltonianQ; TOL=1e-8) = 
            (consistency_check(hq.MAT,TOL=TOL) && consistency_check(hq.OVERLAP,TOL=TOL))


eigen2mat(e,v) = (v .* transpose(e)) * v'
function matrix_sqrt(M::Matrix)
    eig = eigen(M)
    eigen2mat(sqrt.(complex.(eig.values)), eig.vectors)
end

##
#for i=1:200
#    A = rand(100,100) .+ im*rand(100,100)
#    A = A .+ A'
#    sqrtA = matrix_sqrt(A)
#    @info norm(sqrtA*sqrtA.-A)
#end
##

##% -----------------------------------------------------------


# core function
# return (H, S)
function LCAOHamiltonian(
    q::Vector,
    HQ::LCAOHamiltonianQ
    )
    #Sinv     = inv(qMAT(q, HQ.OVERLAP, HQ.LATT.UC.a))
    #sqrtSinv = matrix_sqrt(Sinv)
    #return sqrtSinv * qMAT(q, HQ.MAT, HQ.LATT.UC.a) * sqrtSinv
    return qMAT(q, HQ.MAT, HQ.LATT.UC.a), qMAT(q, HQ.OVERLAP, HQ.LATT.UC.a)
end


##% -----------------------------------------------------------


function kspace_LCAO_hamiltonian(
    T::HoppingParameter,
    S::HoppingParameter,
    V::Dict,
    latt::Lattice;
    inbboxbd = ((x,A)->true),
    vector_identical = (x,y)->norm(x.-y)<1e-5
    )::LCAOHamiltonianQ

    # verify
    @assert check_compat(T)
    @assert check_compat(S)
    @assert check_zeroth_unitcell(latt, inbboxbd) "Zeroth unit cell outside bounding box."

    # block structure for hopping Hamiltonian
    BLK, dimH = generate_blocks_for_HQ(latt.UC)
    Nξ        = orbit_numbers( latt.UC )

    # -----------------------------------------------------------
    # block constructors
    # new versions in HoppingParameter.jl
    # old version below
    #@inline modeA(i,j,Δ,t) = Matrix( t*I, Nξ[i], Nξ[j] )
    #@inline modeB(i,j,Δ,t) = ((@assert size(t)==(Nξ[i],Nξ[j])); t)
    #@inline modeDEF(i,j,Δ,t) = zeros(ComplexF64,Nξ[i],Nξ[j])
    #χ = Dict( :A=>modeA, :B=>modeB, :DEFAULT=>modeDEF )
    #@inline Λ0(ii,jj,Δ,t) = χ[mode_H(t)](ii,jj,Δ,t)
    #Λ(ii,jj,Δ,sp) = Λ0( ii, jj, Δ, get(T.HPLN,sp,0.0) )
    #Λd(ii) = Λ0( ii, ii, 0, get(V,latt.UC.m[ii],0.0) )
    # -----------------------------------------------------------
    # hopping_parameter_matrix(type, ΔIJ, params)
    # type   = hk[1]
    # ΔIJ    = (Δ,Nξ[ii],Nξ[jj])
    # params = hk[2] , 
    #   FORMAT   (:TYPE,    (p1,p2,...) )
    #   DEFAULT  (:DEFAULT, (0.0,)      )
    Λ0(ii,jj,hk) = hopping_parameter_matrix(hk[1], (latt.R0[:,ii],latt.R0[:,jj],Nξ[ii],Nξ[jj]), hk[2])
    Λ(ii,jj,sp)  = Λ0(ii, jj, get(T.HPLN, sp, (:DEFAULT,(0.0,))))
    Λd(ii)       = on_site_potential(Nξ[ii], get(V,latt.UC.m[ii],0.0))
    # -----------------------------------------------------------
    
    # compute matrices for all reciprocal vectors
    mat = kspace_matrix( latt, BLK, Λ, Λd, inbboxbd, vector_identical )
    # -----------------------------------------------------------

    Σ(ii,jj,sp)  = Λ0(ii, jj, get(S.HPLN, sp, (:DEFAULT,(0.0,))))
    Σd(ii)       = on_site_potential(Nξ[ii], :return_zero)
    # -----------------------------------------------------------

    # compute matrices for all reciprocal vectors
    overlap = kspace_matrix( latt, BLK, Σ, Σd, inbboxbd,  )

    return LCAOHamiltonianQ( copy(latt), mat, overlap, copy(T), copy(S), copy(V) )
end


function kspace_LCAO_hamiltonian_monolayer(
    T::HoppingParameter,
    S::HoppingParameter,
    V::Dict,
    latt::Lattice;
    EPS=1e-6
    )::HoppingHamiltonianQ
    @assert verify_monolayer(latt, EPS)
    #!  pbc[1] and pbc[2] are taken as 'true' regardless of their actual settings
    kspace_LCAO_hamiltonian( T, S, V, latt, 
                             inbboxbd=inbbox3,
                             vector_identical = (x,y)->norm(x.-y)<EPS )
end


function kspace_LCAO_hamiltonian_cylinder(
    T::HoppingParameter,
    S::HoppingParameter,
    V::Dict,
    latt::Lattice;
    boundary=1,
    EPS=1e-6
    )::HoppingHamiltonianQ
    @assert verify_cylinder(latt, boundary, EPS)
    inbboxbd_cylinder( v, A ) = ( boundary==1 ? inbbox3(v,A) && inbbox1(v,A) : inbbox3(v,A) && inbbox2(v,A) )
    kspace_LCAO_hamiltonian( T, S, V, latt, 
                             inbboxbd=inbboxbd_cylinder,
                             vector_identical = (x,y)->norm(x.-y)<EPS )
end


function kspace_zero_LCAO_hamiltonian(latt::Lattice)
    return kspace_LCAO_hamiltonian( zero_HoppingParameter(latt), 
                                    one_HoppingParameter(latt), 
                                    zero_onsite_potential(latt), latt )
end


function kspace_zero_LCAO_hamiltonian_monolayer(latt::Lattice)
    return kspace_LCAO_hamiltonian_monolayer( zero_HoppingParameter(latt),
                                              one_HoppingParameter(latt), 
                                              zero_onsite_potential(latt),
                                              latt )
end


function kspace_zero_LCAO_hamiltonian_cylinder(latt::Lattice; boundary=1)
    return kspace_LCAO_hamiltonian_cylinder( zero_HoppingParameter(latt),
                                             one_HoppingParameter(latt), 
                                             zero_onsite_potential(latt),
                                             latt,
                                             boundary=boundary )
end


function kspace_chemical_potential_LCAO(latt::Lattice)
    return kspace_LCAO_hamiltonian( zero_HoppingParameter(latt),
                                    one_HoppingParameter(latt), 
                                    one_onsite_potential(latt),
                                    latt )
end


function kspace_chemical_potential_monolayer_LCAO(latt::Lattice)
    return kspace_LCAO_hamiltonian_monolayer( zero_HoppingParameter(latt),
                                              one_HoppingParameter(latt), 
                                              one_onsite_potential(latt),
                                              latt )
end


function kspace_chemical_potential_cylinder_LCAO(latt::Lattice; boundary=1)
    return kspace_LCAO_hamiltonian_cylinder( zero_HoppingParameter(latt),
                                             one_HoppingParameter(latt),
                                             one_onsite_potential(latt),
                                             latt,
                                             boundary=boundary )
end
