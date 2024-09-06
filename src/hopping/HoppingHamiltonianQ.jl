# for fast construction of
# the momentum space Hopping Hamiltonian
# H(q) = ∑_l exp(i q.r[l]) h[l]

mutable struct HoppingHamiltonianQ
    LATT::Lattice
    MAT::Dict               # mn => Φ
    T::HoppingParameter     # spring constants Dict(:f1=>1.0, :f2=>1.2)
    V::Dict                 # on-site potential
end

@inline is_HoppingHamiltonianQ(hhq) = ( Set(fieldnames(HoppingHamiltonianQ))==Set(fieldnames(typeof(hhq)))
                                     && is_Lattice(          hhq.LATT )
                                     && isa(  hhq.MAT,       Dict     )
                                     && is_HoppingParameter( hhq.T    )
                                     && isa(  hhq.V,         Dict     ) )


check_compat(HQ::HoppingHamiltonianQ) = ( check_compat(HQ.LATT)
                                       && all( size(HQ.MAT[r],1)==size(HQ.MAT[r],2)==total_num_orbits(HQ.LATT)
                                               for r ∈ keys(HQ.MAT) )
                                       && all(length(r)==dimensions(HQ.LATT) for r ∈ keys(HQ.MAT))
                                       && Set(sck(HQ.T.HPLN))==Set(sck(HQ.LATT.LN.SPNB))
                                       && HQ.T.UC == HQ.LATT.UC && check_compat(HQ.T) )


copy(HQ::HoppingHamiltonianQ) = HoppingHamiltonianQ( copy(HQ.LATT), copy(HQ.MAT),
                                                     copy(HQ.T),    copy(HQ.V) )


@inline convert(::Type{T}, x::T) where {T<:HoppingHamiltonianQ} = x


function convert(::Type{T}, hhq) where {T<:HoppingHamiltonianQ}
    @assert is_HoppingHamiltonianQ(hhq)
    HoppingHamiltonianQ( convert( Lattice,          hhq.LATT ),
                         convert( Dict,             hhq.MAT  ),
                         convert( HoppingParameter, hhq.T    ),
                         convert( Dict,             hhq.V    )  )
end


@inline consistency_check(hq::HoppingHamiltonianQ; TOL=1e-8) = consistency_check(hq.MAT,TOL=TOL)


# favourite function
function HoppingHamiltonian(
    q::Vector,
    HQ::HoppingHamiltonianQ
    )
    return qMAT(q, HQ.MAT, HQ.LATT.UC.a)
end


# -----------------------------------------------------------


function kspace_hopping_hamiltonian(
    T::HoppingParameter,
    V::Dict,
    latt::Lattice;
    inbboxbd = ((x,A)->true)
    )::HoppingHamiltonianQ

    # verify
    @assert check_compat(T)
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
    Λ0(ii,jj,Δ,hk) = hopping_parameter_matrix(
                            hk[1], 
                            (Δ, latt.R0[:,ii], latt.R0[:,jj], Nξ[ii], Nξ[jj]), 
                            hk[2]
                     )
    Λ(ii,jj,Δ,sp)  = Λ0(ii, jj, Δ, get(T.HPLN, sp, (:DEFAULT,(0.0,))))
    Λd(ii)         = on_site_potential(Nξ[ii], get(V,latt.UC.m[ii],0.0))
    # -----------------------------------------------------------

    # compute matrices for all reciprocal vectors
    mat = kspace_matrix( latt, BLK, Λ, Λd, inbboxbd )

    return HoppingHamiltonianQ( copy(latt), mat, copy(T), copy(V) )
end


function kspace_hopping_hamiltonian_monolayer(
    T::HoppingParameter,
    V::Dict,
    latt::Lattice;
    EPS=1e-6
    )::HoppingHamiltonianQ
    @assert verify_monolayer(latt, EPS)
    #!  pbc[1] and pbc[2] are taken as 'true' regardless of their actual settings
    kspace_hopping_hamiltonian( T, V, latt, 
                                inbboxbd=inbbox3 )
end


function kspace_hopping_hamiltonian_cylinder(
    T::HoppingParameter,
    V::Dict,
    latt::Lattice;
    boundary=1,
    EPS=1e-6
    )::HoppingHamiltonianQ
    @assert verify_cylinder(latt, boundary, EPS)
    inbboxbd_cylinder( v, A ) = ( boundary==1 ? inbbox3(v,A) && inbbox1(v,A) : inbbox3(v,A) && inbbox2(v,A) )
    kspace_hopping_hamiltonian( T, V, latt, 
                                inbboxbd=inbboxbd_cylinder )
end


function kspace_zero_hopping_hamiltonian(latt::Lattice)
    return kspace_hopping_hamiltonian( zero_HoppingParameter(latt), zero_onsite_potential(latt), latt )
end


function kspace_zero_hopping_hamiltonian_monolayer(latt::Lattice)
    return kspace_hopping_hamiltonian_monolayer( zero_HoppingParameter(latt),
                                                 zero_onsite_potential(latt),
                                                 latt )
end


function kspace_zero_hopping_hamiltonian_cylinder(latt::Lattice; boundary=1)
    return kspace_hopping_hamiltonian_cylinder( zero_HoppingParameter(latt),
                                                zero_onsite_potential(latt),
                                                latt,
                                                boundary=boundary )
end


function kspace_chemical_potential(latt::Lattice)
    return kspace_hopping_hamiltonian( zero_HoppingParameter(latt),
                                       one_onsite_potential(latt),
                                       latt )
end


function kspace_chemical_potential_monolayer(latt::Lattice)
    return kspace_hopping_hamiltonian_monolayer( zero_HoppingParameter(latt),
                                                 one_onsite_potential(latt),
                                                 latt )
end


function kspace_chemical_potential_cylinder(latt::Lattice; boundary=1)
    return kspace_hopping_hamiltonian_cylinder( zero_HoppingParameter(latt),
                                                one_onsite_potential(latt),
                                                latt,
                                                boundary=boundary )
end
