mutable struct HoppingHamiltonianR{TH<:AbstractMatrix}
    LATT::Lattice
    H::TH
    EIG::EigenModes
    T::HoppingParameter     # spring constants Dict(:f1=>1.0, :f2=>1.2)
    V::Dict                 # on-site potential
end

@inline is_HoppingHamiltonianR(hhr) = (  Set(fieldnames(HoppingHamiltonianR))==Set(fieldnames(typeof(hhr)))
                                      && is_Lattice(hhr.LATT)
                                      && isa(hhr.H,  AbstractMatrix)
                                      && is_EigenModes(hhr.EIG)
                                      && is_HoppingParameter(hhr.T)
                                      && isa(hhr.V,  Dict) )


check_compat(HR::HoppingHamiltonianR) = ( check_compat(HR.LATT)
                                       && check_compat(HR.EIG)
                                       && check_compat(HR.T)
                                       && size(HR.H,1)==size(HR.H,2)==length(HR.EIG.ω) )


copy(HR::HoppingHamiltonianR{TH}) where { TH<:AbstractMatrix } = HoppingHamiltonianR{TH}(
                        copy(HR.LATT), copy(HR.H), copy(HR.EIG), copy(HR.T), copy(HR.V) )


@inline convert(::Type{T}, x::T) where {T<:HoppingHamiltonianR} = x


function convert(::Type{T}, hhr) where {T<:HoppingHamiltonianR}
    @assert is_HoppingHamiltonianR(hhr)
    HoppingHamiltonianR{typeof(hhr.H)}( convert(Lattice,            hhr.LATT  ),
                                        convert(AbstractMatrix,     hhr.H     ),
                                        convert(EigenModes,         hhr.EIG   ),
                                        convert(HoppingParameter,   hhr.T     ),
                                        convert(Dict,               hhr.V     )  )
end


@inline function all_real_v1(eigenvalues, ϵ0)
    test(t) = all(abs(imag(t))<ϵ0)
    return all(test.(eigenvalues))
end


function generate_blocks_for_HR(latt::Lattice)
    blocks = block_ranges( orbit_numbers(latt)[inner_site_sublattices(latt)] )
    dimA = sum(length.(blocks))
    return blocks, dimA
end

# -----------------------------------------------------------


function rspace_hopping_hamiltonian(
    T::HoppingParameter,
    V::Dict,
    latt::Lattice;
    solve=false,
    test=true,
    EPS=1e-10
    )::HoppingHamiltonianR

    BLK, dimH = generate_blocks_for_HR(latt)
    Nξ        = orbit_numbers( latt.UC )
    @assert dimH == (total_num_orbits(latt.UC)*num_unitcells(latt))
    #id0EqV, C = position_index( latt.EqV )

    # -----------------------------------------------------------
    # Δ, reserved for Slater-Koster, should be computed
    # instead of assigned (as in kspace)
    @inline sl(i)    = latt.SL[i]
    @inline NξSL(i)  = Nξ[latt.SL[i]]
    @inline MSL(i)   = latt.UC.m[latt.SL[i]]
    Λ0(ii,jj,hk)     = hopping_parameter_matrix(
                            hk[1], 
                            (latt.R0[:,ii].-latt.R0[:,jj], latt.R0[:,ii], latt.R0[:,jj], NξSL(ii), NξSL(jj)), 
                            hk[2]
                       )
    Λ(ii,jj,sp)      = Λ0(ii, jj, get(T.HPLN,sp,(:DEFAULT,(0.0,))))
    Λd(ii,M_dummy)   = on_site_potential(NξSL(ii), get(V,MSL(ii),0.0))

    # -----------------------------------------------------------
    is_Hofstadter = T.HPLN[first(keys(T.HPLN))][1] ∈ [:HOFALLX, :HOFALLY, :HOFALL]
    H = rspace_matrix(latt, BLK, Λ, Λd; correct_boundary_hopping_phase=is_Hofstadter)

    test_eigen_HR( e ) = ( ! test ) || all_real_v1( e, EPS )
    eig = ( solve ? solve_eigenmodes( H, 
                                      EPS=EPS, 
                                      eigen_test=test_eigen_HR, 
                                      eigen_process=(x->x) )
                  : EigenModes(dimH) )

    return HoppingHamiltonianR{typeof(H)}( copy(latt), H, eig, copy(T), copy(V) )
end


function rspace_hopping_hamiltonian_monolayer(
    T::HoppingParameter,
    V::Dict,
    latt::Lattice;
    solve=false,
    test=true,
    EPS=1e-10
    )::HoppingHamiltonianR
    @assert verify_monolayer(latt, EPS)
    # call the ordinary construction
    return rspace_hopping_hamiltonian( T, V, latt,
                                       solve=solve, test=test, EPS=EPS )
end


function solve_eigenmodes!(HR::HoppingHamiltonianR; EPS=1e-8, test=true)
    test_eigen_HR( e ) = all_real_v1( e, EPS ) || ( ! test )
    HR.EIG = solve_eigenmodes( HR.H, EPS=EPS, eigen_test=test_eigen_HR, eigen_process=(x->x) )
    return nothing
end


function rspace_zero_hopping_hamiltonian(latt::Lattice)
    return rspace_hopping_hamiltonian( zero_HoppingParameter(latt),
                                       zero_onsite_potential(latt),
                                       latt,
                                       solve=false, test=false )
end

function rspace_zero_hopping_hamiltonian_monolayer(latt::Lattice)
    return rspace_hopping_hamiltonian_monolayer( zero_HoppingParameter(latt),
                                                 zero_onsite_potential(latt),
                                                 latt,
                                                 solve=false, test=false )
end

function rspace_chemical_potential(latt::Lattice)
    return rspace_hopping_hamiltonian( zero_HoppingParameter(latt),
                                       one_onsite_potential(latt),
                                       latt,
                                       solve=false, test=false )
end

function rspace_chemical_potential_monolayer(latt::Lattice)
    return rspace_hopping_hamiltonian_monolayer( zero_HoppingParameter(latt),
                                                 one_onsite_potential(latt),
                                                 latt,
                                                 solve=false, test=false )
end


function rspace_onsite_potential(latt::Lattice,Vonsite::Vector)
    HR0 = rspace_zero_hopping_hamiltonian(latt)
    @assert size(HR.H,1) == length(Vonsite)
    HR0.H = spdiagm(0=>Vonsite)
    return HR0
end

function rspace_onsite_potential_monolayer(latt::Lattice,Vonsite::Vector)
    HR0 = rspace_zero_hopping_hamiltonian_monolayer(latt)
    @assert size(HR.H,1) == length(Vonsite)
    HR0.H = spdiagm(0=>Vonsite)
    return HR0
end


function rspace_onsite_potential(latt::Lattice,Vonsite_func::Function)
    HR0 = rspace_zero_hopping_hamiltonian(latt)
    BLK, dimH = generate_blocks_for_HR(latt)
    for (i,id) ∈ enumerate(inner_site_id(latt))
        HR0.H[BLK[i],BLK[i]] .= Vonsite_func(latt.R0[:,id], latt.UC.ξ[latt.SL[id]])
    end
    return HR0
end


function rspace_onsite_potential_monolayer(latt::Lattice,Vonsite_func::Function)
    HR0 = rspace_zero_hopping_hamiltonian_monolayer(latt)
    BLK, dimH = generate_blocks_for_HR(latt)
    Vonsite = zeros(Float64,dimH)
    for (i,id) ∈ enumerate(inner_site_id(latt))
        Vonsite[BLK[i]] = Vonsite_func(latt.R0[:,id])
    end
    HR0.H = spdiagm(0=>Vonsite)
    return HR0
end


function rspace_mean_field(latt::Lattice,MF)
    HR0 = rspace_zero_hopping_hamiltonian(latt)
    @assert size(HR.H) == size(MF)
    HR0.H = convert(typeof(HR0.H), MF)
    return HR0
end
