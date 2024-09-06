#TODO EXAMPLES !!!!!!!

mutable struct HoppingHamiltonianSymbolicQ
    LATT::Lattice
    MATS::Dict              # var => ( mn => Φ )
    T::HoppingParameter     # spring constants Dict(:f1=>f1(x1,x2), :f2=>f2(x3,x4))
    Tparams::Dict           # parameters of functions in the above dict :f1 => (:x,:y,)
    V::Dict                 # on-site potential: :m => f(x)
    Vparams::Dict           # parameters of functions in the above dict :m => (:x,)
end


function kspace_hopping_hamiltonian(
    params::Dict,
    HSQ::HoppingHamiltonianSymbolicQ;
    boundary = 0,
    inbboxbd = ((x,A)->true)
    )::HoppingHamiltonianQ
    T = dispatch_params(params, HSQ.T.HPLN, HSQ.Tparams, 0.0)
    V = dispatch_params(params, HSQ.V     , HSQ.Vparams, 0.0)
    kspace_hopping_hamiltonian( HoppingParameter(HSQ.LATT.UC,T), V, HSQ.LATT;
                                inbboxbd=inbboxbd )
end

function kspace_hopping_hamiltonian_monolayer(
    params::Dict,
    HSQ::HoppingHamiltonianSymbolicQ;
    EPS=1e-6
    )::HoppingHamiltonianQ
    @assert verify_monolayer(HSQ.LATT, EPS)
    #  pbc[1] and pbc[2] are taken as 'true' regardless of their actual settings
    kspace_hopping_hamiltonian( params, HSQ, boundary=0, inbboxbd=inbbox3 )
end


function kspace_hopping_hamiltonian_cylinder(
    params::Dict,
    HSQ::HoppingHamiltonianSymbolicQ;
    boundary=1,
    EPS=1e-6
    )::HoppingHamiltonianQ
    @assert verify_cylinder(HSQ.LATT, boundary, EPS)
    inbboxbd_cylinder( v, A ) = ( boundary==1 ? inbbox3(v,A) && inbbox1(v,A) : inbbox3(v,A) && inbbox2(v,A) )
    kspace_hopping_hamiltonian( params, HSQ, boundary=boundary, inbboxbd=inbboxbd_cylinder )
end


function dHQda(
    P1::Tuple,
    P2::Tuple,
    da::T,
    latt::Lattice
    ) where { T<:Number }
    (t1,v1) = P1
    (t2,v2) = P2
    HQ1 = kspace_hopping_hamiltonian(HoppingParameter(latt.UC,t1), v1, latt)
    HQ2 = kspace_hopping_hamiltonian(HoppingParameter(latt.UC,t2), v2, latt)
    MAT0 = zeros(ComplexF64, size(first(values(HQ1.MAT))))
    return Dict( k => ((get(HQ2.MAT,k,MAT0).-get(HQ1.MAT,k,MAT0))./da)
                 for k ∈ sort(sck(HQ1.MAT) ∪ sck(HQ2.MAT)) ) |> delete_empty
end


function kspace_hopping_hamiltonian_symbolic(
    TF::Dict{Symbol,V1},
    Tparams::Dict{Symbol,VP1},
    VF::Dict{Symbol,V2},
    Vparams::Dict{Symbol,VP2},
    latt::Lattice;
    da=1.0
    )::HoppingHamiltonianSymbolicQ where { V1<:Function, V2<:Function, VP1, VP2 }
    # positive random value of all parameters
    p0 = Dict( Dict(k=>rand() for sp ∈ sck(Tparams) for k ∈ Tparams[sp]) ∪
               Dict(k=>rand() for sp ∈ sck(Vparams) for k ∈ Vparams[sp]) )
    # construct T0,V0
    t0 = dispatch_params(p0,TF,Tparams,0.0)
    v0 = dispatch_params(p0,VF,Vparams,0.0)
    # construct T1,V1 from t0,v0
    @inline modify_param(kmod) = setindex!(copy(p0),p0[kmod]+da,kmod)
    @inline mt(kmod) = dispatch_params(modify_param(kmod),TF,Tparams,0.0)
    @inline mv(kmod) = dispatch_params(modify_param(kmod),VF,Vparams,0.0)
    # modify each key χ in p0 and compute the change dHQda
    MATS = Dict( χ => dHQda((t0,v0),(mt(χ),mv(χ)),da,latt) for χ ∈ sck(p0) )
    return HoppingHamiltonianSymbolicQ( copy(latt), MATS,
                                        HoppingParameter(latt.UC,copy(TF)), copy(Tparams),
                                        copy(VF), copy(Vparams) )
end
