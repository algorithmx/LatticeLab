mutable struct DynamicalMatrixSymbolicQ
    LATT::Lattice
    MATS::Dict              # var => ( mn => Φ )
    SP::Dict                # spring constants Dict(:f1=>1.0, :f2=>1.2)
    SPparams::Dict          # parameters of functions in the above dict :f1 => (:x,:y,)
    M::Dict                 # mass Dict(:mA=>1.0, :mB=>1.2)
end


function kspace_dynamical_matrix(
    params::Dict,
    DSQ::DynamicalMatrixSymbolicQ;
    boundary = 0,
    inbboxbd = ((x,A)->true)
    )::DynamicalMatrixQ
    SP = dispatch_params(params, DSQ.SP, DSQ.SPparams, 0.0)
    M  = dispatch_params(params, DSQ.M,  DSQ.Mparams,  1.0)
    kspace_dynamical_matrix( M, SP, DSQ.LATT;
                             boundary=boundary, inbboxbd=inbboxbd )
end


function kspace_dynamical_matrix_monolayer(
    params::Dict,
    DSQ::DynamicalMatrixSymbolicQ;
    EPS=1e-6
    )::DynamicalMatrixQ
    @assert verify_monolayer(DSQ.LATT, EPS)
    #  pbc[1] and pbc[2] are taken as 'true' regardless of their actual settings
    kspace_dynamical_matrix( params, DSQ, boundary=0, inbboxbd=inbbox3 )
end


function kspace_dynamical_matrix_cylinder(
    params::Dict,
    DSQ::DynamicalMatrixSymbolicQ;
    boundary=1,
    EPS=1e-6
    )::DynamicalMatrixQ
    @assert verify_cylinder(DSQ.LATT, boundary, EPS)
    inbboxbd_cylinder( v, A ) = ( boundary==1 ? inbbox3(v,A) && inbbox1(v,A) : inbbox3(v,A) && inbbox2(v,A) )
    kspace_dynamical_matrix( params, DSQ, boundary=boundary, inbboxbd=inbboxbd_cylinder )
end


function dDQda(
    P1::Tuple,
    P2::Tuple,
    da::T,
    latt::Lattice
    ) where { T<:Number }
    (sp1,m1) = P1
    (sp2,m2) = P2
    DQ1 = kspace_dynamical_matrix(sp1, m1, latt)
    DQ2 = kspace_dynamical_matrix(sp2, m2, latt)
    MAT0 = zeros(ComplexF64, size(first(values(DQ1.MAT))))
    return Dict( k => ((get(DQ2.MAT,k,MAT0).-get(DQ1.MAT,k,MAT0))./da)
                 for k ∈ sort(sck(DQ1.MAT) ∪ sck(DQ2.MAT)) ) |> delete_empty
end



function kspace_dynamical_matrix_symbolic(
    SPF::Dict{Symbol,V1},           # force constant functions
    SPparams::Dict{Symbol,VP1},     # force constant parameters
    MASS::Dict,
    latt::Lattice;
    da=1.0
    )::DynamicalMatrixSymbolicQ where { V1<:Function, VP1}
    # positive random value of all parameters
    p0 = Dict(k=>rand() for x ∈ sck(SPparams) for k ∈ SPparams[x])
    # construct sp0, m0 dictionaries
    sp0 = dispatch_params(p0,SPF,SPparams,0.0)
    m0  = Dict(lb=>get(MASS,lb,1.0) for lb ∈latt.UC.m)
    # construct T1,V1 from T0,V0
    @inline modify_param(kmod) = setindex!(copy(p0),p0[kmod]+da,kmod)
    @inline msp(kmod) = dispatch_params(modify_param(kmod),SPF,SPparams,0.0)
    # modify each key χ in p0 and compute the change dHQda
    MATS = Dict( χ => dDQda((sp0,m0),(msp(χ),m0),da,latt) for χ ∈ sck(p0) )
    return DynamicalMatrixSymbolicQ(copy(latt), MATS, copy(SPF), copy(SPparams), m0)
end


#=

mutable struct DynamicalMatrixSymbolicQ_old
    LATT::Lattice
    MATS::Dict              # var => ( mn => Φ )
    SP::Dict                # spring constants Dict(:f1=>1.0, :f2=>1.2)
    SPparams::Dict          # parameters of functions in the above dict :f1 => (:x,:y,)
    M::Dict                 # mass Dict(:mA=>1.0, :mB=>1.2)
    Mparams::Dict           # parameters of functions in the above dict :m => (:x,)
end


function kspace_dynamical_matrix_symbolic_old(
    SPF::Dict{Symbol,V1},           # force constant functions
    SPparams::Dict{Symbol,VP1},     # force constant parameters
    MF::Dict{Symbol,V2},            # mass functions
    Mparams::Dict{Symbol,VP2},      # mass parameters
    latt::Lattice;
    da=1.0
    )::DynamicalMatrixSymbolicQ where { V1<:Function, V2<:Function, VP1, VP2 }
    # positive random value of all parameters
    p0 = Dict( Dict(k=>rand() for x ∈ sck(SPparams) for k ∈ SPparams[x]) ∪
               Dict(k=>rand() for x ∈ sck(Mparams ) for k ∈ Mparams[x] ) )
    # construct sp0, m0 dictionaries
    sp0 = dispatch_params(p0,SPF,SPparams,0.0)
    m0  = dispatch_params(p0,MF, Mparams, 1.0)
    # construct T1,V1 from T0,V0
    @inline modify_param(kmod) = setindex!(copy(p0),p0[kmod]+da,kmod)
    @inline msp(kmod) = dispatch_params(modify_param(kmod),SPF,SPparams,0.0)
    @inline  mm(kmod) = dispatch_params(modify_param(kmod),MF, Mparams, 1.0)
    # modify each key χ in p0 and compute the change dHQda
    MATS = Dict( χ => dDQda((sp0,m0),(msp(χ),mm(χ)),da,latt) for χ ∈ sck(p0) )
    return DynamicalMatrixSymbolicQ( copy(latt), MATS,
                                     copy(SPF), copy(SPparams),
                                     copy(MF), copy(Mparams) )
end

=#
