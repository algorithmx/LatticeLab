mutable struct DynamicalMatrixR{T<:AbstractMatrix}
    LATT::Lattice
    D::T
    EIG::EigenModes
    SP::Dict            # spring constants Dict(:f1=>1.0, :f2=>1.2)
    M::Dict             # mass Dict(:mA=>1.0, :mB=>1.2)
end

#NOTE rules for the spring constant dict
#NB example 1 : force-field along link
# Dict{Symbol,Number}
# Dict(:f1=>1.0, :f2=>1.2)
#NB example 2 : force-field perpendicular to link
# Dict{Symbol,Union{Tuple,Number}}
# Dict(:f1=>(1.0,[0.2,zhat]), :f2=>1.2)
# Tuple ordering : (spring const along link, [spring const perp to link, perp direction] , ignore )

@inline is_DynamicalMatrixR(dmr) = (   Set(fieldnames(DynamicalMatrixR))==Set(fieldnames(typeof(dmr)))
                                    && is_Lattice(dmr.LATT)
                                    && isa(dmr.D,  AbstractMatrix)
                                    && is_EigenModes(dmr.EIG)
                                    && isa(dmr.SP, Dict)
                                    && isa(dmr.M,  Dict) )


check_compat(DR::DynamicalMatrixR) = ( check_compat(DR.LATT)
                                    && check_compat(DR.EIG)
                                    && size(DR.D,1)==length(DR.EIG.ω)==(num_inner_sites(DR.LATT)*dimensions(DR.LATT))
                                    && Set(collect(keys(DR.SP)))==Set(keys(DR.LATT.LN.SPNB))
                                    && all([ ( (typeof(DR.SP[k])<:Number)
                                            || (typeof(DR.SP[k])<:Tuple && isa(first(DR.SP[k]),Number)) )
                                             for k ∈ keys(DR.SP) ])
                                    && Set(collect(keys(DR.M)) )==Set(DR.LATT.UC.m) )


copy(DR::DynamicalMatrixR{T}) where { T<:AbstractMatrix } =
    DynamicalMatrixR{T}( copy(DR.LATT), copy(DR.D), copy(DR.EIG),
                         copy(DR.SP), copy(DR.M) )

@inline convert(::Type{T}, x::T) where {T<:DynamicalMatrixR} = x


function convert(::Type{T}, dmr) where {T<:DynamicalMatrixR}
    @assert is_DynamicalMatrixR(dmr)
    DynamicalMatrixR{typeof(dmr.D)}( convert(Lattice,        dmr.LATT  ),
                                     convert(AbstractMatrix, dmr.D     ),
                                     convert(EigenModes,     dmr.EIG   ),
                                     convert(Dict,           dmr.SP    ),
                                     convert(Dict,           dmr.M     )  )
end


@inline function all_non_negative_real_v1(eigenvalues_squared, ϵ0)
    test(t) = all(real(t)>(-ϵ0)) && all(abs(imag(t))<ϵ0)
    return all(test.(eigenvalues_squared))
end


function generate_blocks_for_DR(latt::Lattice)
    blocks = block_ranges( dimensions(latt) ⨰ num_inner_sites(latt) )
    dimA = sum(length.(blocks))
    return blocks, dimA
end


function rspace_dynamical_matrix(
    K::Dict,
    M::Dict,
    latt::Lattice;
    solve=false,
    test=true,
    EPS=1e-10
    )::DynamicalMatrixR
    BLK, dimD = generate_blocks_for_DR(latt)
    @assert dimD == (dimensions(latt)*num_inner_sites(latt))
    id0EqV, C = position_index( latt.EqV )

    # ---------------------------------------------------------------------
    @inline msqrt(i) = sqrt( get(M,mass_site(latt,i),1.0) )
    @inline dij(i,j) = latt.R0[:,i].-latt.R0[:,j]
    @inline sl(i)    = latt.SL[i]
    Λ0(i0,j0,hk)     = (1.0/(msqrt(i0)*msqrt(j0))).*interatomic_force_matrix(hk[1],(dij(i0,j0),sl(i0),sl(j0)),hk[2])
    Λ(i0,j0,sp)      = Λ0(i0,j0,get(K,sp,(:DEFAULT,(0.0,))))
    Λd(s,A)          = ((-1.0/msqrt(s)) .* sum( msqrt(id0EqV[sprime]).*LinearAlgebra.Matrix(A[BLK[C[s]],BLK[sprime]])
                                                for sprime ∈ 1:num_inner_sites(latt) ))
    # ---------------------------------------------------------------------

    D = rspace_matrix( latt, BLK, Λ, Λd )

    sqrt_eigen( x ) = sqrt(max(0.0,real(x)))
    test_eigen_DR( e ) = all_non_negative_real_v1( e, EPS ) || ( ! test )
    eig = ( solve ? solve_eigenmodes( D, EPS=EPS, eigen_test=test_eigen_DR, eigen_process=sqrt_eigen )
                  : EigenModes(dimD) )

    return DynamicalMatrixR{typeof(D)}( copy(latt), D, eig, copy(K), copy(M) )
end



function rspace_dynamical_matrix_monolayer(
    K::Dict,
    M::Dict,
    latt::Lattice;
    solve=false,
    test=true,
    EPS=1e-10
    )::DynamicalMatrixR
    @assert verify_monolayer(latt, EPS)
    # call the ordinary construction
    return rspace_dynamical_matrix( K, M, latt, solve=solve, test=test, EPS=EPS )
end



function fold_DymamicalMatrixR(
    D::T,
    LATT::Lattice;
    folder = (x->tr(abs.(x)))
    )::Matrix where { T<:AbstractMatrix }
    BLK, dimD = generate_blocks_for_DR(LATT)
    L = length(BLK)
    W = zeros(Float64,L,L)
    for i = 1:L
        for j = 1:L
            W[i,j] = folder(D[BLK[i],BLK[j]])
        end
    end
    return W
end

function fold_DymamicalMatrixR(
    DR::DynamicalMatrixR;
    folder = (x->tr(abs.(x)))
    )::Matrix
    return fold_DymamicalMatrixR(DR.D, DR.LATT, folder=folder)
end
