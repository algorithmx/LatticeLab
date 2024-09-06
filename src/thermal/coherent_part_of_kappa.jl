function ∂sqrtM(∂M, sqrtM)
    @assert size(∂M)==size(sqrtM)
    dimD   = size(∂M,1)
    dimD2  = dimD*dimD
    A      = zeros(ComplexF64, dimD, dimD, dimD, dimD)
    sqrtMt = transpose(sqrtM)
    for i=1:dimD
        A[i,:,i,:] .+= sqrtMt
    end
    for i=1:dimD
        A[:,i,:,i] .+= sqrtM
    end
    #TODO SLOW !
    @time reshape(inv(reshape(A,(dimD2,dimD2))) * ∂M[:], (dimD, dimD))
end

##
# using LinearAlgebra
# pD     = rand(10,10)
# sqrtD  = rand(10,10)
# ∂sqrtD = ∂sqrtM(pD, sqrtD)
# LinearAlgebra.norm(∂sqrtD*sqrtD.+sqrtD*∂sqrtD .- pD)
##


eigen2mat(e,v) = (v .* transpose(e)) * v'

function tensor_V(D, ∂D)
    #% Simoncelli et. al., Nat. Phys 
    #% eq(12), term V
    @info size(D)

    eig = LinearAlgebra.eigen(Matrix(D))
    sqrtD  = eigen2mat(sqrt.(complex.(eig.values)), eig.vectors)
    @time ∂sqrtD = ∂sqrtM(∂D, sqrtD)
    print(""); flush(stdout);
    
    (eig.vectors') * ∂sqrtD * eig.vectors
end


function _part_of_kappa(
    part::Symbol,
    qpoint_weight::Vector{Int},
    T_unit_Kelvin::Float64,
    ω_unit_THz::Matrix,
    Γ_unit_THz::Matrix, 
    V_unit_X::Array{ComplexF64,4}
    )
    @assert all(ω_unit_THz[1,1:3].≈0)
    @assert all(Γ_unit_THz[1,1:3].≈0)

    #* UNIT CONVERSION
    #% Quantity[1, "eV"]/(Quantity[1, "Kelvins"] Quantity[1, "BoltzmannConstant"]) // UnitConvert
    kBT_unit_eV = T_unit_Kelvin * 11604.5177737893316139076
    #% Quantity[1, "eV"]/(Quantity[1, "PlanckConstant"] Quantity[1, "THz"]) // UnitConvert
    ħω_unit_eV  = ω_unit_THz   .* 241.7989222163607992909
    ħΓ_unit_eV  = Γ_unit_THz   .* 241.7989222163607992909
    #%  V_unit_X : X should be adapted to the coefficient outside the summation
    #% here X is sqrt(unit DM)/(unit q) = (2π factor * THz)/(1/Angstrom)
    Nc = sum(qpoint_weight)
    coeff_to_W_per_m_per_K = 0.8768198143268912122/(Nc*(kBT_unit_eV^2))

    #% Simoncelli et. al., Nat. Phys 
    #% eq(12), second part, terms after ∑
    Nqpoints = size(ħω_unit_eV,1)
    Nbands   = size(ħω_unit_eV,2)
    @assert Nqpoints == size(Γ_unit_THz,1)
    @assert Nbands   == size(Γ_unit_THz,2)
    @assert Nqpoints == size(V_unit_X,1)
    @assert Nbands   == size(V_unit_X,2)
    @assert Nbands   == size(V_unit_X,3)
    @assert 3        == size(V_unit_X,4)

    mask_non_diagonal_s = ones(Nqpoints, Nbands, Nbands)

    if part==:C
        for iq=1:Nqpoints
            mask_non_diagonal_s[iq,:,:] .= qpoint_weight[iq]
        end
        for ib=1:Nbands
            mask_non_diagonal_s[:,ib,ib] .= 0
        end
    elseif part==:P
        mask_non_diagonal_s = zeros(Nqpoints, Nbands, Nbands)
        for iq=1:Nqpoints
            for ib=1:Nbands
                mask_non_diagonal_s[iq,ib,ib] = qpoint_weight[iq]
            end
        end
    end

    # Aqs_plus_Aqsp
    @inline P(qs) = reshape(qs,size(qs,1),size(qs,2),1) .+ reshape(qs, size(qs,1),1,size(qs,2))
    # Aqs_minus_Aqsp
    @inline M(qs) = reshape(qs,size(qs,1),size(qs,2),1) .- reshape(qs, size(qs,1),1,size(qs,2))

    expω = exp.(ħω_unit_eV./kBT_unit_eV)
    Nbar_Nbar_p1 =  expω ./ ((expω.-1).^2)  #% occupation -> Nbar

    #% denominator
    DENOM = ((2.0.*M(ħω_unit_eV)).^2 .+ P(ħΓ_unit_eV).^2) #% [E]²

    #%  at q = 0, v = 1,2,3, ω = 0, Γ = 0, DENOM === 0 !!!
    #TODO do it carefully !  It is not simply zero ....
    Nbar_Nbar_p1[1,1:3] .= 0  
    DENOM[1,1:3,1:3] .= 1e40

    X = ( (0.5.*  P(ħω_unit_eV))
              .* (P(ħω_unit_eV.*Nbar_Nbar_p1) ./ DENOM)  
              .*  P(ħΓ_unit_eV)  ) #% [E]

    #% terms inside the summation; without the coefficient
    return [ 
        coeff_to_W_per_m_per_K * sum( 
              X .* V_unit_X[:,:,:,α] 
                .* permutedims(V_unit_X[:,:,:,β], (1,3,2)) 
                .* mask_non_diagonal_s )
        for α=1:3, β=1:3 
    ] 

end


coherent_part_of_kappa(
    qpoint_weight::Vector,
    T_unit_Kelvin::Float64,
    ω_unit_THz::Matrix,  
    Γ_unit_THz::Matrix, 
    V_unit_X::Array{ComplexF64,4}
    ) =  _part_of_kappa(
            :C,
            qpoint_weight,
            T_unit_Kelvin, 
            ω_unit_THz,  
            Γ_unit_THz, 
            V_unit_X
)


coherent_part_of_kappa(
    T_unit_Kelvin::Float64,
    ω_unit_THz::Matrix,  
    Γ_unit_THz::Matrix, 
    V_unit_X::Array{ComplexF64,4}
    ) = coherent_part_of_kappa(
            ones(size(ω_unit_THz,1)),
            T_unit_Kelvin, 
            ω_unit_THz,  
            Γ_unit_THz, 
            V_unit_X
)


peierls_part_of_kappa(
    qpoint_weight::Vector,
    T_unit_Kelvin::Float64,
    ω_unit_THz::Matrix,  
    Γ_unit_THz::Matrix, 
    V_unit_X::Array{ComplexF64,4}
    ) =  _part_of_kappa(
            :P,
            qpoint_weight,
            T_unit_Kelvin, 
            ω_unit_THz,  
            Γ_unit_THz, 
            V_unit_X
)


peierls_part_of_kappa(
    T_unit_Kelvin::Float64,
    ω_unit_THz::Matrix,  
    Γ_unit_THz::Matrix, 
    V_unit_X::Array{ComplexF64,4}
    ) = peierls_part_of_kappa(
            ones(size(ω_unit_THz,1)),
            T_unit_Kelvin, 
            ω_unit_THz,  
            Γ_unit_THz, 
            V_unit_X
)

