mutable struct HeisenbergHamiltonianenR
    LATT::Lattice
    HS::Vector
    JS::Vector{ExchangeCoupling}  #
    ZS::Vector                    # on-site Zeemann
    COEFFS::Vector                # coeffs for linear combination
    G::Vector                     # classical ground state 
end
