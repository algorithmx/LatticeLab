mutable struct HeisenbergHamiltonianR{TH<:AbstractMatrix}
    LATT::Lattice
    H::TH
    J::ExchangeCoupling     #
    Z::Dict                 # on-site Zeemann
    G::Vector               # classical ground state 
end
