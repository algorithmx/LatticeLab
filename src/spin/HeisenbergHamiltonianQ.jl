#TODO  implement  https://arxiv.org/pdf/1012.4166.pdf

mutable struct HeisenbergHamiltonianQ
    LATT::Lattice
    MAT::Dict               # mn => M
    J::ExchangeCoupling     #
    Z::Dict                 # on-site Zeemann
end

