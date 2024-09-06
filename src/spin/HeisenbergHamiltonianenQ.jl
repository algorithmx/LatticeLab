#TODO  implement  https://arxiv.org/pdf/1012.4166.pdf

mutable struct HeisenbergHamiltonianenQ
    LATT::Lattice
    MATS::Vector                  # mn => M
    JS::Vector{ExchangeCoupling}  #
    ZS::Vector                    # on-site Zeemann
    COEFFS::Vector                # coeffs for linear combination
end

