function stacking(LATT1::Lattice,LATT2::Lattice)
    @assert isequal_upto_orbits(LATT1,LATT2)
    LATT = copy(LATT1)
    LATT.UC.ξ = [vcat(ξ1,ξ2) for (ξ1,ξ2) ∈ zip(LATT1.UC.ξ,LATT2.UC.ξ)]
    return LATT
end

⊕(LATT1::Lattice,LATT2::Lattice) = stacking(LATT1,LATT2)

function stacking(H1::HoppingHamiltonianQ,H2::HoppingHamiltonianQ)
    @assert isequal_upto_orbits(H1.LATT, H2.LATT)
    LATT = (H1.LATT ⊕ H2.LATT)
    (m1,n1) = size(first(values(H1.MAT)))
    (m2,n2) = size(first(values(H2.MAT)))
    MAT1 = Dict(k=>rightpad(sparse(v),m2,n2) for (k,v) ∈ H1.MAT)
    MAT2 = Dict(k=>leftpad(m1,n1,sparse(v)) for (k,v) ∈ H2.MAT)
    return HoppingHamiltonianQ(LATT,MAT1,H1.T,H1.V)+HoppingHamiltonianQ(LATT,MAT2,H2.T,H2.V)
end

⊕(H1::HoppingHamiltonianQ,H2::HoppingHamiltonianQ) = stacking(H1,H2)

function stacking(H1::HoppingHamiltonianenQ,H2::HoppingHamiltonianenQ)
    @assert isequal_upto_orbits(H1.LATT, H2.LATT)
    LATT = (H1.LATT ⊕ H2.LATT)
    (m1,n1) = size(first(values(H1.MATS[1])))
    (m2,n2) = size(first(values(H2.MATS[1])))
    MATS1 = [Dict(k=>rightpad(sparse(v),m2,n2) for (k,v) ∈ MAT) for MAT ∈ H1.MATS]
    MATS2 = [Dict(k=>leftpad(m1,n1,sparse(v))  for (k,v) ∈ MAT) for MAT ∈ H2.MATS]
    return HoppingHamiltonianenQ(LATT,MATS1,H1.TS,H1.COEFFS,H1.VS)+HoppingHamiltonianenQ(LATT,MATS2,H2.TS,H2.COEFFS,H2.VS)
end

⊕(H1::HoppingHamiltonianenQ,H2::HoppingHamiltonianenQ) = stacking(H1,H2)


function stacking(H1::SuperconductingGapFunctionQ,H2::SuperconductingGapFunctionQ)
    @assert isequal_upto_orbits(H1.LATT, H2.LATT)
    LATT = (H1.LATT ⊕ H2.LATT)
    (m1,n1) = size(first(values(H1.MAT)))
    (m2,n2) = size(first(values(H2.MAT)))
    MAT1 = Dict(k=>rightpad(sparse(v),m2,n2) for (k,v) ∈ H1.MAT)
    MAT2 = Dict(k=>leftpad(m1,n1,sparse(v)) for (k,v) ∈ H2.MAT)
    return SuperconductingGapFunctionQ(LATT,MAT1)+SuperconductingGapFunctionQ(LATT,MAT2)
end

⊕(H1::SuperconductingGapFunctionQ,H2::SuperconductingGapFunctionQ) = stacking(H1,H2)

function stacking(H1::SuperconductingGapFunctionenQ,H2::SuperconductingGapFunctionenQ)
    @assert isequal_upto_orbits(H1.LATT, H2.LATT)
    LATT = (H1.LATT ⊕ H2.LATT)
    (m1,n1) = size(first(values(H1.MATS[1])))
    (m2,n2) = size(first(values(H2.MATS[1])))
    MATS1 = [Dict(k=>rightpad(sparse(v),m2,n2) for (k,v) ∈ MAT) for MAT ∈ H1.MATS]
    MATS2 = [Dict(k=>leftpad(m1,n1,sparse(v))  for (k,v) ∈ MAT) for MAT ∈ H2.MATS]
    return SuperconductingGapFunctionenQ(LATT,MATS1,H1.COEFFS)+SuperconductingGapFunctionenQ(LATT,MATS2,H2.COEFFS)
end

⊕(H1::SuperconductingGapFunctionenQ,H2::SuperconductingGapFunctionenQ) = stacking(H1,H2)

function stacking(BdG1::BdGHamiltonianQ,BdG2::BdGHamiltonianQ)
    @assert isequal_upto_orbits(BdG1.H.LATT, BdG2.H.LATT)
    @assert BdG1.iσy == BdG2.iσy
    LATT = (BdG1.H.LATT ⊕ BdG2.H.LATT)
    Tiσy = typeof(BdG1.iσy)
    return BdGHamiltonianenQ{Tiσy}((BdG1.H ⊕ BdG2.H), (BdG1.Δ ⊕ BdG2.Δ), (BdG1.iσy ⊕ BdG2.iσy))
end


⊕(BdG1::BdGHamiltonianQ,BdG2::BdGHamiltonianQ) = stacking(BdG1,BdG2)


function stacking(BdG1::BdGHamiltonianenQ,BdG2::BdGHamiltonianenQ)
    @assert isequal_upto_orbits(BdG1.HS.LATT, BdG2.HS.LATT)
    @assert BdG1.iσy == BdG2.iσy
    LATT = (BdG1.HS.LATT ⊕ BdG2.HS.LATT)
    Tiσy = typeof(BdG1.iσy)
    return BdGHamiltonianenQ{Tiσy}((BdG1.HS ⊕ BdG2.HS), (BdG1.ΔS ⊕ BdG2.ΔS), (BdG1.iσy ⊕ BdG2.iσy))
end


⊕(BdG1::BdGHamiltonianenQ,BdG2::BdGHamiltonianenQ) = stacking(BdG1,BdG2)
