function generate_kmesh(dk::Float64, latt::Lattice)
    b = reciprocal_basis(latt)
    invb = inv(b)
    inside_BZ(point) = all(x->(x<1.0 && x>=0.0),(invb*point))
    min_max = map(i->(minimum(b[:,i]),maximum(b[:,i])), collect(1:size(b,2)))
    return [ collect(p)
                for p ∈ Iterators.product((collect((bmin-dk):dk:(bmax+dk)) for (bmin,bmax) ∈ min_max)...)
                    if inside_BZ(collect(p)) ]
end


@inline function EnS0(
    kmsh::Vector,
    ff::Function;
    POST = (x->x),
    TEST = (ω->true)
    )::Vector{Float64}
    ω = vcat(map(ff,kmsh)...)
    @assert TEST(ω)
    return sort(POST.(ω))
end


function EnS(
    dk::Float64,
    DQ::Union{DynamicalMatrixQ,DynamicalMatricesQ};
    EPS = 1e-10,
    test_eigen = true
    )::Vector{Float64}
    kmsh = generate_kmesh(dk, DQ.LATT)
    ff(p) = (eigvals(LinearAlgebra.Matrix(DynamicalMatrix(p,DQ))) .+ 0.0im)
    EnS0( kmsh, ff, POST=(x->real(sqrt(x))),
                    TEST=(w->((!test_eigen) || (all(abs.(imag.(w)).<EPS) && all(real.(w).>-EPS)))) )
end


function EnS(
    dk::Float64,
    HQ::Union{HoppingHamiltonianQ,HoppingHamiltonianenQ};
    EPS = 1e-10,
    test_eigen = true
    )::Vector{Float64}
    kmsh = generate_kmesh(dk, HQ.LATT)
    ff(p) = (eigvals(Matrix(HoppingHamiltonian(p,HQ))) .+ 0.0im)
    EnS0( kmsh, ff, POST=(x->real(x)),
                    TEST=(w->((!test_eigen) || (w->(all(abs.(imag.(w)).<EPS))))) )
end

