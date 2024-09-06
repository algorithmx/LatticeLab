mutable struct EigenModes
    ω::Vector{ComplexF64}
    u::Modes
    function EigenModes( ω::Vector{Float64}, u::Modes )
        return new( ω.+(+0.0+0.0im), u )
    end
    function EigenModes( ω::Vector{ComplexF64}, u::Modes )
        return new( ω, u )
    end
    function EigenModes( Nmodes::Int64 )
        return new(zeros(ComplexF64,Nmodes), zeros(ComplexF64,Nmodes,Nmodes))
    end
end


@inline is_EigenModes(em) = (  Set(fieldnames(EigenModes))==Set(fieldnames(typeof(em)))
                            && isa(em.ω, Vector) && eltype(em.ω)<:Number
                            && is_Modes(em.u) )


check_compat(e0::EigenModes) = ( length(e0.ω)==size(e0.u,2) )


copy(eig::EigenModes) = EigenModes(eig.ω, eig.u)


function DoS(En::Vector{Float64}, η; N=1000)
    min_En = minimum(En)
    max_En = maximum(En)
    d = 1.05*(max_En-min_En)/N
    ω = (min_En - 0.025*(max_En-min_En) - im*η) .+ d.*collect(1:N)
    dos = Float64[ imag(sum(π./(en0.-En))) for en0 ∈ ω ]
    return ω, dos
end


@inline convert(::Type{T}, x::T) where {T<:EigenModes} = x


function convert(::Type{T}, em) where {T<:EigenModes}
    @assert is_EigenModes(em)
    @assert length(em.ω)==size(em.u,2) # check_compat(em)
    EigenModes(Vector{ComplexF64}(em.ω),Modes(em.u))
end


function solve_eigenmodes(
    mat::T;
    EPS = 1e-10,
    eigen_test = (x->true),
    eigen_process = (x->x)
    ) where { T<:AbstractMatrix }

    M  = convert( SparseMatrixCSC{ComplexF64,Int64}, mat )
    nM = norm( (M') .- M )
    @assert nM < EPS  "norm( (M') .- M )  =  $nM"

    print("[DEBUG] norm((M').-M) = "); println(nM); 
    flush(stdout);

    println("[DEBUG] eigen(...) costs:  ")

    @time sol = eigen(LinearAlgebra.Matrix(M))

    flush(stdout);
    print("[DEBUG] Max(imag(eigenvalues)) = ");
    println(maximum(abs.(imag.(sol.values))));  
    flush(stdout);

    @assert eigen_test(sol.values)
    ASS = sortperm(real.(sol.values))

    return EigenModes( map(eigen_process, sol.values[ASS]), sol.vectors[:,ASS] )
end
