
# return a list of k-path => [(kpoint_1, eigenvalues_1), (kpoint_2, eigenvalues_2), ...]
function band_structure_base(
    latt::Lattice,
    kpath::Vector{Pair{String,Vector{T}}},
    COMPUTE_EIGEN::Function,
    TEST_EIGEN::Function;
    Δk = 0.02,
    test_eigen = true,
    ) where {T<:Number}
    # walk on kpath
    bands = []
    for i ∈ 1:length(kpath)-1
        STOP  = eltype(last(kpath[i+1]))<:Real ? Float64.(last(kpath[i+1])) : ComplexF64.(last(kpath[i+1]))
        START = eltype(last(kpath[i]))<:Real   ? Float64.(last(kpath[i]))   : ComplexF64.(last(kpath[i]))
        KEY = first(kpath[i])*"->"*first(kpath[i+1])
        NK  = Int(ceil(norm(relative_reciprocal_vector(STOP.-START,latt))/Δk)) #XXX
        EIG = [ COMPUTE_EIGEN(r//NK, START, STOP) for r ∈ 0:NK-1 ]
        if test_eigen
            @assert TEST_EIGEN(EIG)
        end
        push!( bands, KEY=>EIG )
    end
    return bands
end


function band_structure_base(
    latt,
    kgrid::Array{Float64,2},
    COMPUTE_EIGEN::Function,
    TEST_EIGEN::Function;
    test_eigen = true,
    )
    # COMPUTE_EIGEN()  returns a tuple : (k,en)
    bands = hcat([COMPUTE_EIGEN(kgrid[:,i])[2] for i=1:size(kgrid,2)]...)
    @assert (!test_eigen) || TEST_EIGEN(bands)
    return bands
end
