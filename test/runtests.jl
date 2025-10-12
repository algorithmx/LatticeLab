using LatticeLab
using Test

println("Running LatticeLab test suite...")

# Include all test files
include("test_basic_types.jl")
include("test_lattice_construction.jl")
include("test_hamiltonians.jl")

println("All tests completed!")