using LatticeLab
using Test

@testset "Basic Types" begin

    @testset "UnitCell" begin
        # Test creating a simple 2D square lattice unit cell
        a = [1.0 0.0; 0.0 1.0]  # Square lattice vectors
        δ = [0.0 0.5; 0.0 0.0]  # Two sublattices: (0,0) and (0.5,0)
        m = [:A, :B]  # Atom labels
        ξ = [[:s], [:s]]  # One s-orbital per sublattice

        uc = UnitCell(2, 2, a, δ, m, ξ)

        @test uc.dim == 2
        @test uc.nsubl == 2
        @test size(uc.a) == (2, 2)
        @test size(uc.δ) == (2, 2)
        @test length(uc.m) == 2
        @test length(uc.ξ) == 2
        @test uc.m == [:A, :B]
    end

    @testset "Coordinates" begin
        # Test Coordinates type
        coords = Coordinates([1.0 2.0; 3.0 4.0])
        @test size(coords) == (2, 2)
        @test coords[1,1] == 1.0
        @test coords[2,2] == 4.0
    end

    @testset "Masses" begin
        # Test Masses type
        masses = Masses([:A, :B, :C])
        @test length(masses) == 3
        @test masses[1] == :A
        @test masses[3] == :C
    end

    @testset "Orbits" begin
        # Test Orbits type - should be a flat Vector{Symbol}
        orbits = Orbits([:s, :px, :s])
        @test length(orbits) == 3
        @test orbits[1] == :s
        @test orbits[2] == :px
        @test orbits[3] == :s
    end

    @testset "BoundingBox" begin
        # Test BoundingBox type
        origin = [0.0, 0.0]
        trans = [1 0; 0 1]  # Matrix{Int64} as required by BoundingBox type
        units = [3, 3]
        pbc = [true, true]

        bbox = BoundingBox((origin, trans, units, pbc))

        @test bbox[1] == origin
        @test bbox[2] == trans
        @test bbox[3] == units
        @test bbox[4] == pbc
    end
end