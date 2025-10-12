using LatticeLab
using Test

@testset "Lattice Construction" begin

    @testset "LinkInfo Creation" begin
        # Create a simple square lattice unit cell
        a = [1.0 0.0; 0.0 1.0]
        δ = [0.0 0.5; 0.0 0.0]
        m = [:A, :B]
        ξ = [[:s], [:s]]
        uc = UnitCell(2, 2, a, δ, m, ξ)

        # Define bounding box
        BBOX(m,n) = ([-0.001, -0.001], [1 0; 0 1], [m,n], [true, true])

        # Create link info by distance
        ln = LatticeLab.link_info_by_distance_direction(
            Dict(:t => (1.0, [[1.0, 0.0], [0.0, 1.0]])),
            uc;
            bounding_box = BBOX(2, 2),
            rounding_digits = 6
        )

        @test typeof(ln) <: LatticeLab.LinkInfo
        @test haskey(ln.SPNB, :t)
        @test length(ln.SPNB[:t]) > 0  # Check that there are some links
    end

    @testset "Simple Square Lattice" begin
        # Create a simple square lattice
        a = [1.0 0.0; 0.0 1.0]
        δ = reshape([0.0, 0.0], 2, 1)  # Single site at origin - 2x1 matrix for 2D lattice
        m = [:A]
        ξ = [[:s]]
        uc = UnitCell(2, 1, a, δ, m, ξ)

        # Define bounding box for 2x2 lattice
        BBOX(m,n) = ([-0.001, -0.001], [1 0; 0 1], [m,n], [true, true])

        # Create link info
        ln = LatticeLab.link_info_by_distance_direction(
            Dict(:t => (1.0, [[1.0, 0.0], [0.0, 1.0]])),
            uc;
            bounding_box = BBOX(2, 2),
            rounding_digits = 6
        )

        # Build lattice
        lattice = build_lattice(ln, BBOX(2, 2))

        @test LatticeLab.is_Lattice(lattice)
        @test LatticeLab.num_sites(lattice) == 49  # 7x7 lattice due to margin expansion during construction
        @test LatticeLab.num_sublattice(lattice) == 1
        @test LatticeLab.dimensions(lattice) == 2
        @test LatticeLab.check_compat(lattice)
    end

    @testset "Rectangular Lattice" begin
        # Create a rectangular lattice
        a = [2.0 0.0; 0.0 1.0]  # Rectangular unit cell
        δ = reshape([0.0, 0.0], 2, 1)  # Single site at origin - 2x1 matrix for 2D lattice
        m = [:A]
        ξ = [[:s]]
        uc = UnitCell(2, 1, a, δ, m, ξ)

        # Define bounding box
        BBOX(m,n) = ([-0.001, -0.001], [1 0; 0 1], [m,n], [true, true])

        # Create link info for nearest neighbors
        ln = LatticeLab.link_info_by_distance_direction(
            Dict(:tx => (2.0, [[1.0, 0.0]]), :ty => (1.0, [[0.0, 1.0]])),
            uc;
            bounding_box = BBOX(2, 3),
            rounding_digits = 6
        )

        # Build lattice
        lattice = build_lattice(ln, BBOX(2, 3))

        @test LatticeLab.is_Lattice(lattice)
        @test LatticeLab.num_sites(lattice) == 56  # 8x7 lattice due to margin expansion during construction
        @test LatticeLab.check_compat(lattice)
    end

    @testset "Two-Sublattice System" begin
        # Create a two-sublattice system (like graphene but simpler)
        a = [1.0 0.0; 0.0 1.0]
        δ = [0.0 0.5; 0.0 0.0]  # Two sublattices
        m = [:A, :B]
        ξ = [[:s], [:s]]
        uc = UnitCell(2, 2, a, δ, m, ξ)

        BBOX(m,n) = ([-0.001, -0.001], [1 0; 0 1], [m,n], [true, true])

        # Create link info
        ln = LatticeLab.link_info_by_distance_direction(
            Dict(:t => (1.0, [[0.5, 0.0]])),  # Connection between A and B
            uc;
            bounding_box = BBOX(2, 2),
            rounding_digits = 6
        )

        # Build lattice
        lattice = build_lattice(ln, BBOX(2, 2))

        @test LatticeLab.is_Lattice(lattice)
        @test LatticeLab.num_sites(lattice) == 98  # 7x7 unit cells with 2 sites each, due to margin expansion
        @test LatticeLab.num_sublattice(lattice) == 2
        @test LatticeLab.check_compat(lattice)
    end
end