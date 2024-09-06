using Pkg; Pkg.activate("../");
using LatticeLab
Pkg.activate("../../BandStructures");
using BandStructures
using LinearAlgebra
using PyPlot

## =========================================================================

sigma0 = [1   0; 0  1]
sigmax = [0   1; 1  0]
sigmay = [0 -im; im 0]
sigmaz = [1   0; 0 -1]
I2 = sigma0 ;

## =========================================================================

BBOX(m,n) = ([-0.00021,-0.0001], # origin
             I2,                 # supercell basis
             [m,n],              # supercell shifts
             [true,true])        # P.B.C. conditions

UC_Sqlatt = LatticeLab.UnitCell(
    2,              # dimension
    1,              # number of sublattices
    I2,             # Bravais basis
    zeros(2,1),     # sublattice coordinates  
    [:A, ],         # sublattice symbols
    [[:pzu, :pzd],] # orbits 
) ;
@assert check_compat(UC_Sqlatt)
LN_Sqlatt = LatticeLab.link_info_by_distance_direction(
    Dict(:txp => (1,[[ 1,0]]),:typ => (1,[[0, 1]]),
         :txm => (1,[[-1,0]]),:tym => (1,[[0,-1]]), ),
    UC_Sqlatt;
    bounding_box = BBOX(1,1),
    rounding_digits=12
) ;
@assert check_compat(LN_Sqlatt)

Sqlatt = build_lattice(LN_Sqlatt,BBOX(3,3)) ;

VSTYLE = Dict(:A=>("black",6.0,:dot))
ESTYLE = Dict(
    :txp=>("red",  1.5,:solid),
    :txm=>("red",  3.5,:dashed),
    :typ=>("blue", 1.5,:solid),
    :tym=>("blue", 3.5,:dashed),
)
display("image/svg+xml", 
    show_lattice_svg(Sqlatt, VSTYLE, ESTYLE; upscale=60))

## =========================================================================

# the complicated way of assigning hopping parameters
func1() = (:MAT,(sigma0,))
func(x) = (:MAT,(x,     ))
hBHZF  = Dict(:txp => (a,b)->func((-im*a/2).*sigmax.+(b).*sigmaz),
              :typ => (a,b)->func((-im*a/2).*sigmay.+(b).*sigmaz),
              :txm => (a,b)->func(((-im*a/2).*sigmax.+(b).*sigmaz)'),
              :tym => (a,b)->func(((-im*a/2).*sigmay.+(b).*sigmaz)'),)
hBHZP  = Dict(:txp => (:tA, :tB), :typ => (:tA, :tB),
              :txm => (:tA, :tB), :tym => (:tA, :tB),)
@inline hBHZ(A,B) = LatticeLab.dispatch_params(
    Dict(:tA=>A,:tB=>B), 
    hBHZF, hBHZP, 0.0 )
@inline vBHZ(Z,B) = Dict(:A=>((Z-4B).*sigmaz))

## =========================================================================

b = reciprocal_basis(Sqlatt) ;
kpath = [ "Γ" => (0   ).*b[:,1],
          "X" => (1//2).*b[:,1],
          "M" => (1//2).*b[:,1] .+ (1//2).*b[:,2],
          "Γ" => (0   ).*b[:,1], 
          "Y" => (0   ).*b[:,1] .+ (1//2).*b[:,2],] ;

## =========================================================================

BBOX_open_x(m,n) = BoundingBox((
             [-0.00021,-0.0001], # origin
             I2,                 # supercell basis
             [m,n],              # supercell shifts
             [false,true]))      # P.B.C. conditions
BBOX_open_y(m,n) = BoundingBox((
             [-0.00021,-0.0001], # origin
             I2,                 # supercell basis
             [m,n],              # supercell shifts
             [true,false]))      # P.B.C. conditions

LATT00_y_SqLatt = build_lattice(LN_Sqlatt, BBOX(1,128)) ;
LATTy_open = enlarge(LATT00_y_SqLatt, BBOX_open_y(3,1)) ;

b_y = reciprocal_basis(LATTy_open)
kpath_y = [ "Γ" => (0.0).*b_y[:,1],
            "Γ " => (0.5).*b_y[:,1], ] ;

## =========================================================================

B = 1.0
HQy_open = [kspace_hopping_hamiltonian( 
                HoppingParameter(LATTy_open.UC, hBHZ(1.0, B)), 
                vBHZ(r*B, B), 
                LATTy_open)
            for r = [0.2,0.8,2.0,3.2,4.0] ] ;
BSy_open = [LatticeLab.band_structure(
                kpath_y, HQ, Δk=1e-3
                ) |> BandStructures.LatticeLab_bands_BandStructure
            for HQ in HQy_open ] ;

plot_bands("SquareLattice_BHZ_B=1.0_ZBratio=(0.2,0.8,2.0,3.2,4.0).png", 
           BSy_open,
           dpi   = 800,
           settings = Dict(
                :line_colors=>["red", "blue", "green", "cyan", "violet"],
                :range => (-1.7,1.7),
                :K_sep=>-1,
                :lw=>0.8,
                :aspect_ratio=>0.7,
                :figure_size=>(10,6),
           )
)

## =========================================================================

B = 0.5
HQy_open = [kspace_hopping_hamiltonian( 
                HoppingParameter(LATTy_open.UC, hBHZ(1.0, B)), 
                vBHZ(r*B, B), 
                LATTy_open)
            for r = [0.2,0.8,2.0,3.2,4.0] ] ;
BSy_open = [LatticeLab.band_structure(
                kpath_y, HQ, Δk=1e-3
                ) |> BandStructures.LatticeLab_bands_BandStructure
            for HQ in HQy_open ] ;

plot_bands("SquareLattice_BHZ_B=0.5_ZBratio=(0.2,0.8,2.0,3.2,4.0).png", 
           BSy_open,
           dpi   = 800,
           settings = Dict(
                :line_colors=>["red", "blue", "green", "cyan", "violet"],
                :range => (-1.7,1.7),
                :K_sep=>-1,
                :lw=>0.8,
                :aspect_ratio=>0.7,
                :figure_size=>(10,6),
           )
)

## =========================================================================

Z = 0.0
HQs = [ kspace_hopping_hamiltonian( 
            HoppingParameter(Sqlatt.UC, hBHZ(1.0,B)), 
            vBHZ(Z,B), 
            Sqlatt )
        for B = 0:0.01:0.04 ] ;
BSs = [ LatticeLab.band_structure(
            kpath, HQ, Δk=5e-4
            ) |> BandStructures.LatticeLab_bands_BandStructure
        for HQ in HQs ] ;
plot_bands("SquareLattice_BHZ_Z=0.0.png", 
           BSs,
           dpi  = 800,
           settings = Dict(
                :line_colors=>["black", "red", "blue", "green", "orange"],
                :lw=>0.8,
                :K_sep=>-1,
                :aspect_ratio=>1,
                :figure_size=>(12,6),)
)


## =========================================================================

Z = 0.17

HQs = [ kspace_hopping_hamiltonian( 
            HoppingParameter(Sqlatt.UC, hBHZ(1.0,B)), 
            vBHZ(Z,B), 
            Sqlatt )
        for B = 0:0.01:0.04 ] ;
BSs = [ LatticeLab.band_structure(
            kpath, HQ, Δk=5e-4
            ) |> BandStructures.LatticeLab_bands_BandStructure
        for HQ in HQs ] ;
plot_bands("SquareLattice_BHZ_Z=0.17.png", 
           BSs,
           dpi  = 800,
           settings = Dict(
                :line_colors=>["black", "red", "blue", "green", "orange"],
                :lw=>0.8,
                :K_sep=>-1,
                :aspect_ratio=>1,
                :figure_size=>(12,6),)
)