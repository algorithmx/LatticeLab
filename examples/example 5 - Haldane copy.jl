using Pkg; Pkg.activate("../");
using LatticeLab
Pkg.activate("../../BandStructures");
using BandStructures
using LinearAlgebra
using PyPlot

a1 = [1,   0        ]
a2 = [1/2, sqrt(3)/2]
α = sqrt(3)
Atr = [a1, a2, a2.-a1, -a1, -a2, a1.-a2]
d1 = (α/3).*(Atr[1].+Atr[2])
d2 = (α/3).*(Atr[3].+Atr[4])
d3 = (α/3).*(Atr[5].+Atr[6])
delta = [d1,d2,d3]

BBOX(m,n) = ([-0.00021,-0.0001], # origin
             [1 0 ; 0 1],                 # supercell basis
             [m,n],                       # supercell shifts
             [true,true])                 # P.B.C. conditions

aa = hcat( α.*a1, α.*a2) |> LatticeLab.Coordinates
dd = hcat([0, 0], d1   ) |> LatticeLab.Coordinates
mm = [ :A,       :B      ] |> LatticeLab.Masses
xx = [[:pu,:pd],[:pu,:pd]] .|> LatticeLab.Orbits

UC_Honeycomb = LatticeLab.UnitCell( 2, 2, aa, dd, mm, xx )

nb1a(k) = Spring( (1,2) => [ delta[k],] )
nb1b(k) = Spring( (2,1) => [-delta[k],] )
nb2_cwise_s1  = Spring((1,1) => [α.*Atr[2], α.*Atr[4], α.*Atr[6]])
nb2_cwise_s2  = Spring((2,2) => [α.*Atr[2], α.*Atr[4], α.*Atr[6]])
nb2_ccwise_s1 = Spring((1,1) => [α.*Atr[1], α.*Atr[3], α.*Atr[5]])
nb2_ccwise_s2 = Spring((2,2) => [α.*Atr[1], α.*Atr[3], α.*Atr[5]])
NB_Haldane = [ nb1a(1), nb1a(2), nb1a(3), nb1b(1), nb1b(2), nb1b(3), 
               nb2_cwise_s1, nb2_cwise_s2, nb2_ccwise_s1, nb2_ccwise_s2 ]
SP_Haldane = [ :t1Ar,   :t1Ag,   :t1Ab,   :t1Br,   :t1Bg,   :t1Bb,
               :t2ps1,  :t2ps2,  :t2ms1,  :t2ms2  ]
LN_Haldane = LatticeLab.LinkInfo( 
    UC_Honeycomb, 
    LatticeLab.zict(SP_Haldane, NB_Haldane) )
@assert check_compat(LN_Haldane)

## ===============================================================

VSTYLE = Dict(
    :A=>("black",6.0,:dot),
    :B=>("black",6.0,:circle)
)
ESTYLE = Dict(
    :t1Ar => ("red",    2.0, :solid ),
    :t1Ag => ("green",  2.0, :solid ),
    :t1Ab => ("blue",   2.0, :solid ),
    :t1Br => ("red",    5.0, :dashed),
    :t1Bg => ("green",  5.0, :dashed),
    :t1Bb => ("blue",   5.0, :dashed),
    :t2ps1=> ("violet", 2.0, :solid ),
    :t2ps2=> ("pink",   2.0, :solid ),
    :t2ms1=> ("violet", 5.0, :dashed),
    :t2ms2=> ("pink",   5.0, :dashed),
)
display("image/svg+xml", 
    show_lattice_svg(Honeycomb_Haldane, VSTYLE, ESTYLE; upscale=60))

## ===============================================================

# the complicated way of assigning hopping parameters
func1() = (:ALL,(1,))
func(x) = (:ALL,(x,))
hHaldaneF  = Dict( :t1Ar=>func1, :t1Ag=>func1, :t1Ab=>func1,
                   :t1Br=>func1, :t1Bg=>func1, :t1Bb=>func1,
                   # the Haldanish lines
                   :t2ps1=>(x,y)->func(x*cis(y)),
                   :t2ps2=>(x,y)->func(x*cis(-y)),
                   :t2ms1=>(x,y)->func(x*cis(-y)),
                   :t2ms2=>(x,y)->func(x*cis(y)), )
hHaldaneP  = Dict( :t1Ar=>(), :t1Ag=>(), :t1Ab=>(),
                   :t1Br=>(), :t1Bg=>(), :t1Bb=>(),
                   :t2ps1=>(:tnn,:phi),:t2ms1=>(:tnn,:phi),
                   :t2ps2=>(:tnn,:phi),:t2ms2=>(:tnn,:phi), )
@inline hHaldane(tnn,phi) = LatticeLab.dispatch_params(
    Dict(:tnn=>tnn,:phi=>phi), 
    hHaldaneF, 
    hHaldaneP, 
    0.0
)

## ===============================================================


b = reciprocal_basis(Honeycomb_Haldane) ;
kpath = [ "Γ" => (0   ).*b[:,1],
          "M" => (1//2).*b[:,1] .+ (0//2).*b[:,2],
          "K" => (2//3).*b[:,1] .+ (1//3).*b[:,2],
          "Γ" => (0   ).*b[:,1], ] ;
HQs = [ kspace_hopping_hamiltonian( 
            HoppingParameter(Honeycomb_Haldane.UC, hHaldane(0.3,phi)), 
            zero_onsite_potential(Honeycomb_Haldane), 
            Honeycomb_Haldane )
        for phi = 0:0.1:0.4 ] ;
BSs = [ LatticeLab.band_structure(
            kpath, HQ, Δk=1e-3
            ) |> BandStructures.LatticeLab_bands_BandStructure
        for HQ in HQs ] ;

## ===============================================================

plot_bands("Haldane_tnn=0.3.png", 
           BSs,
           dpi  = 1600,
           settings = Dict(
                :line_colors=>["black", "red", "blue", "green", "orange"],
                :lw=>0.8,
                :K_sep=>-1,
                :aspect_ratio=>0.6,
                :figure_size=>(16,20),)
)

## ===============================================================

t2 = 0.3
K1 = (2//3).*b[:,1] .+ (1//3).*b[:,2]
K2 = (1//3).*b[:,1] .+ (2//3).*b[:,2]
function is_gapless(H,eps=1e-3)
    en1 = eigen(Matrix(HoppingHamiltonian(K1,H))).values
    en2 = eigen(Matrix(HoppingHamiltonian(K2,H))).values
    (  abs(maximum(en1)-minimum(en1))<eps 
    || abs(maximum(en2)-minimum(en2))<eps  )
end
Ham(t2,M,phi) = 
    kspace_hopping_hamiltonian( 
        HoppingParameter(Honeycomb_Haldane.UC, hHaldane(t2,phi)), 
        Dict(:A=>M,:B=>-M),
        Honeycomb_Haldane )
@time Phase = [ is_gapless(Ham(t2,t2*r,phi),4e-2) 
                for r=-6:0.05:6, phi=-2pi:(pi/50):2pi ] ;

matshow(Phase)
xlabel("phi : -2π ... 2π")
ylabel("M/t2 : -6 ... 6")
title("gapless at K or K', t2=0.3")
savefig("Haldane_gapless.png")

## ===============================================================

BBOX_ENL(m,n) = BoundingBox(
    ([-0.00021,-0.0001],  # origin
     [[1,0] [-1,2]],      # supercell basis
     [m,n],               # supercell shifts
     [true,true])         # P.B.C. conditions
)
BBOX_stack_stripe_along_x(Len) = BoundingBox(
    ([-0.00021,-0.0001],  # origin
     [[1,0] [0,1]],       # supercell basis
     [Len,1],             # supercell shifts
     [true,false])        # P.B.C. conditions
)
LATT00_y_Honeycomb_ENL = build_lattice(LN_Haldane, BBOX_ENL(1,16)) ;
LATTy_open = enlarge(LATT00_y_Honeycomb_ENL, BBOX_stack_stripe_along_x(3)) ;
b_y = reciprocal_basis(LATTy_open) ;
kpath_y = [ "Γ" => (0).*b_y[:,1],
            "Γ" => (1).*b_y[:,1], ] ;
HQy_open = [kspace_hopping_hamiltonian( 
                HoppingParameter(LATTy_open.UC, hHaldane(0.1,0.2)), 
                Dict(:A=>M,:B=>-M), 
                LATTy_open)
            for M = 0.06:0.02:0.14 ] ;
BSy_open = [LatticeLab.band_structure(
                kpath_y, HQ, Δk=2e-4
                ) |> BandStructures.LatticeLab_bands_BandStructure
            for HQ in HQy_open ] ;

## ===============================================================

plot_bands("Haldane_t2=0.1_phi=0.2_M=(0.06,0.02,0.14).png", 
           BSy_open,
           dpi   = 1600,
           settings = Dict(
                :line_colors=>["black", "red", "blue", "green", "orange"],
                :K_sep=>-1,
                :lw=>0.8,
                :aspect_ratio=>0.4,
                :figure_size=>(16,24),)
)


## ===============================================================

plot_bands("Haldane_t2=0.1_phi=0.2_M=(0.06,0.02,0.14)__ZOOM.png", 
           BSy_open,
           dpi   = 1600,
           settings = Dict(
                :line_colors=>["black", "red", "blue", "green", "orange"],
                :K_sep => -1,
                :range => (-0.55,-0.05),
                :lw    => 0.8,
                :aspect_ratio => 3,
                :figure_size  => (12,16),)
)