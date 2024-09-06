using LatticeLab
a1 = [1, 0]
a2 = [1/2, sqrt(3)/2] ;
BBOX(m,n) = ([-0.001,-0.001], # origin
             [1 0 ; 0 1],                 # supercell basis
             [m,n],                       # supercell shifts
             [true,true])                 # P.B.C. conditions

aa = hcat( 2 .*a1, 2 .*a2 ) |> Coordinates
dd = hcat( [0,0], a1, a2) |> Coordinates
mm = [ :A,    :B,    :C     ] |> Masses
xx = [ [:pz,],[:pz,],[:pz,] ] .|> Orbits

UC_Kagome = UnitCell(2, 3, aa, dd, mm, xx) ;

LN_Kagome = LatticeLab.link_info_by_distance_direction(
    Dict(:t => (1,[]),), # symbol => (distance, (directions))
    UC_Kagome;
    bounding_box = BBOX(1,1),
    rounding_digits=6
)

Kagome = build_lattice(LN_Kagome,BBOX(3,3)) ;

VSTYLE = Dict(
    :A=>("red",  6.0,:dot),
    :B=>("green",6.0,:dot),
    :C=>("blue", 6.0,:dot)
)
ESTYLE = Dict(
    :t=>("black", 3.0,:solid),
)
display("image/svg+xml", 
    show_lattice_svg(Kagome, VSTYLE, ESTYLE; upscale=60))