using LatticeLab
d1 = [ sqrt(3)/2, 1/2]
d2 = [-sqrt(3)/2, 1/2]
d3 = [         0,  -1] ;

BBOX(m,n) = ([-0.001rand(),-0.001rand()], # origin
             [1 0 ; 0 1],                 # supercell basis
             [m,n],                       # supercell shifts
             [true,true])                 # P.B.C. conditions

UC_Honeycomb = LatticeLab.UnitCell(
    2,                  # dimension
    2,                  # number of sublattices
    [sqrt(3) sqrt(3)/2 
     0             3/2],# Bravais basis
    [0       sqrt(3)/2
     0             1/2],# sublattice coordinates  
    [ :A,    :B  ],     # sublattice symbols
    [[:pz], [:pz]]      # orbits 
) ;

UC_LN = LatticeLab.link_info_by_distance_direction(
    Dict( :tr => (1,[d1,-d1]),
          :tg => (1,[d2,-d2]),
          :tb => (1,[d3,-d3]), ), # symbol => (distance, (directions))
    UC_Honeycomb;
    bounding_box = BBOX(1,1),
    rounding_digits=12
) ;

Graphene = build_lattice(UC_LN,BBOX(3,3)) ;

VSTYLE = Dict(
    :A=>("black",6.0,:dot),
    :B=>("black",6.0,:circle)
)
ESTYLE = Dict(
    :tr=>("red",  3.0,:solid),
    :tg=>("green",3.0,:solid),
    :tb=>("blue", 3.0,:solid),
)

display("image/svg+xml", show_lattice_svg(Graphene, VSTYLE, ESTYLE; upscale=60))    