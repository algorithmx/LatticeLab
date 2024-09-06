ENV["PSEUDO_ROOT"] = "/home/dabajabaza/abinitio/pseudopotentials"
include("../../Biphenylene/include/using.jl")
include("../../Biphenylene/include/tools.jl")
include("../../Biphenylene/include/CIF2LN.jl")
include("../../Biphenylene/include/wannier90_utils.jl")

Bi2Pd_hr_fn  = "Bi2Pd_soc_hr.dat"
Bi2Pd_cif_fn = "Bi2Pd.cif"

## ===========================================================

Porbits = [:pzu,:pzd,:pxu,:pxd,:pyu,:pyd]
Dorbits = [:dz2u,:dz2d,:dzxu,:dzxd,:dzyu,:dzyd,:dx2y2u,:dx2y2d,:dxyu,:dxyd] ;
P_u     = diagm(0=>vcat([[1,0] for i=1:11]...))
P_d     = diagm(0=>vcat([[0,1] for i=1:11]...))
P_Bi    = diagm(0=>vcat([[1 for i=1:12], [0 for i=1:10]]...))
P_Pd    = diagm(0=>vcat([[0 for i=1:12], [1 for i=1:10]]...))
P_Bi_l0 = diagm(0=>vcat([[1,1,0,0,0,0,1,1,0,0,0,0], [0 for i=13:22]]...))
P_Bi_l1 = P_Bi .- P_Bi_l0
P_Pd_l0 = diagm(0=>vcat([[0 for i=1:12], [1,1], [0 for i=15:22]]...))
P_Pd_l1 = diagm(0=>vcat([[0 for i=1:14], [1,1,1,1], [0 for i=19:22]]...))
P_Pd_l2 = P_Pd .- P_Pd_l0 .- P_Pd_l1 
op(OP) = (k,v)->abs(v'*OP*v)

## ===========================================================

function Bi2Pd(fn)
    ORIG = [-8e-5, -1e-4, -7e-5]
    ID3  = [1 0 0; 0 1 0; 0 0 1]
    @inline select_non_empty(D) = Dict(k=>v for (k,v) in D if length(v)>0)
    xx = [Porbits, Porbits, Dorbits] .|> LatticeLab.Orbits
    UC = CIF2UC(CIF(fn), xx)
    @time begin
        ALL_LN     = all_links_between_sublattices( 
                       UC; 
                       maxlen=120, max_distance=200, large_enough_margin=20, 
                       rounding_digits=6 ) |> select_non_empty
        ALL_LN_nnn = pick_nearest_by_ij(
                       ALL_LN, 
                       merge( Dict((i,i)=>15 for i=1:3), 
                              Dict((i,j)=>15 for i=1:3 for j=1:3 if i!=j) ); 
                       digits=6
        )
        # TODO needs finer constructor, with x- and y-directions
        SP = all_links_dict_to_sp_no_directions(ALL_LN_nnn) ;
    end
    LN = LinkInfo(copy(UC), SP) ;
    return build_lattice(LN, BoundingBox((ORIG,ID3,[3,3,1],[true,true,false]))) ;
end

Bi2Pd_LATT = Bi2Pd(Bi2Pd_cif_fn) ;

## ===========================================================

const KPATH_REL = [ 
    "Γ"  =>  [0,    0,   0 ],
    "X"  =>  [1//2, 0,   0 ],
    "M"  =>  [1//2, 1//2,0 ],
    "Γ"  =>  [0,    0,   0 ], 
    "Y"  =>  [0, 1//2,   0 ], ]
@inline kvec_rel2abs(rel,  b) = [k=>(b*vrel) for (k,vrel) in rel]
@inline HSP(b)  = kvec_rel2abs(KPATH_REL,b)

## ===========================================================

Bi2Pd_HQ_Wann = make_HWannQ(
        readlines(Bi2Pd_hr_fn), 
        Bi2Pd_LATT; 
        thr=0.001
) ;

BS_Wann = LatticeLab.band_structure_with_markers(
                HSP(reciprocal_basis(Bi2Pd_LATT)),
                Bi2Pd_HQ_Wann,
                [op(P_Bi_l0),op(P_Bi_l1),op(P_Pd_l0),op(P_Pd_l1),op(P_Pd_l2)],
                Δk=1e-2
) |> BandStructures.LatticeLab_bands_BandStructure ;

## ===========================================================

plot_bands( "Bi2Pd_Wannier_Pd_lz=0(green),1(orange),2(violet).pdf",
            [BS_Wann,];
            SIZE  = [0.0,   0.0,    90.0,    90.0,     90.0     ],
            COLOR = ["red", "blue", "green", "orange", "violet" ],
            settings = Dict(
                  :line_colors=>["black"],
                  :lw=>[0.4],
                  :ref_levels=>[(4.0201,"gray"),(0.0,"gray")],
                  :range=>(2,6),   
                  :K_sep=>-1,
                  :markerstrokealpha=>0.7,
                  :markerstrokewidth=>0.3,
                  :aspect_ratio=>0.6,
                  :figure_size=>(8,16),)
)

plot_bands( "Bi2Pd_Wannier_Bi_lz=0(red),1(blue).pdf",
            [BS_Wann,];
            SIZE  = [90.0,  90.0,   0.0,     0.0,       0.0     ],
            COLOR = ["red", "blue", "green", "orange", "violet" ],
            settings = Dict(
                  :line_colors=>["black"],
                  :lw=>[0.4],
                  :ref_levels=>[(4.0201,"gray"),(0.0,"gray")],
                  :range=>(2,6),   
                  :K_sep=>-1,
                  :markerstrokealpha=>0.7,
                  :markerstrokewidth=>0.3,
                  :aspect_ratio=>0.6,
                  :figure_size=>(8,16),)
)

## ===========================================================

function SP2T_raw(HQWannier; thr=0.02)
    @inline trim_empty(X) = [k=>v for (k,v) ∈ X if length(v)>0]
    @inline trim_imag(X)  = (abs(imag(X))<1e-7 ? real(X) : X)
    LATT = HQWannier.LATT
    HQS = kspace_hopping_hamiltonian_symbolic(
                Dict(k=>(x->(:ALL,(x,))) for k in keys(LATT.LN.SPNB)),  #TODO
                Dict(k=>[k,]             for k in keys(LATT.LN.SPNB)),
                Dict(m=>(x->x) for m ∈ LATT.UC.m),
                Dict(m=>[m,]   for m ∈ LATT.UC.m),
                LATT;
                da = 1.0
    ) ;
    #TODO
    ret = Dict(k=>[] for (k,dic) ∈ HQS.MATS)
    ret_k = []
    MAT = HQWannier.MAT
    for (k,dic) ∈ HQS.MATS
        ret_k = []
        for (q,m) ∈ dic
            X,Y,V = findnz(m)
            push!(ret_k, q=>[ (x,y)=>trim_imag(MAT[q][x,y])
                                   for (x,y) ∈ zip(X,Y)
                                       if abs(MAT[q][x,y])>thr ])
        end
        ret[k] = trim_empty(ret_k)
    end
    return Dict(trim_empty(ret))
end

## ===========================================================

@inline show_hop(v) = println(
                        join([ "   "*string(v[1])*" => ", 
                               sort(["\t"*string(k=>p) for (k,p) ∈ v[2]])...],
                             "\n")
)
@inline show_t(k,v) = (println(string(k)*" => "); show_hop.(v))

R = SP2T_raw(Bi2Pd_HQ_Wann; thr=0.5) ;

[show_t(k,R[k]) for k ∈ sort(collect(keys(R)))] ;