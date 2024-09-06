ENV["PSEUDO_ROOT"] = "/home/dabajabaza/abinitio/pseudopotentials"
include("../../Biphenylene/include/using.jl")
include("../../Biphenylene/include/tools.jl")
include("../../Biphenylene/include/CIF2LN.jl")
include("../../Biphenylene/include/wannier90_utils.jl")

global const WKSPC = "/mnt/dg_hpc/MXene/Mo2N/ps___SG15___kp_16,16,1,0,0,0_kp_16,16,1,0,0,0_cut_100.0,400.0"

projfn           = "$WKSPC/Mo2N_projwfc/Mo2N.proj.projwfc_up"
proj_out_fn      = "$WKSPC/Mo2N_projwfc/Mo2N.projwfc.x.out"
bandfn           = "$WKSPC/Mo2N_nscf/Mo2N.pw.x.out"
band_highsymm_fn = "$WKSPC/Mo2N_band/Mo2N.pw.x.out"
scffn            = "$WKSPC/Mo2N_scf/Mo2N.pw.x.out"
band_w90_fn      = "$WKSPC/Mo2N_wannier90/Mo2N_band.dat"
hr_w90_fn        = "$WKSPC/Mo2N_wannier90/Mo2N_hr.dat"
MoN2_cif_fn      = "$WKSPC/../Mo2N.0.cif"


function pw_bands_from_xml_unit_eV(band_xml_file)
    __hatree_unit_in_eV__ = 27.211324570273 # Electron Volt
    @inline parsenf(s) = parse.(Float64,split(s,r"\s",keepempty=false))
    xml = XMLDict.parse_xml(join(readlines(band_xml_file),"\n")) ;
    ks_energies_dics = xml["output"]["band_structure"]["ks_energies"] ;
    k_points = hcat([parsenf(ks_en["k_point"][""]) 
                     for ks_en in ks_energies_dics]...)
    ks_energies = hcat([__hatree_unit_in_eV__ .* parsenf(ks_en["eigenvalues"][""]) 
                        for ks_en in ks_energies_dics]...)
    return k_points, ks_energies
end

function Mo2N(fn)
    ORIG = [-8e-5, -1e-4, -7e-5]
    ID3  = [1 0 0; 0 1 0; 0 0 1]
    @inline select_non_empty(D) = Dict(k=>v for (k,v) ∈ D if length(v)>0)
    ξξ = [[:s,:d1,:d2,:d3,:d4,:d5],
          [:s,:d1,:d2,:d3,:d4,:d5],
          [:px,:py,:pz]] .|> LatticeLab.Orbits
    UC = CIF2UC(CIF(fn), ξξ)
    ALL_LN     = all_links_between_sublattices( 
                   UC; 
                   maxlen=12, max_distance=22, 
                   rounding_digits=8 ) |> select_non_empty
    ALL_LN_nnn = pick_nearest_by_ij(
                   ALL_LN, 
                   merge( Dict((i,i)=>0 for i=1:6), 
                          Dict((i,j)=>2 for i=1:6 for j=1:6 if i!=j) ); 
                   digits=6
    )
    SP = all_links_dict_to_sp_no_directions(ALL_LN_nnn) ;
    LN = LinkInfo(copy(UC), SP) ;
    return build_lattice(LN, BoundingBox((ORIG,ID3,[3,3,1],[true,true,false]))) ;
end

function band_structure_from_pw_xml_projwfc_up(
    xml_fn::String, 
    projwfc_up_fn::String,
    alat::Float64, 
    kpath_abs,
    atom_orbit_list
    )
    ## construct BS
    kp, en = pw_bands_from_xml_unit_eV(xml_fn)
    kpabs  = ((2π/alat).*kp)
    BS = BandStructures.QE_bandsx_out(
                Pair{String, Tuple{Vector{T} where T, Int64}}[kpath_abs...,], 
                [kpabs[:,i] => en[:,i] for i=1:size(en,2)]
    )
    ## add markers
    proj, states = projwfcx_output_projwfc_up(readlines(projwfc_up_fn))
    @info join(string.(states), "\n")
    nkp = size(proj,1)
    nbd = size(proj,2)
    nop = length(atom_orbit_list)
    @inline select_atom_l(atm, orb) = 
        [ i for (i,t) in enumerate(states) 
            if occursin(atm,t.atom) && t.orbit==orb ]
    @inline weight_atm_l(atm, orbit) = 
        reshape(  sum(proj[:,:,select_atom_l(atm, orbit)], dims=3), nkp, nbd  )
    wt = zeros(Float64, nkp, nbd, nop)
    for (i,(atm,orb)) ∈ enumerate(atom_orbit_list)
        wt[:,:,i] = weight_atm_l(atm,orb)
    end
    ik = 0
    BS.Markers = [ 
        kp_lb => [((ik+=1); (k, Matrix{ComplexF64}(wt[ik,:,:])) ) 
                   for (k,e) ∈ kp_en_pairs]
        for (kp_lb,kp_en_pairs) ∈ BS.Bands
    ]
    ## return 
    return BS
end

const KPATH_REL = [ 
    "Γ"  =>  [0,    0,   0 ],
    "K"  =>  [1//3, 1//3,0 ],
    "M"  =>  [1//2, 0,   0 ],
    "Γ"  =>  [0,    0,   0 ],  ]
const KPATH_REL_I = [ 
    "Γ"  =>  ([0,    0,   0 ], 100),
    "K"  =>  ([1//3, 1//3,0 ], 50 ),
    "M"  =>  ([1//2, 0,   0 ], Int(floor(50*sqrt(3)))),
    "Γ"  =>  ([0,    0,   0 ], 1), ]
@inline kvec_rel2abs(rel,  b) = [k=>(b*vrel) for (k,vrel) in rel]
@inline kveci_rel2abs(rel, b) = [k=>(b*vrel[1],vrel[2]) for (k,vrel) in rel]
@inline kvec_rel(rel, b) = [k=>(b*vrel) for (k,vrel) in rel]
@inline HSP(b)  = kvec_rel2abs(KPATH_REL,b)
@inline HSPi(b) = kveci_rel2abs(KPATH_REL_I,b)


Mo2N_LATT = Mo2N(MoN2_cif_fn) ;

HQw = make_HWannQ(
        readlines(hr_w90_fn), 
        Mo2N_LATT; 
        thr=0.001
) ;

BS_Wann = LatticeLab.band_structure(
    HSP(reciprocal_basis(Mo2N_LATT)), HQw, Δk=1e-3
) |> LatticeLab_bands_BandStructure;

BS_DFT  = band_structure_from_pw_xml_projwfc_up(
    "$WKSPC/Mo2N_scf/Mo2N.band.xml", 
    "$WKSPC/Mo2N_band_projwfc/Mo2N.proj.projwfc_up", 
    Mo2N_LATT.UC.a[1,1], 
    HSPi(reciprocal_basis(Mo2N_LATT)),
    [("Mo","4S"),("Mo","4P"),("Mo","4D"),("Mo","5S"),("N","2P")]
) ;

Mo2N_plot_settings = Dict(:line_colors=>["grey", "black"],
                          :lw=>[0.4,0.7],
                          :ref_levels=>[(-3.8617,"gray"),(0.0,"gray")],
                          :K_sep=>-1,
                          :range=>(-21.0,15.0),   
                          :markerstrokealpha=>0.8,
                          :markerstrokewidth=>0.4,
                          :aspect_ratio=>0.4,
                          :figure_size=>(8,16),)

plot_bands("MoN2_DFT_black_Wann_red.pdf", 
           [BS_DFT, BS_Wann],
           SIZE = [0, 0, 0, 0, 0],
           dpi  = 600,
           settings= merge(Mo2N_plot_settings,
                        Dict(:line_colors=>["black", "red"],
                             :lw=>[1.3,0.7],)))

plot_bands("MoN2_DFT_Wann_Mo_4S+4P.pdf", 
           [BS_DFT, BS_Wann],
           SIZE = [120.0, 120.0, 0, 0, 0],
           dpi  = 600,
           settings= Mo2N_plot_settings)

plot_bands("MoN2_DFT_Wann_Mo_4D+5S.pdf", 
           [BS_DFT, BS_Wann],
           SIZE = [0, 0, 120, 120, 0],
           dpi  = 600,
           settings = Mo2N_plot_settings)

plot_bands("MoN2_DFT_Wann_N_2P.pdf", 
           [BS_DFT, BS_Wann],
           SIZE = [0, 0, 0, 0, 120],
           dpi  = 600,
           settings = Mo2N_plot_settings )