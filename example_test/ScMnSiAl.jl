global const WKSPC = "/home/dabajabaza/jianguoyun/Workspace/TiMnScSi"

#  --------------------------------------------

using BSON
using DrWatson
using LinearAlgebra
using SparseArrays
using Pkg
Pkg.activate("/home/dabajabaza/jianguoyun/Workspace/LatticeLab")
using LatticeLab
Pkg.activate("/home/dabajabaza/jianguoyun/Workspace/BandStructures")
using BandStructures
using PyCall
phonopy = pyimport("phonopy")

# --------------------------------------------

function lattice_dynmat_from_phonopy(
    ph::PyCall.PyObject,
    POSCAR::String, 
    FORSE_SETS::String, 
    CONF::String;
    max_distance=30.0,
    maxlen=20,
    margin=6, 
    bounding_box=([-0.001,-0.0004,-0.0006],[1 0 0; 0 1 0; 0 0 1],[1,1,1],[true,true,true]),
    rounding_digits=5,
    PHONOPY_EPS=1e-6,
    PHONOPY_PREC_CTRL=8
    )

    println("\n-------------------------------\nlattice_dynmat_from_phonopy():")
    confs = [   "max_distance=$max_distance",
                "maxlen=$maxlen",
                "margin=$margin",
                "bounding_box=$bounding_box",
                "rounding_digits=$rounding_digits",
                "PHONOPY_EPS=$PHONOPY_EPS",
                "PHONOPY_PREC_CTRL=$PHONOPY_PREC_CTRL",  ]
    @info "Configurations:\n" * join(confs, "\n\t")
    println("")

    @info "LatticeLab.PhonopyCopy_from_vasp"
    @time PH = LatticeLab.PhonopyCopy_from_vasp(ph, POSCAR, FORSE_SETS, CONF)

    @info "LatticeLab.all_links_between_sublattices"
    @time ALL_LN = LatticeLab.all_links_between_sublattices(
                                    PH.UC; 
                                    maxlen=maxlen,
                                    max_distance=max_distance,
                                    large_enough_margin=margin, 
                                    bounding_box=bounding_box,
                                    rounding_digits=rounding_digits )

    @info "LatticeLab.extract_force_constants"
    @time SP, NB = LatticeLab.extract_force_constants(
                                    PH, ALL_LN;
                                    SVEC_PREC_CTRL=PHONOPY_PREC_CTRL, 
                                    FC_PREC_CTRL=PHONOPY_PREC_CTRL, 
                                    SVEC_NORM_EPS=PHONOPY_EPS,
                                    FC_EPS=PHONOPY_EPS )

    @info "LatticeLab.build_lattice"
    @time LATT = LatticeLab.build_lattice(NB, bounding_box)

    @info "LatticeLab.kspace_dynamical_matrix"
    @time DM = LatticeLab.kspace_dynamical_matrix(SP, PH.MASS, LATT)

    return PH, DM
end


# --------------------------------------------

function process_phonopy_database_element(
    ph::PyCall.PyObject, 
    MP::String, 
    FD::String; 
    rounding_digits=4, 
    margin=6,
    max_distance=20.0,
    maxlen=20,
    )

    POSCAR = "$FD/$MP/POSCAR-unitcell"
    FORcE = "$FD/$MP/" * (isfile("$FD/$MP/FORCE_SETS") ? "FORCE_SETS" : "FORCE_CONSTANTS")
    @assert isfile(FORcE)
    CONF = "$FD/$MP/phonopy.conf"
    PH, DM, err = (nothing, nothing, nothing)
    try
        PH, DM = lattice_dynmat_from_phonopy(
                        ph, POSCAR, FORcE, CONF,
                        rounding_digits=rounding_digits, 
                        max_distance=max_distance,
                        maxlen=maxlen,
                        margin=margin
        )
    catch _e_
        err = string(typeof(_e_))
        rethrow(_e_)
    end
    BSON.bson(  "$FD/$MP/LatticeLab_parse_result_$(MP).bson", 
                RESULTS=(PH=PH,
                         rounding_digits=rounding_digits, margin=margin,
                         err=err)   )
    return DM, err
end


# --------------------------------------------


function plot_band_ScMnSiAl(
    proj_func,
    figure_file_name,
    DM, 
    HighSymmRel, 
    Mass; 
    Δk=0.02,
    COLOR=["red","blue","green","orange"], 
    SIZE = 40.0,
    settings=Dict(  :colors=>["black",],
                    :lw=>0.4,
                    :range=>nothing,
                    :fontsize=>12,
                    :markerstrokealpha=>0.7,
                    :markerstrokewidth=>0.2,
                    :markersizetrimratio=>0.01,
                    :aspect_ratio=>3,
                    :figure_size => (10,4),   )
    )
    M1  = Dict(el=>Mass[el] for (el,m) in DM.M)
    DM1 = LatticeLab.kspace_dynamical_matrix(DM.SP, M1, DM.LATT)
    b   = LatticeLab.reciprocal_basis(DM.LATT.UC)
    HighSymmAbs = [[name=>b*qrel for (name,qrel) ∈ HS] for HS ∈ HighSymmRel]
    BS  = [LatticeLab.band_structure_with_eigenvectors(HSA, DM1; Δk=Δk, test_eigen=false, eps=1e-8)  for HSA ∈ HighSymmAbs]
    OPS1= proj_func(DM1.LATT.UC)
    BS1 = [ compute_band_markers(b, OPS1) |> BandStructures.LatticeLab_bands_BandStructure 
            for b in BS];
    plot_bands(figure_file_name, BS1, dpi=1500, COLOR=COLOR, SIZE=SIZE, settings=settings) ;
    return BS1
end

## --------------------------------------------

LTOPS(UC) = [(k,v)->abs(v'*LongitudinalProjector(k,UC)*v), 
             (k,v)->abs(v'*TransverseProjector(k,UC)*v) ]

LOPS_EL(UC,EL) = [ (k,v)->abs(v'*AtomProjector(EL,UC)*LongitudinalProjector(k,UC)*v) ]
LOPS_Sc(UC) = LOPS_EL(UC,:Sc)
LOPS_Mn(UC) = LOPS_EL(UC,:Mn)
LOPS_Al(UC) = LOPS_EL(UC,:Al)
LOPS_Si(UC) = LOPS_EL(UC,:Si)
LOPS_Na(UC) = TOPS_EL(UC,:Na)
LOPS_Sb(UC) = TOPS_EL(UC,:Sb)
LOPS_Cu(UC) = TOPS_EL(UC,:Cu)
LOPS_W(UC)  = LOPS_EL(UC,:W)
LOPS_Se(UC) = LOPS_EL(UC,:Se)


TOPS_EL(UC,EL) = [ (k,v)->abs(v'*AtomProjector(EL,UC)*TransverseProjector(k,UC)*v) ]
TOPS_Sc(UC) = TOPS_EL(UC,:Sc)
TOPS_Mn(UC) = TOPS_EL(UC,:Mn)
TOPS_Al(UC) = TOPS_EL(UC,:Al)
TOPS_Si(UC) = TOPS_EL(UC,:Si)
TOPS_Na(UC) = TOPS_EL(UC,:Na)
TOPS_Sb(UC) = TOPS_EL(UC,:Sb)
TOPS_Cu(UC) = TOPS_EL(UC,:Cu)
TOPS_W(UC)  = TOPS_EL(UC,:W)
TOPS_Se(UC) = TOPS_EL(UC,:Se)


## --------------------------------------------

sk = pyimport("seekpath")
MP = "mp-999998-20180417"
MP = "mp-999901-20180417" # CHEM = "NaSbCu"

MP = "mp-830-20180417" 
MP = "mp-1821-20180417"

# --------------------------------------------

DM, err = process_phonopy_database_element( phonopy, 
                                            MP, 
                                            WKSPC, 
                                            rounding_digits=4,
                                            max_distance=50.0,
                                            maxlen=40,
) ;

uc = BandStructures.UnitCell(3,DM.LATT.UC.nsubl,DM.LATT.UC.a,DM.LATT.UC.δ,DM.LATT.UC.m,DM.LATT.UC.ξ) ;
_kpath = call_get_path(sk, uc; with_time_reversal=false)
HighSymmRel = kpath(_kpath) ;

# --------------------------------------------

function Gamma_kpath(seekpath_result::Dict; ratio=0.1)
    G = [0,0,0]
    PATHALL = []
    for (k,p) ∈ seekpath_result["point_coords"]
        if k!="GAMMA"
            push!(PATHALL, ["Γ"=>G,k=>ratio.*p])
        end
    end
    PATHALL
end

gamma_path = Gamma_kpath(_kpath)

## --------------------------------------------

CHEM = "NaSbCu"
CHEM = "GaN"
CHEM = "WSe"

plot_band_ScMnSiAl(
    LTOPS,
    "$(CHEM)_LT.pdf",
    DM, 
    gamma_path, 
    DM.M; 
    Δk=0.0005,
    COLOR=["red","blue","green","orange"], 
    SIZE = 40.0,
    settings=Dict(  :colors=>["black",],
                    :lw=>0.4,
                    :range=>(0,0.05),
                    :fontsize=>12,
                    :K_sep=>0.01,
                    :markerstrokealpha=>0.7,
                    :markerstrokewidth=>0.2,
                    :markersizetrimratio=>0.01,
                    :aspect_ratio => 8,
                    :figure_size  => (8,4),   )
) ;

## --------------------------------------------

for (atm, proj, col) in [("Sc",LOPS_Sc,"violet"),("Si",LOPS_Si,"grey"),("Al",LOPS_Al,"red"),("Mn",LOPS_Mn,"blue")]
    plot_band_ScMnSiAl(
        proj,
        "ScMnSiAl_L_$(atm).pdf",
        DM, 
        HighSymmRel, 
        DM.M; 
        Δk=0.01,
        COLOR=[col,], 
        SIZE = 50.0,
        settings=Dict(  :colors=>["black",],
                        :lw=>0.4,
                        :range=>(0,0.3),
                        :fontsize=>12,
                        :markerstrokealpha=>0.9,
                        :markerstrokewidth=>0.5,
                        :markersizetrimratio=>0.01,
                        :aspect_ratio=>6,
                        :figure_size => (8,4),   )
    ) ;
end

## --------------------------------------------

for (atm, proj, col) in [("Sc",TOPS_Sc,"violet"),("Si",TOPS_Si,"grey"),("Al",TOPS_Al,"red"),("Mn",TOPS_Mn,"blue")]
    plot_band_ScMnSiAl(
        proj,
        "ScMnSiAl_T_$(atm).pdf",
        DM, 
        HighSymmRel, 
        DM.M; 
        Δk=0.01,
        COLOR=[col,], 
        SIZE = 50.0,
        settings=Dict(  :colors=>["black",],
                        :lw=>0.4,
                        :range=>(0,0.3),
                        :fontsize=>12,
                        :markerstrokealpha=>0.9,
                        :markerstrokewidth=>0.5,
                        :markersizetrimratio=>0.01,
                        :aspect_ratio=>6,
                        :figure_size => (8,4),   )
    ) ;
end

## --------------------------------------------

## --------------------------------------------


A = [("Sc",LOPS_Sc,"violet"),("Si",LOPS_Si,"grey"),("Al",LOPS_Al,"red"),("Mn",LOPS_Mn,"blue")]

A = [("Na",LOPS_Na,"orange"),("Cu",LOPS_Cu,"green"),("Sb",LOPS_Sb,"blue"),]

A = [("W",LOPS_W,"red"),("Se",LOPS_Se,"orange"),]

CHEM = "ScMnSiAl"

CHEM = "NaSbCu"

CHEM = "WSe"


for (atm, proj, col) in A
    plot_band_ScMnSiAl(
        proj,
        "$(CHEM)_L_$(atm)_gamma.pdf",
        DM, 
        gamma_path, 
        DM.M; 
        Δk=0.001,
        COLOR=[col,], 
        SIZE = 50.0,
        settings=Dict(  :colors=>["black",],
                        :lw=>0.4,
                        :range=>(0,0.05),
                        :K_sep=>0.02,
                        :fontsize=>12,
                        :markerstrokealpha=>0.9,
                        :markerstrokewidth=>0.5,
                        :markersizetrimratio=>0.01,
                        :aspect_ratio=>6,
                        :figure_size => (8,4),   )
    ) ;
end


## --------------------------------------------


A = [("Sc",TOPS_Sc,"violet"),("Si",TOPS_Si,"grey"),("Al",TOPS_Al,"red"),("Mn",TOPS_Mn,"blue")]

A = [("Na",TOPS_Na,"orange"),("Cu",TOPS_Cu,"green"),("Sb",TOPS_Sb,"blue"),]

A = [("W",TOPS_W,"red"),("Se",TOPS_Se,"orange"),]

CHEM = "ScMnSiAl"

CHEM = "NaSbCu"

CHEM = "WSe"


for (atm, proj, col) in A
    plot_band_ScMnSiAl(
        proj,
        "$(CHEM)_T_$(atm)_gamma.pdf",
        DM, 
        gamma_path, 
        DM.M; 
        Δk=0.001,
        COLOR=[col,], 
        SIZE = 50.0,
        settings=Dict(  :colors=>["black",],
                        :lw=>0.4,
                        :range=>(0,0.05),
                        :K_sep=>0.02,
                        :fontsize=>12,
                        :markerstrokealpha=>0.9,
                        :markerstrokewidth=>0.5,
                        :markersizetrimratio=>0.01,
                        :aspect_ratio=>6,
                        :figure_size => (8,4),   )
    ) ;
end
