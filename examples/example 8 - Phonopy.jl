## ===========================================================

using LinearAlgebra
using SparseArrays
using Pkg
Pkg.activate("/home/dabajabaza/jianguoyun/Nutstore/LatticeLab")
using LatticeLab
Pkg.activate("/home/dabajabaza/jianguoyun/Nutstore/BandStructures")
using BandStructures
using PyCall
phonopy = pyimport("phonopy")
sk = pyimport("seekpath")

## ===========================================================

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


## ===========================================================


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
    return DM, err
end


## ===========================================================


function plot_band_mp(
    OPS,
    figure_file_name,
    DM, 
    HighSymmRel; 
    Δk=0.01,
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
    b   = LatticeLab.reciprocal_basis(DM.LATT.UC)
    HighSymmAbs = [[name=>b*qrel for (name,qrel) ∈ HS] for HS ∈ HighSymmRel]
    BS  = [ LatticeLab.band_structure_with_eigenvectors(
                HSA, DM; 
                Δk=Δk, test_eigen=false, eps=1e-8)
            for HSA ∈ HighSymmAbs]
    BS1 = [ compute_band_markers(
                b, OPS
            ) |> BandStructures.LatticeLab_bands_BandStructure
            for b in BS ];
    plot_bands(
            figure_file_name, BS1; 
            dpi=800, COLOR=COLOR, SIZE=SIZE, settings=settings )
    return BS1
end


## ===========================================================

@inline l2b(UC) = BandStructures.UnitCell(3,UC.nsubl,UC.a,UC.δ,UC.m,UC.ξ)

OPS_EL(UC,EL) = [ (k,v)->abs(v'*AtomProjector(EL,UC)*v) ]


## ===========================================================

DM, err = process_phonopy_database_element( phonopy, 
                                            "mp-830-20180417", # CHEM = "GaN"
                                            "./", 
                                            rounding_digits=5,
                                            max_distance=50.0,
                                            maxlen=40 )
HighSymmRel = kpath(call_get_path(sk, l2b(DM.LATT.UC); with_time_reversal=false)) ;

## ===========================================================

ALONG  = [("Ga",LOPS_EL(DM.LATT.UC,:Ga),"red"),
          ("N", LOPS_EL(DM.LATT.UC,:N),"grey") ]

ATTRANS = [("Ga",TOPS_EL(DM.LATT.UC,:Ga),"red"),
           ("N", TOPS_EL(DM.LATT.UC,:N),"grey") ]

for (A,lt) in [(ALONG,"L"), (ATTRANS,"T")]
    for (atm, proj, col) in A
        plot_band_mp(
            proj,
            "GaN_$(lt)_$(atm).pdf",
            DM, 
            HighSymmRel; 
            Δk=0.01,
            COLOR=[col,], 
            SIZE = 50.0,
            settings=Dict(  :colors=>["black",],
                            :lw=>0.3,
                            :fontsize=>6,
                            :K_sep => 0.3,
                            :markerstrokealpha=>0.6,
                            :markerstrokewidth=>0.3,
                            :markersizetrimratio=>0.01,
                            :aspect_ratio=>6,
                            :figure_size => (8,4) ) )
    end
end