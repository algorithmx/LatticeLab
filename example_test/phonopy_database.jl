global const WKSPC = "/home/dabajabaza/Downloads"

#  --------------------------------------------

using BSON
using DrWatson
using LinearAlgebra
using SparseArrays
using Pkg
Pkg.activate("/home/dabajabaza/jianguoyun/Nutstore/LatticeLab")
using LatticeLab
Pkg.activate("/home/dabajabaza/jianguoyun/Nutstore/BandStructures")
using BandStructures
using PyCall
phonopy = pyimport("phonopy")

## --------------------------------------------

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
    margin=6
    )

    POSCAR = "$FD/$MP/POSCAR-unitcell"
    FORSE = "$FD/$MP/" * (isfile("$FD/$MP/FORCE_SETS") ? "FORCE_SETS" : "FORCE_CONSTANTS")
    @assert isfile(FORSE)
    CONF = "$FD/$MP/phonopy.conf"
    PH, DM, err = (nothing, nothing, nothing)
    try
        PH, DM = lattice_dynmat_from_phonopy(
                        ph, POSCAR, FORSE, CONF,
                        rounding_digits=rounding_digits, 
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


function plot_band_variable_mass_(
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
    BSa = [BandStructures.LatticeLab_bands_BandStructure(b.Kpath, b.Bands, b.Markers)  for b in BS]
    UC  =  BandStructures.LatticeLab_UnitCell(DM1.LATT.UC)
    OPS1= proj_func(UC)
    BS1 = [compute_band_markers(b, OPS1) for b in BSa];
    plot_bands(figure_file_name, BS1, dpi=1500, COLOR=COLOR, SIZE=SIZE, settings=settings) ;
    return BS1
end


function plot_band_variable_mass(
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
                    :figure_size => (10,4), )
    )
    plot_band_variable_mass_(
        CommonProjectors_ATOM,
        figure_file_name,
        DM, 
        HighSymmRel, 
        Mass; 
        Δk = Δk,
        COLOR = COLOR,
        SIZE = SIZE,
        settings = settings
    )
end

function plot_band_variable_mass_chirality(
    figure_file_name,
    DM, 
    KREF,
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
                    :figure_size => (10,4), )
    )
    plot_band_variable_mass_(
        UC->CommonProjectors_CH(UC,KREF),
        figure_file_name,
        DM, 
        HighSymmRel, 
        Mass; 
        Δk = Δk,
        COLOR = COLOR,
        SIZE = SIZE,
        settings = settings
    )
end


function plot_band_variable_mass_chirality_atom(
    figure_file_name,
    DM, 
    KREF,
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
                    :figure_size => (10,4),  )
    )
    plot_band_variable_mass_(
        UC->CommonProjectors_CH_ATOM(UC,KREF),
        figure_file_name,
        DM, 
        HighSymmRel, 
        Mass; 
        Δk = Δk,
        COLOR = COLOR,
        SIZE = SIZE,
        settings = settings
    )
end


## --------------------------------------------

settings=Dict(  :colors=>["black",],
                :lw=>0.4,
                :range=>nothing,
                :fontsize=>12,
                :markerstrokealpha=>0.7,
                :markerstrokewidth=>0.2,
                :markersizetrimratio=>0.01,
                :aspect_ratio=>3,
                :figure_size => (10,4),
)

## --------------------------------------------

sk = pyimport("seekpath")

MP = "mp-7084-20180417"
MP = "mp-4495-20180417"
MP = "mp-149-20180417"
MP = "mp-1265-20180417"
MP = "mp-22862-20180417"
MP = "mp-11714-20180417"
MP = "mp-830-20180417"
MP = "mp-380-20180417"
MP = "mp-1821-20180417"
MP = "mp-5014-20180417"
MP = "mp-999999-20180417"
MP = "mp-1572-20180417"
MP = "mp-7140-20180417"

DM, err = process_phonopy_database_element(phonopy, MP, WKSPC, rounding_digits=4) ;
uc = BandStructures.UnitCell(3,DM.LATT.UC.nsubl,DM.LATT.UC.a,DM.LATT.UC.δ,DM.LATT.UC.m,DM.LATT.UC.ξ) ;
_kpath = call_get_path(sk, uc; with_time_reversal=false)
HighSymmRel = kpath(_kpath)

BSz = plot_band_variable_mass_chirality_atom(  
       "$(MP)_$(DrWatson.savename(DM.M)).CH.001.pdf", 
        DM, 
        [0,0,1],
        HighSymmRel, 
        DM.M; 
        Δk=0.02, 
        SIZE=40.0, 
        COLOR=["red","blue","green","orange","violet","brown"], 
        settings=settings   ) ;

BSy = plot_band_variable_mass_chirality_atom(  
       "$(MP)_$(DrWatson.savename(DM.M)).CH.010.pdf", 
        DM, 
        [0,1,0],
        HighSymmRel, 
        DM.M; 
        Δk=0.02, 
        SIZE=40.0, 
        COLOR=["red","blue","green","orange","violet","brown"], 
        settings=settings   ) ;

BSx = plot_band_variable_mass_chirality_atom(  
       "$(MP)_$(DrWatson.savename(DM.M)).CH.100.pdf", 
        DM, 
        [1,0,0],
        HighSymmRel, 
        DM.M; 
        Δk=0.02, 
        SIZE=40.0, 
        COLOR=["red","blue","green","orange","violet","brown"], 
        settings=settings   ) ;

## --------------------------------------------



## --------------------------------------------

"http://lamp.tu-graz.ac.at/~hadley/ss1/bzones/fcc.php"
br1 = [1,-1,1]
br2 = [1,1,-1]
br3 = [-1,1,1]

pL = (1//2).*(br1+br2+br3)
pW = ((1//4)*br1+(3//4)*br2+(1//2)*br3)
pK = ((3//8)*br1+(3//4)*br2+(3//8)*br3)
pX = ((0//1)*br1+(1//2)*br2+(1//2)*br3)
pU = ((1//4)*br1+(5//8)*br2+(5//8)*br3)

HighSymmRel = [
"G" => [0,0,0],
"L" => [1,1,1]//2,
"W" => [2,1,3]//4,
"X" => [1,0,1]//2,
"G" => [0,0,0]
]

## --------------------------------------------

MP = "mp-989189-20180417"
MP = "mp-6661-20180417"
MP = "mp-6661-20180417"
MP = "mp-6948-20180417"
MP = "mp-588-20180417"


MP = "mp-1265-20180417"
MP = "mp-190-20180417"
MP = "mp-672-20180417"







bravais(_kpath)

reduced_unitcell_from_seekpath_result(simple_cubic_kpath, uc0)


## --------------------------------------------

## --------------------------------------------


## --------------------------------------------

settings=Dict(  :colors=>["black",],
                :lw=>2.2,
                :range=>nothing,
                :fontsize=>64,
                :markerstrokealpha=>0.8,
                :markerstrokewidth=>2.4,
                :markersizetrimratio=>0.10,
                #:cycle=>[false,]
            )

# --------------------------------------------

HighSymmRel = [
    "G" =>[0,0,0],
    "M" =>[1,0,0]//2, 
    "K" =>[1,1,0]//3, 
    "G" =>[0,0,0], 
    "A" =>[0,0,1]//2,
    "L" =>[1,0,1]//2, 
    "H" =>[1,1,3//2]//3, 
    "A" =>[0,0,1]//2,
]

MP = "mp-4495-20180417"
MP = "mp-569346-20180417"
MP = "mp-830-20180417"
MP = "mp-380-20180417"
MP = "mp-1821-20180417"
MP = "mp-1572-20180417"

DM, err = process_phonopy_database_element(phonopy, MP, WKSPC,rounding_digits=4) ;

##

M1 = DM.M  # Dict(:Cu => 10, :I => 80)

BS = plot_band_variable_mass(  "$(MP)_$(DrWatson.savename(M1)).pdf", 
                                DM, 
                                HighSymmRel, 
                                M1; 
                                Δk=0.03, 
                                SIZE=40.0, 
                                COLOR=["red","blue","green","orange"], 
                                settings=settings   ) ;

## --------------------------------------------


M1 = Dict(:Cu => 80, :I => 10)

BS = plot_band_variable_mass(  "$(MP)_$(DrWatson.savename(M1)).pdf", 
                                DM, 
                                HighSymmRel, 
                                M1; 
                                Δk=0.04, 
                                SIZE=40.0, 
                                COLOR=["red","blue","green","orange"], 
                                settings=settings   ) ;
## --------------------------------------------
 
M1 = Dict(:Cu => 80, :I => 80)

BS = plot_band_variable_mass(  "$(MP)_$(DrWatson.savename(M1)).pdf", 
                                DM, 
                                HighSymmRel, 
                                M1; 
                                Δk=0.04, 
                                SIZE=40.0, 
                                COLOR=["red","blue","green","orange"], 
                                settings=settings   ) ;
## --------------------------------------------

M1 = Dict( :I  => 130, :Cu => 60)

BS = plot_band_variable_mass(  "$(MP)_$(DrWatson.savename(M1)).pdf", 
                                DM, 
                                HighSymmRel, 
                                M1; 
                                Δk=0.04, 
                                SIZE=40.0, 
                                COLOR=["red","blue","green","orange"], 
                                settings=settings   ) ;

## --------------------------------------------

M1 = DM.M

BS = plot_band_variable_mass(  "$(MP)_$(DrWatson.savename(M1)).pdf", 
                                DM, 
                                HighSymmRel, 
                                M1; 
                                Δk=0.04, 
                                SIZE=40.0, 
                                COLOR=["red","blue","green","orange"], 
                                settings=settings   ) ;

## --------------------------------------------


#sk = pyimport("seekpath")

MP = "mp-5014-20180417"

DM, err = process_phonopy_database_element(phonopy, MP, WKSPC) ;

uc = BandStructures.UnitCell(3,DM.LATT.UC.nsubl,DM.LATT.UC.a,DM.LATT.UC.δ,DM.LATT.UC.m,DM.LATT.UC.ξ) ;

_kpath = call_get_path(sk, uc; with_time_reversal=false)

HighSymmRel = kpath(_kpath)

DM, err = process_phonopy_database_element(phonopy, MP, WKSPC) ;

BS = plot_band_variable_mass_chirality(
    "$(MP)_$(DrWatson.savename(DM.M)).CH.1.pdf", 
    DM, 
    HighSymmRel, 
    DM.M; 
    Δk=0.03, 
    SIZE=40.0, 
    COLOR=["red","blue","green","orange"], 
    settings=settings
) ;

## --------------------------------------------

M1 = Dict(:Cu => 100, :Ag => 1, :S => 10)
plot_band_variable_mass("$(MP)_$(DrWatson.savename(M1)).pdf", DM, HighSymmRel, M1; Δk=0.05, SIZE=5.0, COLOR=["red","blue","green"], settings=settings) ;

M1 = Dict(:Cu => 10, :Ag => 100, :S => 1)
plot_band_variable_mass("$(MP)_$(DrWatson.savename(M1)).pdf", DM, HighSymmRel, M1; Δk=0.05, SIZE=5.0, COLOR=["red","blue","green"], settings=settings) ;

M1 = Dict(:Cu => 1, :Ag => 10, :S => 100)
plot_band_variable_mass("$(MP)_$(DrWatson.savename(M1)).pdf", DM, HighSymmRel, M1; Δk=0.05, SIZE=5.0, COLOR=["red","blue","green"], settings=settings) ;

## --------------------------------------------


## --------------------------------------------

HighSymmRel = [
    "G"=>[0,0,0], 
    "M"=>[1,0,0]//2,
    "X"=>[1,1,0]//2, 
    "Y"=>[1,1,1]//2, 
    "G"=>[0,0,0],
]

MP = "mp-7084-20180417"
DM, err = process_phonopy_database_element(phonopy, MP, WKSPC) ;

plot_band_variable_mass_chirality(
    "$(MP)_$(DrWatson.savename(M1)).CH.pdf", 
    DM, 
    HighSymmRel, 
    DM.M; 
    Δk=0.02, 
    COLOR=["red","blue","green"], 
    settings=settings
)

## --------------------------------------------

M1 = Dict(:Sr => 1, :Ca => 10, :Si => 100)
plot_band_variable_mass("$(MP)_$(DrWatson.savename(M1)).pdf", DM, HighSymmRel, M1; Δk=0.05, COLOR=["red","blue","green"], settings=settings)

M1 = Dict(:Sr => 10, :Ca => 1, :Si => 100)
plot_band_variable_mass("$(MP)_$(DrWatson.savename(M1)).pdf", DM, HighSymmRel, M1; Δk=0.05, COLOR=["red","blue","green"], settings=settings)

M1 = Dict(:Sr => 100, :Ca => 1, :Si => 10)
plot_band_variable_mass("$(MP)_$(DrWatson.savename(M1)).pdf", DM, HighSymmRel, M1; Δk=0.05, COLOR=["red","blue","green"], settings=settings)

## --------------------------------------------

M1 = Dict(:K => 1, :Au => 10, :S => 100)
