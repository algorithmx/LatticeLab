using BSON

using Pkg
Pkg.activate("/home/dabajabaza/jianguoyun/Workspace/LatticeLab")
using LatticeLab

using LinearAlgebra
using SparseArrays

using PyCall
phonopy = pyimport("phonopy")

# --------------------------------------------

function lattice_dynmat_from_phonopy(
    ph,
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

    return LATT, PH, SP, NB, DM
end


function process_phonopy_database_element(ph, MP, FD; rounding_digits=4, margin=6)
    POSCAR = "$FD/$MP/POSCAR-unitcell"
    FORSE_SETS = "$FD/$MP/FORCE_SETS"
    CONF = "$FD/$MP/phonopy.conf"
    LATT, PH, SP, NB, DM, err = (nothing, nothing, nothing, nothing, nothing, nothing)
    try
        LATT, PH, SP, NB, DM = lattice_dynmat_from_phonopy(ph, POSCAR, FORSE_SETS, CONF,
                                    rounding_digits=rounding_digits, 
                                    margin=margin) ;
    catch _e_
        err = string(typeof(_e_))
    end
    BSON.bson(  "$FD/$MP/LatticeLab_parse_result_$(MP).bson", 
                RESULTS=(PH=PH, DM=DM,
                         rounding_digits=rounding_digits, margin=margin,
                         err=err)   )
    return err
end


## --------------------------------------------

#using Profile
#Profile.clear()
#Profile.init(n=10^8,delay=0.0002)
global const WKSPC = "/home/dabajabaza/Downloads"

## --------------------------------------------

MP = "mp-6661-20180417"

## --------------------------------------------

MP = "mp-6678-20180417"
#process_phonopy_database_element(phonopy, MP, WKSPC) ;
result = BSON.load("$WKSPC/$MP/LatticeLab_parse_result_$(MP).bson")

## --------------------------------------------

MP = "mp-989189-20180417"
POSCAR = "/home/dabajabaza/Downloads/$MP/POSCAR-unitcell"
FORSE_SETS = "/home/dabajabaza/Downloads/$MP/FORCE_SETS"
CONF = "/home/dabajabaza/Downloads/$MP/phonopy.conf"
LATT, PH, SP, NB, DM = lattice_dynmat_from_phonopy(phonopy, POSCAR, FORSE_SETS, CONF) ;

## --------------------------------------------

POSCAR = "/home/dabajabaza/Downloads/mp-484-20180417/POSCAR-unitcell"
FORSE_SETS = "/home/dabajabaza/Downloads/mp-484-20180417/FORCE_SETS"
CONF = "/home/dabajabaza/Downloads/mp-484-20180417/phonopy.conf"
LATT, PH, SP, NB, DM = lattice_dynmat_from_phonopy(phonopy, POSCAR, FORSE_SETS, CONF) ;

##

POSCAR = "/home/dabajabaza/Downloads/mp-989629-20180417/POSCAR-unitcell"
FORSE_SETS = "/home/dabajabaza/Downloads/mp-989629-20180417/FORCE_SETS"
CONF = "/home/dabajabaza/Downloads/mp-989629-20180417/phonopy.conf"
LATT, PH, SP, NB, DM = lattice_dynmat_from_phonopy(phonopy, POSCAR, FORSE_SETS, CONF, margin=4, max_distance=35.0) ;

##

## ---------------------------------------------------------------

## ---------------------------------------------------------------

## ---------------------------------------------------------------


np = pyimport("numpy")

symprec = 1e-5
reduced_cell_method = "niggli"
reduced_bases = phonopy.structure.cells.get_reduced_bases(PH.SUC.a, method=reduced_cell_method, tolerance=symprec)
trans_mat_float = np.dot(PH.SUC.a, np.linalg.inv(reduced_bases))
trans_mat = Int.(np.rint(trans_mat_float))
@assert all(abs.(trans_mat_float - trans_mat) .< 1e-8)
trans_mat_inv_float = np.linalg.inv(trans_mat)
trans_mat_inv = np.rint(trans_mat_inv_float) .|> Int
@assert all(abs.(trans_mat_inv_float - trans_mat_inv) .< 1e-8)

# Reduce all positions into the cell formed by the reduced bases.
supercell_fracs = np.dot(PH.SUC.δ', trans_mat)
dR_super = np.rint(supercell_fracs)
@show dR_super |> norm
supercell_fracs -= dR_super
#supercell_fracs = np.array(supercell_fracs, dtype="double", order="C")


primitive_fracs = np.dot(PH.UC.δ', trans_mat)
dR_prim = np.rint(primitive_fracs)
@show dR_prim
primitive_fracs -= dR_prim
#primitive_fracs = np.array(primitive_fracs, dtype="double", order="C")


##

uc = ph.get_unitcell()

suc = ph.get_supercell()

pos_suc = suc.get_positions()

pos_uc = uc.get_positions()

basis_suc = suc.get_cell()
inv(basis_uc) * pos_suc'

