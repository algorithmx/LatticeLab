
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

