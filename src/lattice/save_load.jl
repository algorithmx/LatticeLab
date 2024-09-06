#TODO separate save/load functions 
#TODO remove JLD2 dependency

function save_kw(
    fn_header::AbstractString,
    keyword::String,
    obj
    )
    JLD2.save( File{:JLD2}(fn_header*".$(keyword).jld2"), "$keyword", obj )
    return
end


@inline save_Lattice(fn0::AbstractString, latt::Lattice) = save_kw(fn0, "Lattice", latt)

@inline save_MinimalLattice(fn0::AbstractString, minilatt::MinimalLattice) = save_kw(fn0, "MinimalLattice", minilatt)

function save_MinimalLattice_VL_EL(fn::String, minilatt::MinimalLattice)
    VL = minilatt.V
    EL = minilatt.E
    JLD2.@save  fn  VL  EL
end

@inline save_Lattice_as_Minimal_VL_EL(fn::String, latt::Lattice) = save_MinimalLattice_VL_EL(fn,MinimalLattice(latt))


function load_kw(
    fn0::AbstractString,
    kw
    )
    (keyword, typ) = kw
    fn = fn0*".$(keyword).jld2"
    if ! isfile(fn)
        @error "load_$(keyword)() : file "*fn*" not found !"
        return nothing
    else
        tmp = JLD2.load(File{:JLD2}(fn))
        return convert(typ, tmp["$keyword"])
    end
end

@inline load_Lattice(fn0::AbstractString) = load_kw(fn0, ("Lattice",Lattice))

@inline load_MinimalLattice(fn0::AbstractString) = load_kw(fn0, ("MinimalLattice",Lattice))

