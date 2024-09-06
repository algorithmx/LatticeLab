mutable struct ForceConstant
    UC::LatticeInfo
    FCLN::Dict          #+ name SPLN or ... ?
end


@inline is_ForceConstant(fc) = ( Set(fieldnames(ForceConstant))==Set(fieldnames(typeof(fc)))
                                && is_UnitCell(fc.UC) )


copy(fc::ForceConstant) = ForceConstant( copy(fc.UC), copy(fc.FCLN) )


check_compat(fc::ForceConstant) = true #XXX


@inline convert(::Type{T}, fc::T) where {T<:ForceConstant} = fc


function convert(::Type{T}, fc) where {T<:ForceConstant}
    #XXX other kinds of LatticeInfo ?
    TUC = is_UnitCell(fc.UC) ? UnitCell : LatticeInfo
    ForceConstant( convert(TUC, fc.UC), Dict(fc.FCLN)  )
end
