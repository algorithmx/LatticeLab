struct SmallVecMismatchError <: Exception
    msg::AbstractString
end

struct DiscontinuousDistanceGroupError <: Exception
    msg::AbstractString
end

struct IncorrectUnitCellError <: Exception
    msg::AbstractString
end

struct DistanceListHeadError <: Exception
    msg::AbstractString
end

struct DistanceListTooShortError <: Exception
    msg::AbstractString
end

struct KeysNotUniqueError <: Exception
    msg::AbstractString
end
