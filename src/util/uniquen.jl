
#NOTE see discussion in https://github.com/JuliaLang/julia/issues/19147
#NOTE and https://discourse.julialang.org/t/the-test-isequal-0-0-0-0-returns-false/38404/5
@inline nroundc(myarray,n::Int64) = map(  x->(round( x,digits=n) +(+0.0+0.0im)),  myarray)
@inline nroundr(myarray,n::Int64) = map(  x->(round( x,digits=n) +(+0.0      )),  myarray)
@inline nroundvr(myarray,n::Int64) = map( x->(round.(x,digits=n).+(+0.0      )),  myarray)
@inline nroundvc(myarray,n::Int64) = map( x->(round.(x,digits=n).+(+0.0+0.0im)),  myarray)
@inline uniquen(myarray::Array{T,N},n::Int64) where {T<:Complex,N}   = unique(nroundc(myarray,n))
@inline uniquen(myarray::Array{T,N},n::Int64) where {T<:Real,N}      = unique(nroundr(myarray,n))
@inline uniquen(myarray::Vector{Vector{T}},n::Int64) where {T<:Real}    = unique(nroundvr(myarray,n))
@inline uniquen(myarray::Vector{Vector{T}},n::Int64) where {T<:Complex} = unique(nroundvc(myarray,n))
