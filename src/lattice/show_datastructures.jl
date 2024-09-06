function print(io::IO, UC::UnitCell)
    @inline _tostr_(mat) = join([string(mat[:,i]) for i=1:size(mat,2)],"\n         ")    
    @inline _chop_(l,n_per_l) = [l[(i-1)*n_per_l+1:i*n_per_l] for i=1:(length(l)÷n_per_l)]
    @inline _tostr2_(l) = join([join(string.(l[i]), " , ") for i=1:length(l)],"\n         ")
    str = join([ 
            "dim   :: $(UC.dim)",
            "nsubl :: $(UC.nsubl)",
            "a     :: $(_tostr_(UC.a))",
            "δ     :: $(_tostr_(UC.δ))",
            "m     :: $(_tostr2_(_chop_(UC.m, 8)))",
            "ξ     :: $(_tostr2_(_chop_(UC.ξ, 8)))", ], "\n")
    Base.print(io, str)
    return nothing
end

print(UC::UnitCell) = print(Base.stdout, UC)

#TODO show automatically