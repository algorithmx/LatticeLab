@inline ⨰(l::Vector,n::Int) = repeat(l,n)

@inline ⨰(x::Number,n::Int) = fill(x,n)

@inline ⨰(s::Symbol,n::Int) = fill(s,n)

@inline ⨰(s::String,n::Int) = repeat([s,],n)

@inline ⋅(w,v) = dot(w,v)


@inline f1e(minus_n)  = parse(Float64,  "1e-"*string(minus_n))

@inline bf1e(minus_n) = parse(BigFloat, "1e-"*string(minus_n))


@inline cross3(u::Vector,v::Vector) = [
    u[2]*v[3]-u[3]*v[2],
    u[3]*v[1]-u[1]*v[3],
    u[1]*v[2]-u[2]*v[1],
]

@inline ispermute(P::Vector{Int64}) = sort(unique(P))==collect(1:length(P))

@inline findnz( dense_mat::Matrix ) = findnz(sparse(dense_mat))

@inline findnz( dense_vec::Vector ) = findnz(sparse(dense_vec))

@inline findz(  dense_vec::Vector ) = sort(findall(dense_vec.==0))

@inline speye(m) = SparseArrays.spdiagm(0=>repeat([1.0,],m))

@inline   eye(m) = LinearAlgebra.Matrix(1.0I,m,m)

@inline Eij(i,j,N) = sparse([i,],[j,],[1,], N, N)

function BIJ(I,J,B,N)
    #% N is the dimension of the sparse matrix
    @assert length(I)==size(B,1) && length(J)==size(B,2)
    IJ = CartesianIndices((I,J))[:]
    return sparse( map(x->x[1],IJ), map(x->x[2],IJ), B[:], N, N )
end

#=
A = [[1,0] [0,2]]
A1 = kron(Eij(1,2,2),A)
A2 = BIJ(1:2,3:4,A,4)
=#


function leftpad(m::Int64,n::Int64,B::SparseMatrixCSC)
    Ib, Jb, Vb = findnz(B)
    return sparse( vcat(m.+Ib), vcat(n.+Jb), vcat(Vb),
                   m+size(B,1), n+size(B,2) )
end

function rightpad(A::SparseMatrixCSC,m::Int64,n::Int64)
    Ia, Ja, Va = findnz(A)
    return sparse( vcat(Ia), vcat(Ja), vcat(Va),
                   m+size(A,1), n+size(A,2) )
end

function ⊕(A::SparseMatrixCSC, B::SparseMatrixCSC)
    sizeA = size(A)
    Ia, Ja, Va = findnz(A)
    sizeB = size(B)
    Ib, Jb, Vb = findnz(B)
    return sparse( vcat(Ia,sizeA[1].+Ib),
                   vcat(Ja,sizeA[2].+Jb),
                   vcat(Va,Vb),
                   sizeA[1]+sizeB[1],
                   sizeA[2]+sizeB[2] )
end


# delete `nothing` from a list
#@inline delete_nothing(LIST::Vector) = [ v for v ∈ LIST if v!=nothing ]

@inline delete_non_positive(LIST::Vector{Int64}) = [ v for v ∈ LIST if v>0 ]

@inline symbol_join(symb1,symb2) = Symbol( string(symb1) *"_"* string(symb2) )

@inline symbol_suffix(sym,suf::Int) = Symbol(String(sym)*((suf==1) ? "" : string(suf-1)))


@inline function test_registered_pairs(registered_pairs)
    (  all([(jj,ii) ∈ registered_pairs for (ii,jj) ∈ registered_pairs if ii<=jj]) &&
       all([(jj,ii) ∈ registered_pairs for (ii,jj) ∈ registered_pairs if ii>=jj])  )
end


function normalize_dim1(M,EPS)
    MM = M'*M
    gram = diag(MM)
    println("[DEBUG --- normalize_dim1()] norm = ", norm(MM .- diagm(0=>gram)))
    @assert norm(MM .- diagm(0=>gram))<EPS
    return M * spdiagm(0=>sqrt.(1.0./gram))
end


function normalize_dim2(M,EPS)
    MM = M*(M')
    gram = diag(MM)
    @assert norm(MM .- diagm(0=>gram))<EPS
    return spdiagm(0=>sqrt.(1.0./gram)) * M
end


## ========================================================
# hopping

@inline dispatch_params(params,Fdic,Pdic,defval) =
    Dict( kf => Fdic[kf](Tuple(get(params,kp,defval) for kp ∈ Pdic[kf])...)
          for kf ∈ sck(Fdic) )


@inline ⇒(X,Y) = dispatch_params(X,Y...)

