# rotate vectors, to correct a[:,1] (the first real space basis)

# --------------   rotate   --------------

# this function is revised from Mathematica code
# returns 3x3 matrix that rotate [x1,x2,x3] to [1,0,0]
function rotate_3d(basis::Coordinates)::Coordinates
    List(X...) = (eltype(X)<:Vector) ? hcat(X...) : Real[x for x ∈ X]
    Power(x,n) = x^n
    Sqrt(x) = sqrt(x)
    x1,x2,x3 = basis[:,1]
    transz = List( List(x1/Sqrt(Power(x1,2) + Power(x2,2)),-(x2/Sqrt(Power(x1,2) + Power(x2,2))),0),
                    List(x2/Sqrt(Power(x1,2) + Power(x2,2)),x1/Sqrt(Power(x1,2) + Power(x2,2)),0),
                    List(0,0,1) )
    transy = List(  List( Sqrt(Power(x1,2) + Power(x2,2))/Sqrt(Power(x1,2) + Power(x2,2) + Power(x3,2)),0,
                          -(x3/Sqrt(Power(x1,2) + Power(x2,2) + Power(x3,2))) ),
                    List(0,1,0),
                    List( x3/Sqrt(Power(x1,2) + Power(x2,2) + Power(x3,2)),0,
                          Sqrt(Power(x1,2) + Power(x2,2))/Sqrt(Power(x1,2) + Power(x2,2) + Power(x3,2)) )  )

    return transy * transz
end

# returns 2x2 part of the 3x3 matrix that rotate [x1,x2,x3=0] to [1,0,0]
function rotate_2d(basis::Coordinates)::Coordinates
    basis3 = zeros(eltype(basis), 3, 3)
    basis3[1:2,1:2] = basis[:,:]
    basis3[3,3] = one(eltype(basis))
    return rotate_3d(basis3)[1:2,1:2]
end

rotate(basis::Coordinates)::Coordinates = ( (size(basis,1)==2) ? rotate_2d(basis) : rotate_3d(basis) )
# --------------   rotate   --------------
# ----------------------------------------
