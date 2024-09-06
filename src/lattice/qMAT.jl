@inline qMAT(q,M,BASIS) = sum( cis(dot( q, BASIS * mn )) .* M[mn] for mn ∈ keys(M) )

@inline qMATS(q,BASIS,MS,C) = sum( C[i] .* qMAT(q,MS[i],BASIS) for i=1:length(C) ) #TODO optimize!

@inline ∂qMAT∂q(q,α::Int,M,BASIS) = im.*sum( (BASIS*mn)[α] .* cis(dot(q, BASIS*mn)) .* M[mn] for mn ∈ keys(M) )
