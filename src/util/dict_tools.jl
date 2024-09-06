@inline sck(x) = sort(collect(keys(x)))


@inline get_val(k,dic,default) = ((k ∈ keys(dic)) ? dic[k] : default)


@inline zict(x,y) = Dict(zip(x,y))


function list_of_pairs_to_dict_merge_v(lp::Vector{P}) where {P<:Pair} 
    ret = Dict(k=>[] for k ∈ unique(first.(lp)))
    for (k,v) ∈ lp
        ret[k] = [ret[k]..., v...]
    end
    return ret
end


function extend_dict(
    dic_to_ext,
    key_ref;
    default_value = 0
    )
    dic_ret = copy(dic_to_ext)
    dic_to_ext_key = keys(dic_to_ext)
    for key ∈ key_ref
        if key ∉ dic_to_ext_key
            dic_ret[key] = default_value
        end
    end
    return dic_ret
end


function delete_empty_sparse(
    dic::Dict{K,V}
    ) where { K<:Any, V<:SparseMatrixCSC }
    return Dict( k=>dic[k] 
                 for k ∈ keys(dic)
                    if length(dic[k].nzval)>0 && sum(abs.(dic[k]))>1e-20 )
end


function delete_empty(
    dic::Dict{K,V}
    ) where { K<:Any, V<:SparseMatrixCSC }
    return Dict( k=>dic[k] for k ∈ keys(dic)
                            if length(dic[k].nzval)>0 && sum(abs.(dic[k]))>1e-20 )
end


function delete_empty(
    dic::Dict{K,V}
    ) where { K<:Any, V<:Vector }
    return Dict( k=>unique(v) for (k,v) ∈ dic if length(v)>0 )
end
