function check_pbc(latt::Lattice)::Bool
    a = latt.UC.a
    Nsites = num_sites(latt)
    origin = latt.BBOX[1]
    invS = inv(a*(latt.BBOX[2]*diagm(0=>latt.BBOX[3])))
    for iv ∈ 1:Nsites
        v = latt.R0[:,iv]
        ieqv = latt.EqV[iv]
        is_in_b = inbbox(v.-origin,invS)
        if (ieqv>0)
            if is_in_b || (! inbbox(latt.R0[:,ieqv].-origin,invS))
                return false
            end
        elseif (ieqv==0)
            if !(is_in_b) return false end
        elseif (ieqv==-1)
            if is_in_b return false end
        end
    end
    return true
end
