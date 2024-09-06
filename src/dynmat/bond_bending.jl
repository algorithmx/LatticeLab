
#NOTE CONVENTION ∂²Φ / ∂r(i) ∂r(j) ==> H(i,j)
#NOTE H(1,2) C†(1) C(2)
#NOTE Δ === r_end - r_start
#NOTE example
# ((:r12, [0,0,1], [0,1,0]), 1.0, 1e-8)
@inline function check_mode(
    mode::Symbol,
    istart::Int64, iknee::Int64, iend::Int64,
    si::Int64, sj::Int64,
    r10, r20, r21, Δ,
    eps::Float64
    )
    #NOTE the routine check_mode() rules out the incompatible modes
    # ensure that si, sj fits the convention
    # istart(1) --> iknee(0) -> iend(2)
    #NOTE CONVENTION ∂²Φ / ∂r(i) ∂r(j) ==> D(i,j)
    #NOTE Δ === r_end - r_start
    return ( (mode==:r10 && (iknee!=sj || istart!=si || norm(r10.-Δ)>eps))
          || (mode==:r01 && (iknee!=si || istart!=sj || norm(r10.+Δ)>eps))
          || (mode==:r20 && (iknee!=sj || iend  !=si || norm(r20.-Δ)>eps))
          || (mode==:r02 && (iknee!=si || iend  !=sj || norm(r20.+Δ)>eps))
          || (mode==:r12 && (istart!=si || iend !=sj || norm(r21.+Δ)>eps))
          || (mode==:r21 && (istart!=sj || iend !=si || norm(r21.-Δ)>eps)) )
end

# --------------------------------------------------------------

#NOTE
# a bond angle consists three atoms 1-0-2
# to the double counting problem, fix the bond as reference, such that 1-0-2 ≡ 2-0-1
# this function gives ∂² ϕ / ∂r ∂r', where
# (r,r') can be one of the following 6 combinations (called mode in the code):
# (r1,r2)-->:r12, (r2,r1)-->:r21,
# (r1,r0)-->:r10, (r0,r1)-->:r01
# (r2,r0)-->:r20, (r0,r2)-->:r02
# the routine check_mode() rules out the incompatible modes
# the rest part computes the force constant according to the mode

#NOTE equation : kθ*(cosθ-cosθ0) with cosθ0 ≡ -1
function straight_bond_bending(
    Δ::Vector, si::Int64, sj::Int64,
    mode_ske_r10_r20::Tuple,
    kθ::Real,
    eps::Float64
    )::Matrix
    #NOTE istart, iknee, iend are sublattice ID of the 1-0-2 bond bending nodes
    (mode::Symbol, istart::Int64, iknee::Int64, iend::Int64, r10::Vector, r20::Vector) = mode_ske_r10_r20
    r21 = r20 .- r10

    #NOTE enforce the knee position
    if check_mode(mode, istart, iknee, iend, si, sj, r10, r20, r21, Δ, eps)
        return default_force(Δ)
    end

    (r10hat::Vector, r20hat::Vector, N10::Real, N20::Real) = (normalize(r10), normalize(r20), norm(r10), norm(r20))
    γ = N20/N10
    kθhat = kθ/(N10*N20)
    cosθ0 = real(dot(r10hat,r20hat))
    IDEN = eye(length(Δ))
    T11 = r10hat * r10hat'
    T22 = r20hat * r20hat'
    T12 = r10hat * r20hat'
    if mode == :r12 # ∂² ϕ / ∂r1 ∂r2
        return kθhat .* (IDEN .- T11 .- T22 .+ cosθ0.*T12)
    elseif mode == :r21 # ∂² ϕ / ∂r2 ∂r1
        return kθhat .* (IDEN .- T11 .- T22 .+ cosθ0.*T12) |> transpose
    elseif mode == :r10 # ∂² ϕ / ∂r1 ∂r0
        return kθhat .* ((cosθ0*γ-1)  .* IDEN
                      .+ (1-3cosθ0*γ) .* T11
                      .+                 T22
                      .+ (γ-cosθ0)    .* T12
                      .+ (γ)          .* T12' )
    elseif mode == :r01 # ∂² ϕ / ∂r0 ∂r1 # transpose of :r10
        return kθhat .* ((cosθ0*γ-1)  .* IDEN
                      .+ (1-3cosθ0*γ) .* T11
                      .+                 T22
                      .+ (γ-cosθ0)    .* T12
                      .+ (γ)          .* T12' ) |> transpose
    elseif mode == :r20 # exchange 1<-->2 from :r10
        return straight_bond_bending( Δ, si, sj, (:r10,iend,iknee,istart,r20,r10), kθ, eps )
    elseif mode == :r02 # exchange 1<-->2 from :r01
        return straight_bond_bending( Δ, si, sj, (:r01,iend,iknee,istart,r20,r10), kθ, eps )
    else
        throw("ERROR: straight_bond_bending() mode error : "*string(mode))
    end
    return default_force(Δ)
end


function straight_bond_bending(
    Δ::Vector, si::Int64, sj::Int64,
    mode_ske_r10_r20_list::Vector{T},
    kθ::Real,
    eps::Float64
    )::Matrix where {T<:Tuple}
    # mode_ske_r10_r20_list contains a list of different bond specs so that
    # (Δ, si, sj) is compatible to one or more bonds
    sum( straight_bond_bending(Δ,si,sj,mkrr,kθ,eps) for mkrr ∈ mode_ske_r10_r20_list )
end

# --------------------------------------------------------------


#NOTE the double counting problem
# a bond angle consists three atoms 1-0-2 ≡ 2-0-1
# this function gives ∂² ϕ / ∂r ∂r', where
# (r,r') can be one of the following 6 combinations (called mode in the code):
# (r1,r2)-->:r12, (r2,r1)-->:r21,
# (r1,r0)-->:r10, (r0,r1)-->:r01
# (r2,r0)-->:r20, (r0,r2)-->:r02
# the routine check_mode() rules out the incompatible modes
# the rest part computes the force constant according to the mode

#NOTE equation : (kθ/2)*(cosθ-cosθ0)²
function bond_bending(
    Δ::Vector, si::Int64, sj::Int64,
    mode_ske_r10_r20::Tuple,
    kθ::Real,
    eps::Float64
    )::Matrix
    #NOTE istart, iknee, iend are sublattice ID of the 1-0-2 bond bending nodes
    (mode::Symbol, istart::Int64, iknee::Int64, iend::Int64, r10::Vector, r20::Vector) = mode_ske_r10_r20
    r21 = r20 .- r10
    #NOTE enforce the knee position
    if check_mode(mode, istart, iknee, iend, si, sj, r10, r20, r21, Δ, eps)
        return default_force(Δ)
    end
    (r10hat::Vector, r20hat::Vector, N10::Real, N20::Real) = (normalize(r10), normalize(r20), norm(r10), norm(r20))
    kθhat = kθ/(N10*N20)
    cosθ0 = real(dot(r10hat,r20hat))
    dcosθdr1 = r20hat .- cosθ0.*r10hat
    dcosθdr2 = r10hat .- cosθ0.*r20hat
    if mode == :r12
        return kθhat .* (dcosθdr1 * dcosθdr2') # B/2 cosθ²
    elseif mode == :r21
        return kθhat .* (dcosθdr2 * dcosθdr1') # B/2 cosθ²
    elseif mode == :r10
        dcosθdr0_1 = ((-N20/N10).*dcosθdr1) .+ (-dcosθdr2)
        return kθhat .* (dcosθdr1 * dcosθdr0_1') # B/2 cosθ²
    elseif mode == :r01 # transpose of :r10
        dcosθdr0_1 = ((-N20/N10).*dcosθdr1) .+ (-dcosθdr2)
        return kθhat .* (dcosθdr0_1 * dcosθdr1') # B/2 cosθ²
    elseif mode == :r20 # exchange 1<-->2 from :r10
        dcosθdr0_2 = (-dcosθdr1) .+ ((-N10/N20).*dcosθdr2)
        return kθhat .* (dcosθdr2 * dcosθdr0_2') # B/2 cosθ²
    elseif mode == :r02 # exchange 1<-->2 from :r01
        dcosθdr0_2 = (-dcosθdr1) .+ ((-N10/N20).*dcosθdr2)
        return kθhat .* (dcosθdr0_2 * dcosθdr2') # B/2 cosθ²
    else
        throw("ERROR: bond_bending() mode error : "*string(mode))
    end
    return default_force(Δ)
end


function bond_bending(
    Δ::Vector, si::Int64, sj::Int64,
    mode_ske_r10_r20_list::Vector{T},
    kθ::Real,
    eps::Float64
    )::Matrix where {T<:Tuple}
    sum( bond_bending(Δ,si,sj,mkrr,kθ,eps) for mkrr ∈ mode_ske_r10_r20_list )
end


# ----------------------------- SUMMARY ---------------------------------
