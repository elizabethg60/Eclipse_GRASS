function calc_proj_dist(p1, p2)
    return acos(dot(p1, p2)/(norm(p1)*norm(p2)))
end

function NL94_limb_darkening(μ::T, wavelength::T) where T
    μ < zero(T) && return 0.0

    index = findmin(x->abs(x-wavelength), lambda_nm)[2]

    return a0[index] + a1[index]*μ + a2[index]*μ^2 + a3[index]*μ^3 + a4[index]*μ^4 + a5[index]*μ^5
end

function quad_limb_darkening(μ::T, u1::T, u2::T) where T
    μ < zero(T) && return 0.0
    return !iszero(μ) * (one(T) - u1*(one(T)-μ) - u2*(one(T)-μ)^2)
end
