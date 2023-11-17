function get_xyz(ρ::T, ϕ::T, θ::T) where T
    # pre-compute trig quantitites
    sinϕ, cosϕ = sincos(ϕ)
    sinθ, cosθ = sincos(θ)

    # now get cartesian coords
    x = ρ * cosϕ * cosθ 
    y = ρ * cosϕ * sinθ
    z = ρ * sinϕ

    return [x, y, z]
end

function get_xyz_for_surface(ρ::T; num_lats::Int=100, num_lons::Int=100) where T
    # get grids of polar and azimuthal angles
    ϕ = deg2rad.(range(-90.0, 90.0, length=num_lats)) 
    θ = deg2rad.(range(0.0, 360.0, length=num_lons))'
    return get_xyz.(ρ, ϕ, θ)
end 

function frame_transfer(A::Matrix, b::Matrix, out::Matrix)
    for i in 1:length(b)
        out[i] = A*b[i]
    end
    return
end

function earth2patch_vectors(A::Matrix, b::Vector, out::Matrix)
    for i in 1:length(A)	
        out[i] = (A[i] .- b)
    end
    return 
end 

function calc_mu(SP::Vector, OP::Vector)
    return dot(OP, SP) / (norm(OP) * norm(SP))
end

function calc_mu_grid!(A::Matrix, B::Matrix, out::Matrix)
    for i in 1:length(A)
        out[i] = calc_mu(A[i], B[i])
        end
    return
end	

function calc_dA(radius::T, ϕ::T, dϕ::T, dθ::T) where T
    return radius^2 * sin(π/2.0 - ϕ) * dϕ * dθ
end