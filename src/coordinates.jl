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

function get_grid_centers(grid::StepRangeLen)
    start = first(grid) + 0.5 * step(grid)
    stop = last(grid) - 0.5 * step(grid)
    return range(start, stop, length=length(grid)-1)
end
function get_Nθ(ϕc, dϕ)
    return ceil(Int, 2π * cos(ϕc) / dϕ)
end

# function get_xyz_for_surface(ρ::T; num_lats::Int=100, num_lons::Int=100) where T
#     """
#     get grid of polar and azimuthal angles

#     ρ: object radius
#     num_lats: number of latitude for grid
#     num_lons: number of longitude for grid
#     """
#     ϕ = deg2rad.(range(-90.0, 90.0, length=num_lats)) 
#     ϕe = range(deg2rad(-90.0), deg2rad(90.0), length=num_lats)
#     ϕc = get_grid_centers(ϕe)

#     θ = deg2rad.(range(0.0, 360.0, length=num_lons))'
#     θe = range(deg2rad(0.0), deg2rad(360.0), length=num_lons)
#     θc = get_grid_centers(θe)'
#     return get_xyz.(ρ, ϕc, θc)
# end 

function frame_transfer(A::Matrix, b::Matrix, out::Matrix)
    """
    transforms between frames - serial
    
    A: transformation matrix
    b: matrix in initial frame
    out: matrix in desired frame 
    """
    for i in 1:length(b)
        out[i] = A*b[i]
    end
    return
end

function earth2patch_vectors(A::Matrix, b::Vector, out::Matrix)  
    """
    determines line of sight vector from observer to each cell on grid - serial

    A: matrix of vectors of barycenter to cells
    b: vector from barycenter to observer
    out: matrix of vectors of observer to each patch 
    """
    for i in 1:length(A)	
        out[i] = A[i] .- b
    end
    return 
end 

function calc_mu(SP::Vector, OP::Vector)
    #determine mu value for each cell
    return dot(OP, SP) / (norm(OP) * norm(SP)) 
end

function calc_mu_grid!(A::Matrix, B::Vector, out::Matrix)
    """
    create matrix of mu values for each cell - serial

    A: matrix of vectors from sun center to cell
    B: matrix of vectors from observer to cell
    out: mu value between A and B
    """
    for i in 1:length(A)
        out[i] = calc_mu(A[i], B)
        end
    return
end	

function calc_dA(radius, ϕ, dϕ, dθ)
    """
    determines area of each cell

    radius: radius of object
    ϕ: cell's latitude value
    dϕ: step size of latitude
    dθ: step size of longitude
    """
    return radius^2 * sin(π/2.0 - ϕ) * dϕ * dθ
end