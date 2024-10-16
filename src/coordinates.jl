function get_grid_centers(grid::StepRangeLen)
    #determine centers of grid lat / long vectors
    start = first(grid) + 0.5 * step(grid)
    stop = last(grid) - 0.5 * step(grid)
    return range(start, stop, length=length(grid)-1)
end

function get_Nθ(ϕc, dϕ)
    #determine number of lons for each lat slice
    return ceil(Int, 2π * cos(ϕc) / dϕ)
end

function calc_mu(SP::AbstractArray{T,1}, OP::AbstractArray{T,1}) where T
    #determine mu value for each cell
    return dot(OP, SP) / (norm(OP) * norm(SP))
end

function calc_mu_grid!(A::Matrix, B::Matrix, out::Matrix)
    """
    create matrix of mu values for each cell - serial

    A: matrix of vectors from sun center to cell
    B: matrix of vectors from observer to cell
    out: mu value between A and B
    """
    for i in eachindex(A)
        out[i] = calc_mu(view(A[i], 1:3), view(B[i], 1:3))
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
