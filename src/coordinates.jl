function get_grid_centers(grid::StepRangeLen)
    # determine centers of grid lat / long vectors
    start = first(grid) + 0.5 * step(grid)
    stop = last(grid) - 0.5 * step(grid)
    return range(start, stop, length=length(grid)-1)
end

function get_Nθ(ϕc, dϕ)
    # determine number of lons for each lat slice
    return ceil(Int, 2π * cos(ϕc) / dϕ)
end

function calc_mu(SP::AbstractArray{T,1}, OP::AbstractArray{T,1}) where T
    # determine mu value for each cell
    return dot(OP, SP) / (norm(OP) * norm(SP))
end

function calc_mu_grid!(A::Matrix, B::Matrix, out::Matrix)
    for i in eachindex(A)
        out[i] = calc_mu(view(A[i], 1:3), view(B[i], 1:3))
        end
    return
end

function calc_dA(radius, ϕ, dϕ, dθ)
    return radius^2 * sin(π/2.0 - ϕ) * dϕ * dθ
end
