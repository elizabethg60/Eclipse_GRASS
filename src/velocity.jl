function rotation_period(ϕ::T) where T 
    #prescription for differential rotation given a latitude value
    sinϕ = sin(ϕ)
    return 360.0/(14.713 - 2.396*sinϕ^2 - 1.787*sinϕ^4)
end

function v_scalar(lat, lon)
    #determines scalar velocity of each cell - serial
    return (2π * sun_radius * cos(lat)) / rotation_period(lat)
end

function pole_vector_grid!(A::Matrix, out::Matrix)
    """
    remove the z component of each cell - serial

    A: matrix of xyz orientation of each cell
    out: matrix of xyz orientation with z removed
    """ 
    for i in 1:length(A)
        out[i] = A[i] - [0.0, 0.0, A[i][3]]
    end
    return
end  

function v_vector(A::Matrix, B::Matrix, C::Matrix, out::Matrix)
    """
    determine velocity vector (direction and magnitude) of each cell - serial 

    A: xyz position of cell
    B: xyz position of cell with z removed
    C: scalar velocity of each cell
    out: matrix with xyz and velocity of each cell
    """
    for i in 1:length(A)
        cross_product = cross(B[i], [0.0,0.0,sun_radius]) 
        cross_product /= norm(cross_product)
        cross_product *= C[i]
        out[i] = [A[i];cross_product]
    end
    return
end

function projected!(A::Matrix, B:: Matrix, out_no_cb::Matrix, out_cb::Matrix, cb_velocity::Matrix) 
    """
    determine projected velocity of each cell onto line of sight to observer - serial

    A: matrix with xyz and velocity of each cell
    B: matrix with line of sight from each cell to observer
    out: matrix of projected velocities
    """
    for i in 1:length(A)
        vel = [A[i][4],A[i][5],A[i][6]]
        angle = cos(π - acos(dot(B[i], vel) / (norm(B[i]) * norm(vel)))) 

        out_no_cb[i] = (norm(vel) * angle) 
        out_cb[i] = (norm(vel) * angle) + cb_velocity[i] 
    end
    return 
end
