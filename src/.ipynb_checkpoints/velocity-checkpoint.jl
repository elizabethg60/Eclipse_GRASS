function rotation_period(ϕ::T) where T 
    #prescription for differential rotation given a latitude value
    sinϕ = sin(ϕ)
    return 360.0/(14.713 - 2.396*sinϕ^2 - 1.787*sinϕ^4)
end

function v_scalar(lat, lon)
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
    for i in eachindex(A)
        cross_product = cross([0.0,0.0,sun_radius], B[i])
        cross_product ./= norm(cross_product)
        cross_product .*= C[i]
        out[i] = cross_product
    end
    return
end

function projected!(A::Matrix, B::Matrix, out_no_cb::Matrix, out_cb::Matrix, out_cb_new::Matrix, mu::Matrix, cb_velocity::Matrix)
    """
    determine projected velocity of each cell onto line of sight to observer - serial

    A: matrix with xyz and velocity of each cell
    B: matrix with line of sight from each cell to observer
    out: matrix of projected velocities
    """
    for i in 1:length(A)
        #[0.27, 0.325, 12., 0.11]
        #a/(1+exp(−b ∗ (x−c)))+d 0.27/(1+exp(-0.325*(mu[i]-12)))+0.11 
        new_cb_model = 0.27 * (2/(1+exp((mu[i]-0.325)*12))-1) + 0.11
        
        vel = A[i][4:6]
        angle = dot(B[i][1:3], vel) / (norm(B[i][1:3]) * norm(vel))
        
#         vel_old = A[i][4:6] .+ cb_velocity[i]
#         angle_old = dot(B[i][1:3], vel_old) / (norm(B[i][1:3]) * norm(vel_old))
        
#         vel_new = A[i][4:6] .+ new_cb_model
#         angle_new = dot(B[i][1:3], vel_new) / (norm(B[i][1:3]) * norm(vel_new))

        out_no_cb[i] = (norm(vel) * angle)
        out_cb[i] = (norm(vel) * angle) + cb_velocity[i]
        out_cb_new[i] = (norm(vel) * angle) + new_cb_model
    end
    return
end

#updating velocity AND angle gives same results as original with patch vel projected than cb vel added
#diff results found when angle is found using patch vel but the norm vel as cb vel added
    #BUT has worse residuals than without CB SO CB model being included wrong? 

#difference comes from mu - nope makes it super worse 