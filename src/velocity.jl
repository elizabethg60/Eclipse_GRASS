function lat_grid_fc(num_lats::Int=100, num_lon::Int=100)
    #creates matrix of latitude values reflecting solar grid size - serial 
    ϕ = deg2rad.(range(-90.0, 90.0, length=num_lats))
    A = [ϕ for idx in 1:num_lon]
    return hcat(A...)
end 

function rotation_period(ϕ::T) where T 
    #prescription for differential rotation given a latitude value
    sinϕ = sin(ϕ)
    return 360/(0.9324*(14.713 - 2.396*sinϕ^2 - 1.787*sinϕ^4))
end

function v_scalar!(A:: Matrix, out:: Matrix)
    """
    determines scalar velocity of each cell - serial

    A: matrix of latitudes
    out: matrix of scalar velocity value
    """
	for i in 1:length(A)
		out[i] = (2*π*sun_radius*cos(A[i]))/(rotation_period(A[i])) 
	end
	return
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

function projected!(A::Matrix, B:: Matrix, out_no_cb::Matrix, out_cb::Matrix, cb_velocity::Matrix, epoch, BP_bary)
    """
    determine projected velocity of each cell onto line of sight to observer - serial

    A: matrix with xyz and velocity of each cell
    B: matrix with line of sight from each cell to observer
    out: matrix of projected velocities
    """
    for i in 1:length(A)
        vel = [A[i][4],A[i][5],A[i][6]] 
        angle = cos(π - acos(dot(B[i], vel) / (norm(B[i]) * norm(vel)))) 
        out_no_cb[i] = (norm(vel) * angle)  #+ dot(spkssb(10,epoch,"J2000")[1:3], spkssb(399,epoch,"J2000")[4:6])/86.4 
        out_cb[i] = (norm(vel) * angle) + cb_velocity[i] 
    end
    return 
end
