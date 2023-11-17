function lat_grid_fc(num_lats::Int=100, num_lon::Int=100)
    ϕ = deg2rad.(range(-90.0, 90.0, length=num_lats))
    A = [ϕ for idx in 1:num_lon]
    return hcat(A...)
end 

function rotation_period(ϕ::T; A::T=14.713, B::T=-2.396, C::T=-1.787) where T 
    sinϕ = sin(ϕ)
    return 360/(14.713 - 2.396*sinϕ^2 - 1.787*sinϕ^4)
end

function v_scalar!(A:: Matrix, out:: Matrix)
	for i in 1:length(A)
		out[i] = (2*π*sun_radius*cos(A[i]))/(rotation_period(A[i]))
	end
	return
end

function pole_vector_grid!(A::Matrix, b::Vector, out::Matrix)
    for i in 1:length(A)
        out[i] = A[i] - [0.0, 0.0, A[i][3]]
    end
    return
end  

function v_vector(A::Matrix, B::Matrix, C::Matrix, out::Matrix)
    for i in 1:length(A)
        cross_product = cross(B[i], [0.0,0.0,sun_radius]) 
        cross_product /= norm(cross_product)
        cross_product *= C[i]
        out[i] = [A[i];cross_product]
    end
    return
end

function projected!(A::Matrix, B:: Matrix, out::Matrix, cb_velocity::Matrix)
    for i in 1:length(A)
        vel = [A[i][4],A[i][5],A[i][6]]
        angle = dot(B[i], vel) / (norm(B[i]) * norm(vel))
        out[i] = -(norm(vel) * angle) + cb_velocity[i] 
    end
    return 
end