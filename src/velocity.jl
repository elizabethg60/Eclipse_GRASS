function rotation_period(ϕ::T) where T 
    # prescription for differential rotation given a latitude value
    sinϕ = sin(ϕ)
    return 360.0/(14.713 - 2.396*sinϕ^2 - 1.787*sinϕ^4)
end

function v_scalar(lat, lon)
    return (2π * sun_radius * cos(lat)) / rotation_period(lat)
end

function pole_vector_grid!(A::Matrix, out::Matrix)
    for i in 1:length(A)
        out[i] = A[i] - [0.0, 0.0, A[i][3]]
    end
    return
end  

function v_vector(A::Matrix, B::Matrix, C::Matrix, out::Matrix)
    for i in eachindex(A)
        cross_product = cross([0.0,0.0,sun_radius], B[i])
        cross_product ./= norm(cross_product)
        cross_product .*= C[i]
        out[i] = cross_product
    end
    return
end

function projected!(A::Matrix, B::Matrix, out_no_cb::Matrix, out_cb::Matrix, out_cb_new::Matrix, mu::Matrix, cb_velocity::Matrix)
    for i in 1:length(A)
        new_cb_model = 0.27 * ((2/(1+exp((mu[i]-0.325)*12))) - 1) + 0.11
                       
        vel = A[i][4:6]
        angle = dot(B[i][1:3], vel) / (norm(B[i][1:3]) * norm(vel))

        out_no_cb[i] = (norm(vel) * angle)
        out_cb[i] = (norm(vel) * angle) + cb_velocity[i]
        out_cb_new[i] = (norm(vel) * angle) + new_cb_model
    end
    return
end