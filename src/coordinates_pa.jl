function frame_transfer_pa(A::Matrix, b::Matrix, out::Matrix)
    Threads.@threads for i in 1:length(b)
        out[i] = A*b[i]
    end
    return
end
 
function earth2patch_vectors_pa(A::Matrix, b::Vector, out::Matrix)
    Threads.@threads for i in 1:length(A)	
        out[i] = A[i] .- b
    end
    return
end 

function calc_mu_grid_pa!(A::Matrix, B::Matrix, out::Matrix)
    Threads.@threads for i in 1:length(A)
        out[i] = calc_mu(A[i], B[i])
        end
    return
end	

function matrix_multi(A::Matrix, b::Matrix)
    out = Matrix{Float64}(undef,size(A)...)
    Threads.@threads for i in 1:length(b)
        out[i] = A[i]*b[i]
    end
    return out
end
