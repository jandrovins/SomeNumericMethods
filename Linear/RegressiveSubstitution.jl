module RegressiveSubstitution
function RegressiveSubstitution(A)
    X = Array{Float64}[]
    for i in size(A)[1]:-1:1
        j = size(A)[1] - i
        rᵢ= A[i,:]
        bᵢ= rᵢ[end]
        Cₓ= rᵢ[end-1]
        for k in 1:j
            bᵢ-= rᵢ[end-k]
            Cₓ=rᵢ[end-k-1]
        end
        xᵢ= bᵢ/Cₓ
        append!(X, xᵢ)
        A[1:i-1,end-1-j] = A[1:i-1,end-1-j] .* xᵢ
    end
    return X
end
