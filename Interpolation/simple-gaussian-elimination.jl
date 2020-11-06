function RegressiveSubstitution(A)
    X = Rational[]
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

function PartialPivoting(Ab)
    position = findmax(abs.(Ab[:,1]))[2]
    if position != 1
        Ab[position, :], Ab[1,:] = Ab[1,:], Ab[position, :]
    end
end

function TotalPivoting(Ab, col_moves) # Ab is a @view from complete Ab, created recursively
    cAb = parent(Ab)
    Ab_r = size(Ab)[1] # number of rows and columns of Ab
    cAb_r = size(cAb)[1] # number of rows and columns of complete-Ab
    step = (cAb_r - Ab_r) + 1 # difference in rows equals the number of the step

    A = cAb[step:end,step:end-1] # The sub-matrix within which we will find the largest number
    maxnumber, row_col = findmax(abs.(A))
    row, col = Tuple(row_col) # row and column indices of the largest value
    col_moves[step,:] = [col+step-1 step]
    if col != 1
        cAb[:, col+step-1], cAb[:,step] = cAb[:,step], cAb[:, col+step-1]
    end
    if row != 1
        cAb[row+step-1, :], cAb[step,:] = cAb[step,:], cAb[row+step-1, :]
    end
end

function Escalade(A, pivoting="false", col_moves=nothing)
    if size(A)[1] == 1
        return
    end
    #A = convert(Array{Rational}, A)
    A₁₁ = A[1,1]

    if convert(Float64, A₁₁) == 0. && pivoting=="false" # if there is a 0 on the principal diagonal
        for i in 2:size(A)[1]
            if A[i,1] != 0
                A[1,:], A[i,:] = A[i,:], A[1,:]
                break
            end
        end
    end

    if pivoting == "partial"
        PartialPivoting(A)
    elseif pivoting == "total"
        TotalPivoting(A, col_moves)
    end
    
    A₁₁ = A[1,1]
    for i in 2:size(A)[1]
        if A[i,1] == 0
            continue
        end
        rᵢ = A[i,:]
        factor = A[i,1] / A₁₁
        new_rᵢ= rᵢ - A[1,:] * factor
        A[i,:] = new_rᵢ
    end
    println("STEP $(size(parent(A),1) - size(A, 1) + 1):")
    show(stdout, "text/plain", convert(Array{Float64}, parent(A)))
    println()
    println("END STEP \n\n")

    Escalade(@view(A[2:end,2:end]), pivoting, col_moves)
end

function swap_X(X, col_moves)
    col_moves = reverse(col_moves,dims=1)
    for i in 1:size(col_moves,1)
        i₁,i₂ = col_moves[i, :]
        X[i₁],X[i₂] = X[i₂],X[i₁]
    end
end

function Gauss(A, pivoting="false")
    # Print input matrix
    println("Input matrix:")
    show(stdout, "text/plain", convert(Array{Float64},parent(A)))
    println("\n\n")

    col_moves = nothing
    if pivoting == "total"
        pivotes = size(A,1) - 1 # number of times we have to pivote
        col_moves = Array{Int64}(undef, pivotes, 2)
        Escalade(@view(A[:,:]), pivoting, col_moves)
    else
        Escalade(@view(A[:,:]), pivoting)
    end
    X = RegressiveSubstitution(A[:,:])
    X = convert(Array{Float64}, reverse(X))
    if pivoting == "total"
        println("Columns moves historial: $col_moves")
        swap_X(X, col_moves)
    end
    print("X values: ")
    println(X)
    println("Extended escaladed matrix:")
    show(stdout, "text/plain", convert(Array{Float64},parent(A)))
    println()
    return X
end


#A = [4 2 1;
#     0.25 0.5 1;
#     1 1 0]
#b = [1; 0; 0]
#
#Ab = convert(Array{Rational}, hcat(A,b))
##Ab = hcat(A,b)
#
#X = Gauss(@view(Ab[:,:]), "total")
