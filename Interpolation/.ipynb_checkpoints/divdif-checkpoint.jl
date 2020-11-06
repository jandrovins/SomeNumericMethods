include("simple-gaussian-elimination.jl")
using LinearAlgebra
#
#function diff(A::AbstractVector)
#    return [ (A[i+1] - A[i]) for i = 1:length(A)-1 ]
#end
#

function divdif(X,Y)
    n = length(X)
    D = zeros(n, n)

    D[:,1]=Y'
    for i in 2:n
        aux0 = D[i-1:n, i-1]
        aux = diff(aux0)
        aux2 = X[i:n] - X[1:n-i+1]
        D[i:n, i] = aux./aux2
    end
    show(stdout, "text/plain", diag(D))
end


function newtonform(x, d, xx)
    n = length(d)
    result = d[n]
    for i=n-1:-1:1
        result = result * (xx - x[i]) + d[i]
    end
end
