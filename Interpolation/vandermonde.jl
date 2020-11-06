include("simple-gaussian-elimination.jl")

function vandermonde(X, Y)
    n=size(X, 1)
    A=zeros(n,n)

    for i in 1:n
        A[:,i]=(X.^(n-i))
    end
        
    Ab = hcat(A, Y)
    Coef = Gauss(@view(Ab[:,:]), "total")
end

#x = [1; 2; 3; 4; 5; 6; 7]
#y = [1.1247; -0.8540; 0.5864; 9; -0.9062; 0.9081; -0.2700]

x = [-2; -1; 2; 3]
y = [12.13533528; 6.367879441; -4.610943901; 2.085536923]

#vandermonde(x,y)
