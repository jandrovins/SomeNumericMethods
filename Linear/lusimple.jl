using LinearAlgebra

function progsust(M)
    n = size(M,1)
    x = zeros(n,1)

    x[1]=M[1,n+1]/M[1,1]
    
    for i in 2:n
        aux=hcat(1, x[1:i-1]')
        aux1=hcat(M[i,n+1], -M[i,1:i-1]')
        x[i]=dot(aux, aux1)/M[i,i]
    end
    return x
end

function regrsust(M)
    n=size(M,1)
    x=zeros(n,1)

    x[n]=M[n,n+1]/M[n,n]
    for i in n-1:-1:1
        aux=hcat(1, x[i+1:n]')
        aux1=hcat(M[i,n+1], -M[i,i+1:n]')
        x[i]=dot(aux, aux1)/M[i,i]
    end
    return x
end

function LUsimple(A,b)
    n=size(A,1)
    L=Matrix{Rational}(I, n, n)
    U=zeros(n,n)
    M=copy(A)

    for i in 1:n-1
        for j in i+1:n
            if M[j,i] != 0
                L[j,i] = M[j,i]/M[i,i]
                M[j,i:end]=M[j,i:end] - (M[j,i]/M[i,i])*M[i,i:end]
            end
        end
        U[i,i:end]=M[i,i:end]
        U[i+1,i+1:end]=M[i+1,i+1:end]
    end
    U[n,n]=M[n,n]

    z = progsust(hcat(L, b))
    x = regrsust(hcat(U, z))
    return x
end

A = [4 3 -2 -7;
     3 12 8 -3;
     2 3 -9 3;
     1 -2 -5 6]

b = [20; 18; 31; 12]

A = convert(Array{Rational}, A)
b = convert(Array{Rational}, b)
x = LUsimple(A, b)

println(x)
