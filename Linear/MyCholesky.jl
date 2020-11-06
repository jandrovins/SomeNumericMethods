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

function MyCholesky(A, b)
    n=size(A,1)
    L=Matrix{Rational}(I, n, n)
    U=Matrix{Rational}(I, n, n)
    # Factorización
    for i in 1:n-1
        L[i,i]=sqrt(A[i,i]-dot(L[i,1:i-1],U[1:i-1,i]'))
        U[i,i]=L[i,i]
        for j in i+1:n
            L[j,i]=(A[j,i]-dot(L[j,1:i-1],U[1:i-1,i]'))/U[i,i];
        end
        for j in i+1:n
            U[i,j]=(A[i,j]-dot(L[i,1:i-1],U[1:i-1,j]'))/L[i,i];
        end
    end
    L(n,n)=sqrt(A[n,n]-dot(L[n,1:n-1],U[1:n-1,n]'));
    U[n,n]=L[n,n];

    # Entrega de resultados
    z=progsust(hcat(L, b))
    x=regrsust(hcat(U, z))  
end

A = [36 3 -4 5;
     5 -45 10 -2;
     6 8 57 5;
     2 3 -8 -42]

b = [-20; 69; 96; -32]

A = convert(Array{Rational}, A)
b = convert(Array{Rational}, b)
x = MyCholesky(A,b)

println(x)
