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

function LUpartial(A,b)
    n=size(A,1)
    L=Matrix{Rational}(I, n, n)
    U=zeros(n,n)
    P=Matrix{Rational}(I, n, n)
    M=copy(A)

    # Factorization
    for i in 1:n-1
        # Change of rows
        aux0, aux = findmax(abs.(M[i+1:n,i]))

        if aux0>abs(M[i,i])
            aux2=M[i+aux,i:n];
            aux3=P[i+aux,:];
            M[aux+i,i:n]=M[i,i:n];
            P[aux+i,:]=P[i,:];
            M[i,i:n]=aux2;
            P[i,:]=aux3;
            if i>1
               aux4=L[i+aux,1:i-1];
               L[i+aux,1:i-1]=L[i,1:i-1];
               L[i,1:i-1]=aux4;
            end
        end

        for j in i+1:n
            if M[j,i] != 0
                L[j,i] = M[j,i]/M[i,i]
                M[j,i:n]=M[j,i:n] - (M[j,i]/M[i,i])*M[i,i:n]
            end
        end
        U[i,i:n]=M[i,i:n]
        U[i+1,i+1:n]=M[i+1,i+1:n]
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
x = LUpartial(A, b)

println(x)
