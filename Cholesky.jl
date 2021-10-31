using LinearAlgebra

# https://algowiki-project.org/en/Cholesky_decomposition
#      Calcula la factorización de Cholesky de una matriz simétrica
# definida positiva A.
#Le das como input A y te regresa (L,L^T)
function FacChol(A)
    (n,m) = size(A)
    L = Matrix(1.0*I,n,n)
    # Primero definimos la columna 1 de L
    L[1,1] = sqrt(A[1,1])
    for i = 2:n
        L[i,1] = A[i,1]/L[1,1]
    end
    # Ahora definimos el resto de columnas de L
    for j = 2:n
        for i = j:n
            if i == j #caso en la diagonal
                L[i,j] = sqrt(abs( A[i,j]- dot( L[j,1:j-1], L[j, 1:j-1] ) ))
            else #caso debajo de la diagonal
                L[i,j] = ( A[i,j]-dot( L[i,1:j-1],L[j,1:j-1]) )/L[j,j]
            end
        end
    end
    (L,L')
end
