using LinearAlgebra

function Householder(A)
    (m,n)=size(A)

    x = A[:,1]

    e1 = zeros(m)
    e1[1] = 1.0

    v = x+sign(A[1,1])*norm(x)*e1
    v = v/norm(v)

    H = I -2*v*v'
end

#Esta función toma como argumento una matriz A de mxn. Regresa una tupla (Q,R) donde Q es una matriz ortogonal y
# R es una matriz triangular superior tales que QR = A.
function QRHouseholder(A)
    (m,n) = size(A)
    Q = Matrix(1.0*I,m,m)
    R = A

    for i = 1:n
        HTild = Householder(R[i:m,i:n] )

        H = Matrix(1.0*I,m,m)
        H[i:m,i:m] = HTild
    
        R = H*R
        Q = H*Q
    end
    Q = Q'
    return((Q,R))
end

A = [1 2; 3 4; 5 6]

(Q,R) = QRHouseholder(A)
