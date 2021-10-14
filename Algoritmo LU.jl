using LinearAlgebra

#Matriz de eliminación
#La función M(A,i) toma como input una matriz A de nxn y un entero i.
#M(A,i) regresa la matriz de eliminación de A haciendo pivote en A[i,i]
function M(A,i0)
    (n,m) = size(A)
    if n != m
       throw(DimensionMismatch)
    elseif i0 >=n
    else
      M = Matrix(1.0*I,n,n)  #Pivote es A[i0,i0]
      for j = i0+1 : n
          M[j,i0] = - A[j, i0]/A[i0,i0]
      end
      M
    end
 end

 # Inversa de la matriz de eliminación
function MInv(M)
   (n,m) = size(A)
   L = Matrix(1.0I,n,n)
   for j = i0+1 : n
       M[j,i0] =  A[j, i0] / A[i0,i0]
   end
end


# Matriz de permutación. Permuta los renglones 1 y j de la matriz identidad
# de nxn
function P(n,j,l)
   P = Matrix(1.0*I,n,n)
   (P[l,l], P[l,j],  P[j,j], P[j,l]) = (0.0, 1.0, 0.0, 1.0)
   P
end


#Solución hacia adelante de una ecuación lineal triangular inferior
function SolFwd(A,b)
   #A es una matriz de nxn triangular inferior sin ceros en la diagonal
   (n,m) = size(A)
   (C,d) = (A,b)
   x = zeros(n)
   x[1] = d[1]/C[1,1]
   for i = 2:n
      s = 0.0
      for j = 1:(i-1)
         s = s + C[i,j]*x[j]
      end
      x[i] = ( d[i]-s )/C[i,i]
   end
   x
end

# Solución hacia atrás de una ecuación lineal triangular superior
function SolBwd(A,b)
   (n,m) = size(A)
   (C,d) = (A,b)
   x = zeros(n)
   x[n] = d[n]/C[n,n]
   for j = 2:n
      s = 0
      for k = (n+1)-j+1:n
         s = s+C[n+1-j,k]x[k]
      end
      x[n+1-j] =  (d[n+1-j]-s)/C[n+1-j,n+1-j]
      end
   x
   end

# Método LU sin pivoteo parcial.
function ALU(A,b)
      (n,m) = size(A)
      if n != m
         throw(DimensionMismatch)
      else
         C = A
         d = b
         for i = 1:n-1
            d = M(C,i)*d
            C = M(C,i)*C
         end
         (C,d)
      end
end

#Método LU con pivoteo parcial
function LUPP(A,b)
   (n,m) = size(A)
   (U,d) = (A,b)
   #perm = zeros(Int8,1,2)
   for i = 1: n-1
      k = findfirst(x -> x == maximum(abs.(U[i:n,i])),abs.(U[i:n,i]))+(i-1)
      #perm[1,i] = k
      d = P(n,k,i)*d
      U = P(n,k,i)*U
      d = M(U,i)*d
      U = M(U,i)*U
   end
   (U,d)
end
