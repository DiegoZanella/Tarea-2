function PolyModel(n,t,y)
    #Número de observaciones
    m = length(t)

    #Calculamos la matriz de diseño A. Esta es una matriz de  mxn. La entrada A_{ij} = t_i^j
    A = zeros(m,n);
    for i = 1:m
        for j = 1:n
            A[i,j] = t[i]^j
        end
    end


    #Buscamos la factorización de Cholesky de A^T A
    (L,U) = FacChol(A'*A)

    # Encontramos una solución w_0 a la ecuación lineal  Lw = A^Ty mediante la función SolFwd
    #   (porque L es triangular inferior).
    w0 = SolFwd(L,A'*y)

    # Encontramos una solución c_0 a la ecuación lineal L^Tc = w_0 mediante la función SolBwd
    # (porque L^T es triangular superior).
    c0 = SolBwd(L',w0)

    #Aquí construimos la función polinomial PolyMod(x) = c0^T * (1,x, ... , x^n)^T
    function PolyFun(x)
        y = zeros(n+1)
        for i = 0:n
            y[i+1] = x^i
        end
        return c'*y
    end

    return(PolyFun)
end
