using LinearAlgebra
using ForwardDiff


"""
Goal: computes the Jacobian Matrix

#### Input:
    A - data matrix
    t - initial vector
    n - number of observations
    dim - dimension of model
    f - model to fit

#### Output:
    J - jacobian matrix
"""
jacobian(A, f, t, n, dim, J)=
begin
    #J = zeros(n, dim)
    for j = 1 : n
        #calculate r_i
        r(t1) = A[j, 2] - f(A[j, 1], t1)
        for k = 1 : dim
            #find the gradient
            J[j, k] = ForwardDiff.gradient(r, t)[k]
        end
    end
    return J
end


"""
Goal: computes the Residual Matrix

#### Input:
    A - data matrix
    f - model to fit
    t - initial vector
    n - number of observations
    R(optional) - vector to be overscript


#### Output:
    R - residual matrix

"""
residual(A, f, t, n, R)=
begin
    #R = zeros(n)
    for i = 1 : n
        #calculate r_i
        R[i] = A[i, 2] - f(A[i, 1], t)
    end
    return R
end

residual(A, f, t, n)=
begin
    R = zeros(n)
    for i = 1 : n
        #calculate r_i
        R[i] = A[i, 2] - f(A[i, 1], t)
    end
    return R
end


"""
Goal: compute J^T * J + λI

#### Input:
    J - jacobian matrix
    dim - dimension of model
    λ - damping parameter 
    J_1 - matrix to be overscript

#### Output:
    J_1 = J^T * J + λI

"""
auxiliar(J, dim, λ, J_1)=
begin
    J_1 = transpose(J) * J
    for j = 1 : dim
        #add the damping parameter in the main diagonal
        J_1[j, j] = J_1[j, j] + λ
    end
    return J_1
end


"""
Goal: calculate the damping parameter

#### Input:
    A - data matrix
    R - residual matrix
    J - jacobian matrix
    f - model to fit
    t - initial vector
    dp - type of daming parameter(fixed number, DP1, DP2, DP3, DP4, DP5) 
    n - number of observations

#### Output:
    λ - damping parameter
"""
damping(A, R, J, f, t, dp, n)=
begin    
    if typeof(dp) == Float64
        return dp
        
    elseif dp == "DP1"
        return norm(transpose(J) * R) ^ 2 / (0.5 * norm(residual(A, f, t, n)) ^ 2)
        
    elseif dp == "DP2"
        return norm(transpose(J) * R) ^ 2
        
    elseif dp == "DP3"
        return norm(transpose(J) * R)
        
    elseif dp == "DP4"
        return norm(residual(A, f, t, n)) ^ 2
        
    elseif dp == "DP5"
        return norm(residual(A, f, t, n)) 
        
    end
end


"""
Goal: computes an approximate solution to the problem using the Levenberg Maquardt method

#### Input:
    A - data matrix
    f - model function
    t - initial estimation
    dp- type of damping parameter(fixed number, DP1, DP2, DP3, DP4, DP5) 
    ϵ(optional - 1.0e-9) - precision
    itmax(optional - 1000) - maximum of iterations

#### Output:
    t - approximate solution
    k - number of iteractions
"""
levenberg_maquardt(A, f, t, dp, ϵ = 1.0e-10, itmax = 1000)=
begin
    
    n, m = size(A)
    dim = length(t)
    k = 0
    
    R = zeros(n)
    J = zeros(n, dim)
    J_1 = zeros(dim, dim)
    
    J = jacobian(A, f, t, n, dim, J)
    R = residual(A, f, t, n, R)
    
    while norm(transpose(J) * R) > ϵ && k < itmax
        
        #calculate the damping parameter as the user says
        λ = damping(A, R, J, f, t, dp, n)
        #calulate J^T * J + \lambda * I
        J_1 = auxiliar(J, dim, λ, J_1)
        
        #solve the system and find the direction
        d = J_1 \ (- transpose(J) * R)
        
        a = 1.0
        
        #calculate the function
        f1(t1) = 0.5 * norm(residual(A, f, t1, n)) ^ 2
        
        i = 0
        #condição de armijo
        while  f1(t + a * d) > f1(t) + 0.5 * a * transpose(transpose(J) * R) * d && i < itmax
            a = 0.5 * a
        end
        
        #update the data
        t = t + (a * d)
        k = k + 1
        
        J = jacobian(A, f, t, n, dim, J)
        R = residual(A, f, t, n, R)
    end
    return t, k
end

