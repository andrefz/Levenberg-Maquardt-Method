using DelimitedFiles
using Plots
using Statistics

## TESTES UNITÁRIOS
struct sample_type
    name::String
    data::Array{Float64,2}
    n::Int64
    dim::Int64
    model::Function
    solution::Vector{Float64}
end


"""
Goal: reads the file and separates the information

#### Input:
    file

#### Output:
    name of the problem
    data matrix
    number of measurements
    dimension of the problem
    model function
    real solution
"""
sample_parsing(file::String)=
begin
    data_sample = readdlm(file, ':')
    return sample_type(
                      data_sample[1,2],                   # name
                      eval(Meta.parse(data_sample[2,2])), # matrix
                      data_sample[3,2],                   # num pts
                      data_sample[6,2],                   # num parameters
                      eval(Meta.parse(data_sample[5,2])), # function
                      eval(Meta.parse(data_sample[9,2])), # solution
                      )
end


"""
Goal: reads the file and separates the information

#### Input:
    file

#### Output:
    name of the problem
    data matrix
    number of measurements
    dimension of the problem
    model function
    real solution
"""
plotting(A, aprox_solution, title1::String)=
begin
    aprox_matrix = copy(A)

    for i = 1 : n
        aprox_matrix[i, 2] = func(A[i, 1], aprox_solution)
    end

    #sort using the first column
    aprox_matrix = aprox_matrix[sortperm(aprox_matrix[:, 1]), :]
    
    #plot the points using the actual solution
    scatter(A[:, 1], A[:, 2], title = title1, color = "lightgreen", label = "Dados de entrada")
    #plot the points using the found solution by the method
    plot!(aprox_matrix[:, 1], aprox_matrix[:, 2], color = "black", label = "Aproximação")
  
end


## TESTES GERAIS


"""
Goal: avaliate the method 

#### Input:
    model type - model_string(parabola|cubic|gaussian|log)
    initial_value - initial value
    dp - type of damping parameter(fixed number, DP1, DP2, DP3, DP4, DP5)
    num_files(optional) - the number of files to be used. If not given, it will be used all files os the problem

#### Output:
    k - the medium number of iteractions using the damping parameter given
    e - the medium error between the solution found and the solution given using the damping parameter given
"""
allfiles(model_string, initial_value, dp)=
begin
    
    #directory path
    dir_path = "/home/julia/Documents/4° ano/Otimização Não Linear/problems/"*model_string*"/"
    filelist = readdir(dir_path)

    #number of files
    num_files = size(filelist)[1]

    #counter of iteractions
    iteractions = 0

    k = 0
    e = 0
    for i = 1 : num_files 
        sample_test = sample_parsing(dir_path*filelist[i])

        name = sample_test.name
        A = sample_test.data
        n = sample_test.n
        dim = sample_test.dim
        func = sample_test.model
        sol = sample_test.solution
       
        #covert the type of f
        f(x, t) = Base.invokelatest(func, x, t)
        
        #M = cleaning(A, sol, func)   
        #calculate the solution by the method
        e1, k1 = levenberg_maquardt(A, f, initial_value, dp)
        
        
        k = k + k1 #number of iteractions
        e = e + norm(e1 - sol) #error
        
    end
    
    k = k / num_files
    e = e / num_files
   
    

    println("Número médio de iterações o parâmetro: ", k)
    println("Valor médio do erro para o parâmetro: ", e)
end
allfiles(model_string, initial_value, dp, num_files)=
begin
   
    dir_path = "./problems/"*model_string*"/"
    filelist = readdir(dir_path)

    #num_files = size(filelist)[1]

    iteractions = 0

    k = 0
    e = 0
    for i = 1 : num_files 
        sample_test = sample_parsing(dir_path*filelist[i])

        name = sample_test.name
        A = sample_test.data
        n = sample_test.n
        dim = sample_test.dim
        func = sample_test.model
        sol = sample_test.solution
       
        f(x, t) = Base.invokelatest(func, x, t)
          
        e1, k1 = levenberg_maquardt(A, f, initial_value, dp)
        
        k = k + k1

        e = e + norm(e1 - sol)
        
    end
    
    k = k / num_files
    e = e / num_files
   
    

    println("Número médio de iterações o cada parâmetro: ", k)
    println("Valor médio do erro para o parâmetro: ", e)
end


"""
Goal: clean the perturbations of the data matrix

#### Input:
    A - data matrix
    t - initial vector
    f - model to fit

#### Output:
    M - data matrix with no perturbations
"""
cleaning(A, t, f)=
begin
    k = 0
    n, m = size(A)
    for i = 1 : n
        #find the amount of information without the perturbation
        if abs(A[i, 2] - f(A[i, 1], t)) < 1.0e-8 && abs(A[i, 1]) > 1.0e-8  && abs(A[i, 2]) > 1.0e-8
            k = k + 1
        end
    end
    
    M = zeros(k, m)
    d = 1
    for j = 1: n 
        #put this information in a new data matrix
        if abs(A[j, 2] - f(A[j, 1], t)) < 1.0e-8 && abs(A[j, 1]) > 1.0e-8  && abs(A[j, 2]) > 1.0e-8
            M[d, :] = A[j, :]
            d = d + 1        
        end
    end
    
    return M
end
