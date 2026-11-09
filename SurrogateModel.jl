"""
In this code will there be codes for creating a surrgoate model descibing HTO on the anode side.
    It is going to use data from the electrolyzer model to create a surrogate model using polynomial
    regression using JuMP and Ipopt. The surrogate model will be used to predict HTO based 
    on the power input, temperature, and pressure to the electrolyzer.

"""

using JuMP
using Ipopt
using CSV
using DataFrames

include("SaveSurrogateModel.jl")


function scaled(x)
    return (x .- ((maximum(x)+minimum(x))/2))./(maximum(x)-minimum(x))/2
end

function scaled_minmax(x, xmin, xmax)
    return (x .- xmin)./(xmax - xmin)
end

function scaled_maxAbs(x)
    return x ./ maximum(abs.(x))
end



# data = CSV.read(joinpath("/home/annkar/NewCode/Data", "surrogate_data.csv"), DataFrame)
# 
# data_validation = CSV.read(joinpath("/home/annkar/NewCode/Data", "surrogate_data_Pc_validation.csv"), DataFrame)

function fit_surrogate(;csv_path = joinpath("/home/annkar/NewCode/Data", "surrogate_data.csv"))
    # gets the data from the CSV file
    data = CSV.read(csv_path, DataFrame)

    ##=====================================================================
    #defining different optimizer models
    ##=====================================================================
     
    
        #Linear
    model_linear = Model(Ipopt.Optimizer)
    set_attribute(model_linear, "linear_solver", "mumps")
        #Quadratic
    model_quadratic = Model(Ipopt.Optimizer)
    set_attribute(model_quadratic, "linear_solver", "mumps")

        #Hyperbolic
    model_hyperbolic_1 = Model(Ipopt.Optimizer)
    set_attribute(model_hyperbolic_1, "linear_solver", "mumps")
    
    model_hyperbolic_2 = Model(Ipopt.Optimizer)
    set_attribute(model_hyperbolic_2, "linear_solver", "mumps")

    model_hyperbolic_3 = Model(Ipopt.Optimizer)
    set_attribute(model_hyperbolic_3, "linear_solver", "mumps")

    model_hyperbolic_4 = Model(Ipopt.Optimizer)
    set_attribute(model_hyperbolic_4, "linear_solver", "mumps")
    
        #Logarithmic
    model_logarithmic = Model(Ipopt.Optimizer)
    set_attribute(model_logarithmic, "linear_solver", "mumps")

        #Combination of logarithmic and hyperbolic
    model_log_hyper1 = Model(Ipopt.Optimizer)
    set_attribute(model_log_hyper1, "linear_solver", "mumps")

    model_log_hyper2 = Model(Ipopt.Optimizer)
    set_attribute(model_log_hyper2, "linear_solver", "mumps")

    model_log_hyper3 = Model(Ipopt.Optimizer)
    set_attribute(model_log_hyper3, "linear_solver", "mumps")


    ##=====================================================================
    # extracts the values values from csv data
    ##=====================================================================
    Pw = data.Power
    HTO_exp = data.HTO
    N = length(Pw)

    Pw_scaled = scaled_minmax(Pw, minimum(Pw), maximum(Pw))


    ##=====================================================================
    #defining variables for the differet polynomial regression models
    ##=====================================================================
        #linear regression model
    # @variable(model_linear, beta00_L )
    # @variable(model_linear, beta01_L , start=1.0)
    # @variable(model_linear, HTO_L[1:N])
    # #create the constraints for the linear regression model
    # @constraint(model_linear, [i=1:N], beta00_L + beta01_L*Pw_scaled[i]== HTO_L[i])
    # # defining the objective function for the linear regression model
    # @objective(model_linear, Min, sum((HTO_L[i] - HTO_exp[i])^2 for i in 1:N))


        #quadratic regression model
    @variable(model_quadratic, beta00_Q )
    @variable(model_quadratic, beta01_Q )
    @variable(model_quadratic, beta02_Q )
    @variable(model_quadratic, HTO_Q[1:N]) 
    #create the constraints for the polynomial regression model
    @constraint(model_quadratic, [i=1:N], beta00_Q + beta01_Q*Pw_scaled[i] + beta02_Q*Pw_scaled[i]^2 == HTO_Q[i])
    # defining the objective function for the polynomial regression model
    @objective(model_quadratic, Min, sum((HTO_Q[i] - HTO_exp[i])^2 for i in 1:N))

        #hyperbolic regression model
    
        #First version of hyperbolic regression model
    @variable(model_hyperbolic_1, beta00_H1)
    @variable(model_hyperbolic_1, beta01_H1 )
    @variable(model_hyperbolic_1, beta02_H1 ) 
    @variable(model_hyperbolic_1, c00_H1 >= 0.0)
    @variable(model_hyperbolic_1, HTO_H1[1:N]) 
    #create the constraints for the polynomial regression model
    @constraint(model_hyperbolic_1, [i=1:N], beta00_H1 + beta01_H1*Pw_scaled[i] + beta02_H1*(1/(Pw_scaled[i] + c00_H1)) == HTO_H1[i])
    # defining the objective function for the polynomial regression model
    @objective(model_hyperbolic_1, Min, sum((HTO_H1[i] - HTO_exp[i])^2 for i in 1:N))
    
    
        #Second version of hyperbolic regression model
    # @variable(model_hyperbolic_2, beta00_H2)
    @variable(model_hyperbolic_2, beta01_H2 )
    @variable(model_hyperbolic_2, beta02_H2 )
    @variable(model_hyperbolic_2, c00_H2 >= 0.0)
    @variable(model_hyperbolic_2, HTO_H2[1:N]) 
    #create the constraints for the polynomial regression model
    @constraint(model_hyperbolic_2, [i=1:N], beta01_H2*Pw_scaled[i] + beta02_H2*(1/(Pw_scaled[i]+c00_H2)) == HTO_H2[i])
    # defining the objective function for the polynomial regression model
    @objective(model_hyperbolic_2, Min, sum((HTO_H2[i] - HTO_exp[i])^2 for i in 1:N))



        #third version of hyperbolic regression model
    @variable(model_hyperbolic_3, beta00_H3)
    @variable(model_hyperbolic_3, beta01_H3 )
    @variable(model_hyperbolic_3, beta02_H3 )
    @variable(model_hyperbolic_3, beta03_H3 )
    @variable(model_hyperbolic_3, c00_H3 >= 0.0)
    @variable(model_hyperbolic_3, HTO_H3[1:N]) 
    #create the constraints for the polynomial regression model
    @constraint(model_hyperbolic_3, [i=1:N], beta00_H3 + beta01_H3*Pw_scaled[i] + beta02_H3*(1/(beta03_H3*Pw_scaled[i]+c00_H3)) == HTO_H3[i])
    # defining the objective function for the polynomial regression model
    @objective(model_hyperbolic_3, Min, sum((HTO_H3[i] - HTO_exp[i])^2 for i in 1:N))

        #fourth version of hyperbolic regression model
    @variable(model_hyperbolic_4, beta01_H4 )
    @variable(model_hyperbolic_4, beta02_H4 )
    @variable(model_hyperbolic_4, beta03_H4, start=1.0)
    @variable(model_hyperbolic_4, c00_H4 >= 0.0)
    @variable(model_hyperbolic_4, HTO_H4[1:N]) 
    #create the constraints for the polynomial regression model
    @constraint(model_hyperbolic_4, [i=1:N], beta01_H4*Pw_scaled[i] + beta02_H4*(1/(beta03_H4*Pw_scaled[i]+c00_H4)) == HTO_H4[i])
    # defining the objective function for the polynomial regression model
    @objective(model_hyperbolic_4, Min, sum((HTO_H4[i] - HTO_exp[i])^2 for i in 1:N))

        #logarithmic regression model
    @variable(model_logarithmic, beta00_LG)
    @variable(model_logarithmic, beta01_LG)
    @variable(model_logarithmic, beta02_LG)
    @variable(model_logarithmic, c00_LG >= 0.0)
    @variable(model_logarithmic, c01_LG >= 0.0)
    @variable(model_logarithmic, HTO_LG[1:N])
    #create the constraints for the logarithmic regression model
    @constraint(model_logarithmic, [i=1:N], beta00_LG  + beta01_LG*log(Pw_scaled[i] + c00_LG) + beta02_LG*log(Pw_scaled[i] + c01_LG)^2 == HTO_LG[i])
    # defining the objective function for the logarithmic regression model
    @objective(model_logarithmic, Min, sum((HTO_LG[i] - HTO_exp[i])^2 for i in 1:N))

        #combination of logarithmic and hyperbolic regression model

    @variable(model_log_hyper1, beta00_LH1)
    @variable(model_log_hyper1, beta01_LH1)
    @variable(model_log_hyper1, beta02_LH1)
    @variable(model_log_hyper1, beta03_LH1)
    @variable(model_log_hyper1, c01_LH1 >= 0.0)
    @variable(model_log_hyper1, c00_LH1 >= 0.0)
    @variable(model_log_hyper1, c02_LH1 >= 0.0)
    @variable(model_log_hyper1, HTO_LH1[1:N])
    #create the constraints for the combination model
    @constraint(model_log_hyper1, [i=1:N], beta00_LH1 + beta01_LH1*log(Pw_scaled[i] + c00_LH1) + beta02_LH1*log(Pw_scaled[i] + c01_LH1)^2 + beta03_LH1*(1/(Pw_scaled[i]+c02_LH1)) == HTO_LH1[i])
    # defining the objective function for the combination model
    @objective(model_log_hyper1, Min, sum((HTO_LH1[i] - HTO_exp[i])^2 for i in 1:N))

    @variable(model_log_hyper2, beta00_LH2)
    @variable(model_log_hyper2, beta01_LH2)
    @variable(model_log_hyper2, beta02_LH2)
    @variable(model_log_hyper2, c00_LH2 >= 0.0)
    @variable(model_log_hyper2, c01_LH2 >= 0.0)
    @variable(model_log_hyper2, HTO_LH2[1:N])
    #create the constraints for the combination model
    @constraint(model_log_hyper2, [i=1:N], beta00_LH2 + beta01_LH2*log(Pw_scaled[i] + c00_LH2) + beta02_LH2*(1/(Pw_scaled[i]+c01_LH2)) == HTO_LH2[i])
    # defining the objective function for the combination model
    @objective(model_log_hyper2, Min, sum((HTO_LH2[i] - HTO_exp[i])^2 for i in 1:N))

    @variable(model_log_hyper3, beta00_LH3)
    @variable(model_log_hyper3, beta01_LH3)
    @variable(model_log_hyper3, beta02_LH3)
    @variable(model_log_hyper3, beta03_LH3)
    @variable(model_log_hyper3, beta04_LH3)
    @variable(model_log_hyper3, c00_LH3 >= 0.0)
    @variable(model_log_hyper3, c01_LH3 >= 0.0)
    @variable(model_log_hyper3, c02_LH3 >= 0.0)
    @variable(model_log_hyper3, HTO_LH3[1:N])
    #create the constraints for the combination model
    @constraint(model_log_hyper3, [i=1:N], beta00_LH3 + beta01_LH3*log(Pw_scaled[i] + c00_LH3) + beta02_LH3*log(Pw_scaled[i] + c01_LH3)^2 + beta03_LH3*(1/(beta04_LH3*Pw_scaled[i]+c02_LH3)) == HTO_LH3[i])
    
    # defining the objective function for the combination model
    @objective(model_log_hyper3, Min, sum((HTO_LH3[i] - HTO_exp[i])^2 for i in 1:N))


    ##=====================================================================
    #solving the optimization probelm
    ##=====================================================================
    # println("Optimizing the linear model...")
    # optimize!(model_linear)
    
    println("Optimizing the quadratic model...")
    optimize!(model_quadratic)

    println("Optimizing the hyperbolic model...")
    optimize!(model_hyperbolic_1)   

    println("Optimizing the hyperbolic model...")
    optimize!(model_hyperbolic_2)

    println("Optimizing the hyperbolic model...")
    optimize!(model_hyperbolic_3)

    println("Optimizing the hyperbolic model...")
    optimize!(model_hyperbolic_4)

    println("Optimizing the logarithmic model...")
    optimize!(model_logarithmic)

    println("Optimizing the combination model...")
    optimize!(model_log_hyper1)

    println("Optimizing the combination model...")
    optimize!(model_log_hyper2)
    
    println("Optimizing the combination model...")
    optimize!(model_log_hyper3)





    # println("Solver status for linear model: ", termination_status(model_linear))
    println("Solver status for quadratic model: ", termination_status(model_quadratic))
    println("Solver status for hyperbolic model 1: ", termination_status(model_hyperbolic_1))
    println("Solver status for hyperbolic model 2: ", termination_status(model_hyperbolic_2))
    println("Solver status for hyperbolic model 3: ", termination_status(model_hyperbolic_3))
    println("Solver status for hyperbolic model 4: ", termination_status(model_hyperbolic_4))
    println("Solver status for logarithmic model: ", termination_status(model_logarithmic))
    println("Solver status for combination model 1: ", termination_status(model_log_hyper1))
    println("Solver status for combination model 2: ", termination_status(model_log_hyper2))
    println("Solver status for combination model 3: ", termination_status(model_log_hyper3))

    # if is_solved_and_feasible(model_linear)
    #     println("Objective value for linear model: ", objective_value(model_linear))
    #     beta_lin = [value(beta00_L), value(beta01_L)]
    #     SaveModel_oneInput("Power_linear", beta_lin, [], objective_value(model_linear))
    #     # println("\nLinear model:")
    #     # println("beta0 = ", value(beta00_L))
    #     # println("beta1 = ", value(beta01_L))
    #     # println("beta2 = ", value(beta02_L))

    # end

    if is_solved_and_feasible(model_quadratic)
        println("Objective value for quadratic model: ", objective_value(model_quadratic))
        beta_quad = [value(beta00_Q), value(beta01_Q), value(beta02_Q)]
        SaveModel_oneInput("Power_quadratic", beta_quad, [], objective_value(model_quadratic))
        # println("\nQuadratic model:")
        # println("beta0 = ", value(beta00_Q))
        # println("beta1 = ", value(beta01_Q))
        # println("beta2 = ", value(beta02_Q))

    end

    if is_solved_and_feasible(model_hyperbolic_1)
        println("Objective value for hyperbolic 1 model: ", objective_value(model_hyperbolic_1))
        beta_hyp1 = [value(beta00_H1), value(beta01_H1), value(beta02_H1)]
        c_hyp1 = [value(c00_H1)]
        SaveModel_oneInput("Power_hyperbolic_1", beta_hyp1, c_hyp1, objective_value(model_hyperbolic_1))
        # println("\nHyperbolic model 1:")
        # println("beta0 = ", value(beta00_H1))
        # println("beta1 = ", value(beta01_H1))
        # println("beta2 = ", value(beta02_H1))
        # println("c = ", value(c_H1))

    end


    if is_solved_and_feasible(model_hyperbolic_2)
        println("Objective value for hyperbolic 2 model: ", objective_value(model_hyperbolic_2))
        beta_hyp2 = [value(beta01_H2), value(beta02_H2)]
        c_hyp2 = [value(c00_H2)]
        
        SaveModel_oneInput("Power_hyperbolic_2", beta_hyp2, c_hyp2, objective_value(model_hyperbolic_2))
        # println("\nHyperbolic model 2:")
        # println("beta0 = ", value(beta00_H2))
        # println("beta1 = ", value(beta01_H2))
        # println("beta2 = ", value(beta02_H2))
        # println("c = ", value(c_H2))
    end


    if is_solved_and_feasible(model_hyperbolic_3)
        println("Objective value for hyperbolic 3 model: ", objective_value(model_hyperbolic_3))
        c_hyp3 = [value(c00_H3)]
        beta_hyp3 = [value(beta00_H3), value(beta01_H3), value(beta02_H3), value(beta03_H3)]
        SaveModel_oneInput("Power_hyperbolic_3", beta_hyp3, c_hyp3, objective_value(model_hyperbolic_3))
        # println("\nHyperbolic model 3:")
        # println("beta0 = ", value(beta00_H3))
        # println("beta1 = ", value(beta01_H3))
        # println("beta2 = ", value(beta02_H3))
        # println("beta3 = ", value(beta03_H3))
        # println("beta4 = ", value(beta04_H3))

    end

    if is_solved_and_feasible(model_hyperbolic_4)
        println("Objective value for hyperbolic 4 model: ", objective_value(model_hyperbolic_4))
        c_hyp4 = [value(c00_H4)]
        beta_hyp4 = [value(beta01_H4), value(beta02_H4), value(beta03_H4)]
        SaveModel_oneInput("Power_hyperbolic_4", beta_hyp4, c_hyp4, objective_value(model_hyperbolic_4))
        # println("\nHyperbolic model 4:")
        # println("beta1 = ", value(beta01_H4))
        # println("beta2 = ", value(beta02_H4))
        # println("beta3 = ", value(beta03_H4))
        # println("beta4 = ", value(beta04_H4))

    end

    if is_solved_and_feasible(model_logarithmic)
        println("Objective value for logarithmic model: ", objective_value(model_logarithmic))
        beta_log = [value(beta00_LG), value(beta01_LG), value(beta02_LG)]
        c_log = [value(c00_LG), value(c01_LG)]
        SaveModel_oneInput("Power_logarithmic", beta_log, c_log, objective_value(model_logarithmic))
        # println("\nLogarithmic model:")
        # println("beta0 = ", value(beta00_LG))
        # println("beta1 = ", value(beta01_LG))
        # println("beta2 = ", value(beta02_LG))

    end

    if is_solved_and_feasible(model_log_hyper1)
        println("Objective value for combination model: ", objective_value(model_log_hyper1))
        c_log_hyper1 = [value(c00_LH1), value(c01_LH1), value(c02_LH1)]
        beta_log_hyper1 = [value(beta00_LH1), value(beta01_LH1), value(beta02_LH1), value(beta03_LH1)]
        SaveModel_oneInput("Power_log_hyperbolic_1", beta_log_hyper1, c_log_hyper1, objective_value(model_log_hyper1))
        # println("\nCombination model:")
        # println("beta0 = ", value(beta00_LH1))
        # println("beta1 = ", value(beta01_LH1))
        # println("beta3 = ", value(beta03_LH1))
        # println("beta4 = ", value(beta04_LH1))
        # println("c = ", c_log_hyper1)
    end

    if is_solved_and_feasible(model_log_hyper2)
        println("Objective value for combination model: ", objective_value(model_log_hyper2))
        c_log_hyper2 = [value(c00_LH2), value(c01_LH2)]
        beta_log_hyper2 = [value(beta00_LH2), value(beta01_LH2), value(beta02_LH2)]
        SaveModel_oneInput("Power_log_hyperbolic_2", beta_log_hyper2, c_log_hyper2, objective_value(model_log_hyper2))
    end

    if is_solved_and_feasible(model_log_hyper3)
        println("Objective value for combination model: ", objective_value(model_log_hyper3))
        c_log_hyper3 = [value(c00_LH3), value(c01_LH3), value(c02_LH3)]
        beta_log_hyper3 = [value(beta00_LH3), value(beta01_LH3), value(beta02_LH3), value(beta03_LH3), value(beta04_LH3)]
        SaveModel_oneInput("Power_log_hyperbolic_3", beta_log_hyper3, c_log_hyper3, objective_value(model_log_hyper3))
    end


    println(Pw_scaled)


end

using Plots
using CSV
using DataFrames

function Plot_surrogateModel()

    # ============================================
    # Training and validation data
    # ============================================

    data_training = CSV.read(
        "/home/annkar/NewCode/Data/surrogate_data.csv",
        DataFrame
    )

    data_validation = CSV.read(
        "/home/annkar/NewCode/Data/surrogate_data_Pw_validation.csv",
        DataFrame
    )





    Pw_training = scaled_minmax(data_training.Power, minimum(data_training.Power), maximum(data_training.Power))
    HTO_training = data_training.HTO

    Pw_validation = scaled_minmax(data_validation.Power, minimum(data_training.Power), maximum(data_training.Power))
    HTO_validation = data_validation.HTO

    Pw = range(0.0, 1.0, length=100)

    

    # ============================================
    # Load models from dictionary
    # ============================================
    # model_linear = LoadModel_OneInput("Power_linear")

    model_quad = LoadModel_OneInput("Power_quadratic")

    model_hyp1 = LoadModel_OneInput("Power_hyperbolic_1")
    model_hyp2 = LoadModel_OneInput("Power_hyperbolic_2")
    model_hyp3 = LoadModel_OneInput("Power_hyperbolic_3")
    # model_hyp4 = LoadModel_OneInput("Power_hyperbolic_4")

    model_log = LoadModel_OneInput("Power_logarithmic")

    # model_comb1 = LoadModel_OneInput("Power_log_hyperbolic_1")
    model_comb2 = LoadModel_OneInput("Power_log_hyperbolic_2")
    # model_comb3 = LoadModel_OneInput("Power_log_hyperbolic_3")

    # ============================================
    # Extract beta, c and SSE
    # ============================================
    # beta_lin = model_linear["beta"]
    beta_quad = model_quad["beta"]

    beta_hyp1 = model_hyp1["beta"]
    beta_hyp2 = model_hyp2["beta"]
    beta_hyp3 = model_hyp3["beta"]
    # beta_hyp4 = model_hyp4["beta"]

    c_hyp1 = model_hyp1["c"]
    c_hyp2 = model_hyp2["c"]
    c_hyp3 = model_hyp3["c"]
    # c_hyp4 = model_hyp4["c"]

    beta_log = model_log["beta"]

    c_log = model_log["c"]
    

    # beta_comb1 = model_comb1["beta"]
    beta_comb2 = model_comb2["beta"]
    # beta_comb3 = model_comb3["beta"]

    # c_comb1 = model_comb1["c"]
    c_comb2 = model_comb2["c"]
    # c_comb3 = model_comb3["c"]

    # ============================================
    # Define regression models
    # ============================================

    # Linear
    HTO_lin(Pw) = 
        beta_lin[1] +
        beta_lin[2]*Pw 

    # Quadratic
    HTO_quad(Pw) =
        beta_quad[1] +
        beta_quad[2]*Pw +
        beta_quad[3]*Pw^2

    # Hyperbolic 1
    HTO_hyp1(Pw) =
        beta_hyp1[1] +
        beta_hyp1[2]*Pw +
        beta_hyp1[3]/(Pw+c_hyp1[1])

    # Hyperbolic 2
    HTO_hyp2(Pw) =
        beta_hyp2[1]*Pw +
        beta_hyp2[2]/(Pw + c_hyp2[1])

    # Hyperbolic 3
    HTO_hyp3(Pw) =
        beta_hyp3[1] +
        beta_hyp3[2]*Pw +
        beta_hyp3[3]/(beta_hyp3[4]*Pw + c_hyp3[1])

    # Hyperbolic 4
    HTO_hyp4(Pw) =
        beta_hyp4[1]*Pw +
        beta_hyp4[2]/(beta_hyp4[3]*Pw + c_hyp4[1])

    # Logarithmic
    HTO_log(Pw) =
        beta_log[1] +
        beta_log[2]*log(Pw + c_log[1]) +
        beta_log[3]*log(Pw + c_log[2])^2

    # Combination 1
    HTO_comb1(Pw) =
        beta_comb1[1] +
        beta_comb1[2]*log(Pw + c_comb1[1]) +
        beta_comb1[3]*log(Pw + c_comb1[2])^2 +
        beta_comb1[4]/(Pw + c_comb1[3])

    # Combination 2
    HTO_comb2(Pw) =
        beta_comb2[1] +
        beta_comb2[2]*log(Pw + c_comb2[1]) +
        beta_comb2[3]/(Pw + c_comb2[2])

    # Combination 3
    HTO_comb3(Pw) =  beta_comb3[1] + beta_comb3[2]*log(Pw + c_comb3[1]) +
        beta_comb3[3]*log(Pw + c_comb3[2])^2 +
        beta_comb3[4]/(beta_comb3[5]*Pw + c_comb3[3])


    # ============================================
    # Plot
    # ============================================

    p = plot( Pw_training, HTO_training, seriestype = :scatter,  label = "Training Data",
        xlabel = "Power (W)",
        ylabel = "HTO",
        title = "Surrogate Models",
        legend = :topright
    )

    scatter!(Pw_validation, HTO_validation, seriertype = :scatter,  label = "Validation Data"
    )

    # plot!(Pw, HTO_lin.(Pw),  label = "Linear, SSE = $(round(model_linear["SSE"], digits=10))")
    # plot!(Pw, HTO_quad.(Pw),  label = "Quadratic, SSE = $(round(model_quad["SSE"], digits=10))")
    plot!(Pw, HTO_log.(Pw),  label = "Logarithmic, SSE = $(round(model_log["SSE"], digits=5))")
    
    plot!(Pw, HTO_hyp1.(Pw),  label = "Hyperbolic 1, SSE = $(round(model_hyp1["SSE"], digits=6))")
    plot!(Pw, HTO_hyp2.(Pw),  label = "Hyperbolic 2, SSE = $(round(model_hyp2["SSE"], digits=4))")
    plot!(Pw, HTO_hyp3.(Pw),  label = "Hyperbolic 3, SSE = $(round(model_hyp3["SSE"], digits=6))")
    
    # plot!(Pw, HTO_comb1.(Pw), label = "Combination 1, SSE = $(round(model_comb1["SSE"], digits=10))")
    plot!(Pw, HTO_comb2.(Pw), label = "Combination 2, SSE = $(round(model_comb2["SSE"], digits=8))")
    # plot!(Pw, HTO_comb3.(Pw), label = "Combination 3, SSE = $(round(model_comb3["SSE"], digits=10))")

    display(p)

    return p
end

# using Plots
# function Plot_surrogateModel()

#     #points from the test data

#     Pw_test = data.Power
#     HTO_test = data.HTO

#     Pw_validation = data_validation.Power
#     HTO_validation = data_validation.HTO
#     Pw = range(6000.0, 60000.0, length=100)

#     #hyperbolic regression model
#     HTO_Hyper3(Pw) = 3.66 -5.4782619755281916e-5*Pw + -4.207060459518844e8/(-0.15603331492245331*Pw +1.532258273952951e9 )
#     HTO_hyper2(Pw) = 0.33008542216723874 -1.887900957654336e-6*Pw +  28435.177319662478/(Pw + 691.6690671617287 )
#     HTO_hyper1(Pw)= 0.5346020469710464 + -4.523651454068847e-6*Pw + 24433.01113165403/(Pw)
#     HTO_hyper4(Pw) = 2.4151793460947953e-6*Pw + 18.430160021445616/(0.0005265071104831629*Pw +0.8956009805889193 )
    
#     #other regression models
#     HTO_quadratic(Pw) = 5.008270813508028 -0.0001899823712054724*Pw + 2.0484806766357358e-9*Pw^2
#     HTO_logarithmic(Pw) = 17.56897725627651 -1.566209537640515*log(Pw)
#     HTO_logaritmic2(Pw) = 81.77128521397523 -14.606922777139543*log(Pw) + 0.6583110050070654*log(Pw)^2
    
#     HTO_combination_1(Pw) = 128.59044007335902 -20.51525326050656*log(Pw)+ 0.8263621013181012*log(Pw)^2 -165258.5333034125/(Pw + 14427.659893591028)
#     HTO_combination_2(Pw) = 68.11086908898574 -5.403568223608232*log(Pw) + -835163.8565489864/(Pw + 44492.47663937794)
#     HTO_combination_3(Pw) = -2.1567897861465037 + 0.5381549477976186*log(Pw) -0.029258512740254403*log(Pw)^2 + 1538.1436805093158/( 0.054725548267628134*Pw + 32.91486137054424)

#     # plot(Pw, HTO_first.(Pw), label="HTO = b0 + b1Pw + b2/(b3Pw + b4)", xlabel="Power (W)", ylabel="HTO (mol/s)", title="Surrogate Models", legend=:topright)
#     plot(Pw_test, HTO_test, seriestype=:scatter, label="Test Data", xlabel="Power (W)",title="Surrogate Models", ylabel="HTO (mol/s)", legend=:topright)
#     plot!(Pw_validation, HTO_validation, seriestype=:scatter, label="Validation Data", xlabel="Power (W)",title="Surrogate Models", ylabel="HTO (mol/s)", legend=:topright)
#     # plot!(Pw, HTO_Hyper3.(Pw), label="HTO = b0 + b1Pw + b2/(b3Pw + b4), SSE = 3.97", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
#         #Hyperbolic
#     plot!(Pw, HTO_hyper2.(Pw), label="HTO_hyp2 = b0 + b1Pw + b2/(Pw + b3), SSE = 3.7e-6", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
#     plot!(Pw, HTO_hyper1.(Pw), label="HTO_hyp1 = b0 + b1Pw + b2/(Pw), SSE = 0.0019 ", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
#     plot!(Pw, HTO_hyper4.(Pw), label="HTO_hyp4 = b0 + b1Pw + b2/(b3Pw + b4), SSE = 0.0037", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
#     # plot!(Pw, HTO_quadratic.(Pw), label="HTO = b0 + b1Pw + b2Pw^2, SSE = 1.1", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
#         #logarithmic
#     # plot!(Pw, HTO_logarithmic.(Pw), label="HTO_log = b0 + b1*log(Pw), SSE = 1.0238", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
#     # plot!(Pw, HTO_logaritmic2.(Pw), label="HTO_log2 = b0 + b1*log(Pw) + b2*log(Pw)^2, SSE = 0.026", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
#         #combination
#     # plot!(Pw, HTO_combination_1.(Pw), label="HTO_com_1 = b0 + b1*log(Pw) + b2*log(Pw)^2 + b3/(Pw + b4), SSE = 7.52e-6", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
#     # plot!(Pw, HTO_combination_2.(Pw), label="HTO_com_2 = b0 + b1*log(Pw) + b2/(Pw + b3), SSE = 0.00223", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
#     # plot!(Pw, HTO_combination_3.(Pw), label="HTO_com_3 = b0 + b1*log(Pw) + b2*log(Pw)^2 + b3/(b4*Pw + b5), SSE = 8.87e-9", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)

# end

# display(Plot_surrogateModel())

function SurrogteModel_validationData()
    HTO_Hyper3(Pw) = 3.66 -5.4782619755281916e-5*Pw + -4.207060459518844e8/(-0.15603331492245331*Pw +1.532258273952951e9 )
    HTO_hyper2(Pw) = 0.33008542216723874 -1.887900957654336e-6*Pw +  28435.177319662478/(Pw + 691.6690671617287 )
    HTO_hyper1(Pw)= 0.5346020469710464 + -4.523651454068847e-6*Pw + 24433.01113165403/(Pw)
    HTO_hyper4(Pw) = 2.4151793460947953e-6*Pw + 18.430160021445616/(0.0005265071104831629*Pw +0.8956009805889193 )
    
    #other regression models
    HTO_quadratic(Pw) = 5.008270813508028 -0.0001899823712054724*Pw + 2.0484806766357358e-9*Pw^2
    HTO_logarithmic(Pw) = 17.56897725627651 -1.566209537640515*log(Pw)
    HTO_logaritmic2(Pw) = 81.77128521397523 -14.606922777139543*log(Pw) + 0.6583110050070654*log(Pw)^2
    
    HTO_combination_1(Pw) = 128.59044007335902 -20.51525326050656*log(Pw)+ 0.8263621013181012*log(Pw)^2 -165258.5333034125/(Pw + 14427.659893591028)
    HTO_combination_2(Pw) = 68.11086908898574 -5.403568223608232*log(Pw) + -835163.8565489864/(Pw + 44492.47663937794)
    HTO_combination_3(Pw) = -2.1567897861465037 + 0.5381549477976186*log(Pw) -0.029258512740254403*log(Pw)^2 + 1538.1436805093158/( 0.054725548267628134*Pw + 32.91486137054424)

    functions = [HTO_Hyper3, HTO_hyper2, HTO_hyper1, HTO_hyper4, HTO_quadratic, HTO_logarithmic, HTO_logaritmic2, HTO_combination_1, HTO_combination_2, HTO_combination_3]

    for i in 1:length(functions)
        f = functions[i]
        predictions = f.(data_validation.Power)
        residuals = predictions - data_validation.HTO
        SSE = sum(residuals.^2)
        RSME = sqrt(SSE / length(data_validation.Power))
        println("SSE for function $f: $SSE")
        println("RMSE for function $f: $RSME")
        println("--------------------------------------------------")
    end


end

function SurrogateModel_SeveralInputs()
    data = CSV.read(joinpath("/home/annkar/NewCode/Data", "surrogate_data_Pw_T_training.csv"), DataFrame)

    
    Pw_train = data.Power
    T_train = data.Temperature
    HTO_train = data.HTO

    Pw_train_scaled = scaled_minmax(Pw_train, minimum(Pw_train), maximum(Pw_train))
    T_train_scaled = scaled_minmax(T_train, minimum(T_train), maximum(T_train))
    HTO_train_scaled = scaled_minmax(HTO_train, minimum(HTO_train), maximum(HTO_train))

    N = length(Pw_train_scaled)

    ##=====================================================================
    #defining different optimizer models
    ##=====================================================================
    model_hyperbolic_Pw_T_01 = Model(Ipopt.Optimizer)
    set_attribute(model_hyperbolic_Pw_T_01, "linear_solver", "mumps")

    model_quadratic_Pw_T_01 = Model(Ipopt.Optimizer)
    set_attribute(model_quadratic_Pw_T_01, "linear_solver", "mumps")

    model_logarithmic_Pw_T_01 = Model(Ipopt.Optimizer)
    set_attribute(model_logarithmic_Pw_T_01, "linear_solver", "mumps")

    ##=====================================================================
    #defining variables for the differet polynomial regression models
    ##=====================================================================

    @variable(model_hyperbolic_Pw_T_01, beta00_1)
    @variable(model_hyperbolic_Pw_T_01, beta01_1)
    @variable(model_hyperbolic_Pw_T_01, beta02_1)
    @variable(model_hyperbolic_Pw_T_01, beta03_1)
    @variable(model_hyperbolic_Pw_T_01, beta04_1)
    @variable(model_hyperbolic_Pw_T_01, beta05_1)
    @variable(model_hyperbolic_Pw_T_01, beta06_1)
    @variable(model_hyperbolic_Pw_T_01, c0 >= 0.0)
    @variable(model_hyperbolic_Pw_T_01, c1 >= 0.0)
    @variable(model_hyperbolic_Pw_T_01, c2 >= 0.0)

    @variable(model_hyperbolic_Pw_T_01, HTO_01_1[1:N])
    #create the constraints for the polynomial regression model
    # @constraint(model_hyperbolic_Pw_T_01, [i=1:N], beta00_1 + beta01_1*1/(Pw_train_scaled[i]+c0)
    #  + beta02_1*(1/(T_train_scaled[i]+c1))
    #  + beta03_1*Pw_train_scaled[i]*T_train_scaled[i] 
    #  + beta04_1*Pw_train_scaled[i]^2*T_train_scaled[i]
    #  + beta05_1*T_train_scaled[i]/(Pw_train_scaled[i]+c2)  == HTO_01_1[i])

    @constraint(model_hyperbolic_Pw_T_01, [i=1:N], beta00_1 + beta01_1*1/(Pw_train_scaled[i]+c0)
     + beta02_1*T_train_scaled[i]
     + beta03_1*T_train_scaled[i]^2
     + beta04_1*Pw_train_scaled[i]*T_train_scaled[i] 
     + beta05_1*Pw_train_scaled[i]^2*T_train_scaled[i]
     + beta06_1*T_train_scaled[i]/(Pw_train_scaled[i]+c2)  == HTO_01_1[i])
    
    
    # defining the objective function for the polynomial regression model
    @objective(model_hyperbolic_Pw_T_01, Min, sum((HTO_01_1[i] - HTO_train_scaled[i])^2 for i in 1:length(Pw_train_scaled)))

    ##=====================================================================
    #Quadratic regression model
    ##====================================================================

        #new model with different scaling for the inputs
    @variable(model_quadratic_Pw_T_01, beta00_2)
    @variable(model_quadratic_Pw_T_01, beta01_2)
    @variable(model_quadratic_Pw_T_01, beta02_2)
    @variable(model_quadratic_Pw_T_01, beta03_2)
    @variable(model_quadratic_Pw_T_01, beta04_2)
    @variable(model_quadratic_Pw_T_01, HTO_01_2[1:N])
    #create the constraints for the polynomial regression model
    @constraint(model_quadratic_Pw_T_01, [i=1:N], beta00_2 + beta01_2*Pw_train_scaled[i]^2 + beta02_2*T_train_scaled[i]^2 + 
    beta03_2*Pw_train_scaled[i]*T_train_scaled[i] + beta04_2*Pw_train_scaled[i]^2*T_train_scaled[i] == HTO_01_2[i])

    # defining the objective function for the polynomial regression model
    @objective(model_quadratic_Pw_T_01, Min, sum((HTO_01_2[i] - HTO_train_scaled[i])^2 for i in 1:N))


    ##=====================================================================
    @variable(model_logarithmic_Pw_T_01, beta00_3)
    @variable(model_logarithmic_Pw_T_01, beta01_3)
    @variable(model_logarithmic_Pw_T_01, beta02_3)
    @variable(model_logarithmic_Pw_T_01, beta03_3)
    @variable(model_logarithmic_Pw_T_01, beta04_3)
    @variable(model_logarithmic_Pw_T_01, beta05_3)

    @variable(model_logarithmic_Pw_T_01, c >= 0.0)
    @variable(model_logarithmic_Pw_T_01, HTO_01_3[1:N])
    #create the constraints for the polynomial regression model
    @constraint(model_logarithmic_Pw_T_01, [i=1:N], beta00_3 + beta01_3*log(Pw_train_scaled[i]+c) + beta02_3*log(Pw_train_scaled[i]+c)^2
                                                     + beta03_3*log(T_train_scaled[i]+c) + beta04_3*log(T_train_scaled[i]+c)^2
                                                     + beta05_3*log(Pw_train_scaled[i]+c)*log(T_train_scaled[i]+c)== HTO_01_3[i])
    # defining the objective function for the polynomial regression model
    @objective(model_logarithmic_Pw_T_01, Min, sum((HTO_01_3[i] - HTO_train_scaled[i])^2 for i in 1:N))

    ##=====================================================================
    #solving the optimization probelm
    ##=====================================================================

    println("Optimizing the hyperbolic model...")
    optimize!(model_hyperbolic_Pw_T_01)

    println("Optimizing the model...")
    optimize!(model_quadratic_Pw_T_01)

    println("Optimizing the model...")
    optimize!(model_logarithmic_Pw_T_01)


    ##=====================================================================
    #Cheking the solver status and printing the results
    ##=====================================================================

    println("Solver status for combination model 1: ", termination_status(model_hyperbolic_Pw_T_01))
    println("Solver status for combination model 2: ", termination_status(model_quadratic_Pw_T_01))
    println("Solver status for combination model 3: ", termination_status(model_logarithmic_Pw_T_01))





    if is_solved_and_feasible(model_hyperbolic_Pw_T_01)
        println("Objective value for hyperbolic 1 model: ", objective_value(model_hyperbolic_Pw_T_01))
        beta_hyp = [value(beta00_1),
            value(beta01_1),
            value(beta02_1),
            value(beta03_1),
            value(beta04_1),
            value(beta05_1),
            value(beta06_1)
        ]
        c_hyp = [
            value(c0),
            value(c1),
            value(c2)
        ]

        SSE_hyp = [objective_value(model_hyperbolic_Pw_T_01)]

        println(beta_hyp)
        println(c_hyp)

        SaveModel_TwoInput("Pw_T_hyperbolic_v1", beta_hyp, c_hyp, SSE_hyp)        

        println("\nHyperbolic model 1:")
        println("beta0 = ", value(beta00_1))
        println("beta1 = ", value(beta01_1))
        println("beta2 = ", value(beta02_1))
        println("beta3 = ", value(beta03_1))
        println("beta4 = ", value(beta04_1))
        println("beta5 = ", value(beta05_1))
        println("beta6 = ", value(beta06_1))
        println("c = ", value(c0))
        println("c = ", value(c1))
        println("c = ", value(c2))

    end

    if is_solved_and_feasible(model_quadratic_Pw_T_01)
        println("Objective value for quadratic model: ", objective_value(model_quadratic_Pw_T_01))
        beta_quad = [
        value(beta00_2),
        value(beta01_2),
        value(beta02_2),
        value(beta03_2),
        value(beta04_2)
        ]
        SSE_quad = [objective_value(model_quadratic_Pw_T_01)]

        SaveModel_TwoInput("Pw_T_quadratic_v1", beta_quad, Float64[], SSE_quad)
        println("\nQuadratic model 2:")
        println("beta0 = ", value(beta00_2))
        println("beta1 = ", value(beta01_2))
        println("beta2 = ", value(beta02_2))
        println("beta3 = ", value(beta03_2))
        println("beta4 = ", value(beta04_2))

    end
    
    if is_solved_and_feasible(model_logarithmic_Pw_T_01)
        println("Objective value for logarithmic model: ", objective_value(model_logarithmic_Pw_T_01))
        beta_log = [
             value(beta00_3),
             value(beta01_3),
             value(beta02_3),
             value(beta03_3),
             value(beta04_3),
             value(beta05_3)
         ]

         c_log = [value(c)]

         SSE_log = [objective_value(model_logarithmic_Pw_T_01)]

         SaveModel_TwoInput("Pw_T_logarithmic_v1", beta_log, c_log, SSE_log)


        println("\nLogarithmic model:")
        println("beta0 = ", value(beta00_3))
        println("beta1 = ", value(beta01_3))
        println("beta2 = ", value(beta02_3))
        println("beta3 = ", value(beta03_3))
        println("beta4 = ", value(beta04_3))
        println("beta5 = ", value(beta05_3))
        println("c = ", value(c))

    end
end

# using Plots

# data0 = CSV.read(joinpath("/home/annkar/NewCode/Data", "Plot_HTO_vs_varyT.csv"), DataFrame)

# HTO_9 = data0.HTO[1:217]
# HTO_30 = data0.HTO[218:434]
# HTO_60 = data0.HTO[435:651]

# Pc_9 = data0.Power[1:217]*1e-3
# Pc_30 = data0.Power[218:434]*1e-3
# Pc_60 = data0.Power[435:651]*1e-3

# plot(Pc_9, HTO_9, label="T=343.15 K, Pc = 30 bar")
# plot!(Pc_30, HTO_30, label="T=313.15 K, Pc = 30 bar")
# plot!(Pc_60, HTO_60, label="T=293.15 K, Pc = 30 bar", xlabel="Power (kW)", ylabel="HTO (%)", title="HTO vs Power, const Pressure", legend=:topright)



using PlotlyJS, CSV, DataFrames
function Plot_surrogateModel_SeveralInputs()
    data_train = CSV.read(joinpath("/home/annkar/NewCode/Data", "surrogate_data_Pw_T_training.csv"), DataFrame)

    Pw_train = data_train.Power
    T_train = data_train.Temperature
    HTO_train = data_train.HTO

    Pw_train_scaled = scaled_minmax(Pw_train, minimum(Pw_train), maximum(Pw_train))
    T_train_scaled = scaled_minmax(T_train, minimum(T_train), maximum(T_train))
    HTO_train_scaled = scaled_minmax(HTO_train, minimum(HTO_train), maximum(HTO_train))

 # Create 3D scatter plot

    plotScatter = false

    if plotScatter
        trace = PlotlyJS.scatter3d(x = Pw_train, y = T_train, z = HTO_train, mode = "markers",
            marker = attr(size = 5, color = HTO_train, colorscale = "Viridis",  showscale = true, colorbar = attr(title = "HTO")
            ),
            name = "Training data"
        )

        trace2 = PlotlyJS.scatter3d(x = Pw_train_scaled, y = T_train_scaled, z = HTO_train_scaled, mode = "markers",
            marker = attr(size = 5, color = HTO_train_scaled, colorscale = "Viridis",  showscale = true, colorbar = attr(title = "HTO")
            ),
            name = "Training scaled data"
        )
        # Layout
        layout = PlotlyJS.Layout(title = "HTO as a function of Power and Temperature",
            scene = attr(
                xaxis = attr(title = "Power"),
                yaxis = attr(title = "Temperature"),
                zaxis = attr(title = "HTO")
            )
        )


        scatterplot = PlotlyJS.Plot(trace, layout)
        display(scatterplot)
    
    
    
    end

    model_hyp = LoadModel_TwoInput("Pw_T_hyperbolic_v1")
    beta_hyp = model_hyp["beta"]
    c_hyp = model_hyp["c"]

    model_quad = LoadModel_TwoInput("Pw_T_quadratic_v1")
    beta_quad = model_quad["beta"]
    c_quad = model_quad["c"]

    model_log = LoadModel_TwoInput("Pw_T_logarithmic_v1")
    beta_log = model_log["beta"]
    c_log = model_log["c"]



    println(beta_hyp)
    plotPolynomial = true

    if plotPolynomial
        Pw_grid = range(-0.3, 0.3, length=50)
        T_grid = range(-0.3, 0.3, length=50)

        Pw_grid_1 = range(0.0, 1.0, length=50)
        T_grid_1 = range(0.0, 1.0, length=50)


        # Define the polynomial regression model
                
        HTO_model_hyp(Pw,T)= beta_hyp[1] + beta_hyp[2]/(Pw + c_hyp[1]) + beta_hyp[3]/(T + c_hyp[2]) + beta_hyp[4]*Pw*T + beta_hyp[5]*Pw^2*T + beta_hyp[6]*T/(Pw + c_hyp[3])
        HTO_model_quad(Pw, T) =  beta_quad[1] + beta_quad[2]*Pw^2 + beta_quad[3]*T^2 + beta_quad[4]*Pw*T + beta_quad[5]*Pw^2*T
        
        HTO_model_log(Pw, T) = beta_log[1] + beta_log[2]*log(Pw + c_log[1]) + beta_log[3]*log(Pw + c_log[1])^2 + beta_log[4]*log(T + c_log[1]) + beta_log[5]*log(T + c_log[1])^2 + beta_log[6]*log(Pw + c_log[1])*log(T + c_log[1])
        #create a meshgrid for HTO

        
        HTO_grid_hyp = [HTO_model_hyp(Pw, T) for Pw in Pw_grid_1, T in T_grid_1]

        HTO_grid_quad = [HTO_model_quad(Pw, T) for Pw in Pw_grid_1, T in T_grid_1]

        HTO_grid_log = [HTO_model_log(Pw, T) for Pw in Pw_grid_1, T in T_grid_1]
        
        # Create 3D surface plots
        surface_trace01 = PlotlyJS.surface(x = Pw_grid_1, y = T_grid_1, z = HTO_grid_log,
            colorscale = "Viridis", showscale = true, colorbar = attr(title = "HTO"),
            name = "Surrogate Model"
        )     
        
        # Layout
        trace_train = PlotlyJS.scatter3d(x = Pw_train_scaled, y = T_train_scaled, z = HTO_train_scaled, mode = "markers",
            marker = attr(size = 5, color = "red"
            ),
            name = "Training scaled data"
        )        
        
        layout_surface = PlotlyJS.Layout(title = "Surrogate Model: HTO as a function of Power and Temperature",
            scene = attr(
                xaxis = attr(title = "Power"),
                yaxis = attr(title = "Temperature"),
                zaxis = attr(title = "HTO")
            )
        )

        surface_plot = PlotlyJS.Plot([surface_trace01, trace_train],  layout_surface)
        display(surface_plot)
    end


end

export fit_surrogate, Plot_surrogateModel, SurrogateModel_validationData, Plot_surrogateModel_SeveralInputs, SurrogateModel_SeveralInputs


