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

data = CSV.read(joinpath(@__DIR__, "surrogate_data.csv"), DataFrame)

function fit_surrogate( ;csv_path = joinpath(@__DIR__, "surrogate_data.csv"))
    # gets the data from the CSV file
    data = CSV.read(csv_path, DataFrame)

    ##=====================================================================
    #defining different optimizer models
    ##=====================================================================

        #linear

    model_linear = Model(Ipopt.Optimizer)
    set_attribute(model_linear, "linear_solver", "mumps")

        #quadratic

    model_quadratic = Model(Ipopt.Optimizer)
    set_attribute(model_quadratic, "linear_solver", "mumps")

        #concave

    model_concave = Model(Ipopt.Optimizer)
    set_attribute(model_concave, "linear_solver", "mumps")

    ##=====================================================================
    # extracts the values values from csv data
    ##=====================================================================
    Pw = data.Power
    HTO_exp = data.HTO
    N = length(Pw)

    Pw_center = (60000.0+6000.0)/2 #funker kun på lineære, finne 
    Pw_range = (60000.0-6000.0)/2


    Pw_scaled = Pw 
    ##=====================================================================
    #defining variables for the differet polynomial regression models
    ##=====================================================================
        #linear regression model
    @variable(model_linear, beta0_L )
    @variable(model_linear, beta1_L)
    @variable(model_linear, HTO_L[1:N])
    
    #create the constraints for the polynomial regression model
    @constraint(model_linear, [i=1:N], beta0_L + beta1_L*Pw_scaled[i] == HTO_L[i])

    # defining the objective function for the polynomial regression model
    @objective(model_linear, Min, sum((HTO_L[i] - HTO_exp[i])^2 for i in 1:N))
    
        #quadratic regression model
    @variable(model_quadratic, beta0_Q )
    @variable(model_quadratic, beta1_Q )
    @variable(model_quadratic, beta2_Q )
    @variable(model_quadratic, HTO_Q[1:N]) 
    #create the constraints for the polynomial regression model
    @constraint(model_quadratic, [i=1:N], beta0_Q + beta1_Q*Pw_scaled[i] + beta2_Q*Pw_scaled[i]^2 == HTO_Q[i])
    # defining the objective function for the polynomial regression model
    @objective(model_quadratic, Min, sum((HTO_Q[i] - HTO_exp[i])^2 for i in 1:N))

        #concave regression model
    @variable(model_concave, beta01_C )
    @variable(model_concave, beta02_C )
    @variable(model_concave, beta03_C, start=1.0)
    @variable(model_concave, beta04_C, start=1.0)
    @variable(model_concave, HTO_C[1:N]) 
    #create the constraints for the polynomial regression model
    # @constraint(model_concave, [i=1:N], beta00_C + beta01_C*Pw_scaled[i] + beta02_C*Pw_scaled[i] 
                # + beta03_C*(Pw_center/Pw[i]) == HTO_C[i])
    @constraint(model_concave, [i=1:N], beta01_C*Pw_scaled[i] + beta02_C*(1/(beta03_C*Pw_scaled[i]+beta04_C)) == HTO_C[i])
    @constraint(model_concave, [i=1:N], beta01_C*Pw_scaled[i] + beta02_C*(1/(beta03_C*Pw_scaled[i]+beta04_C)) == HTO_C[i])
    # defining the objective function for the polynomial regression model
    @objective(model_concave, Min, sum((HTO_C[i] - HTO_exp[i])^2 for i in 1:N))



    ##=====================================================================
    #solving the optimization probelm
    ##=====================================================================
    
    println("Optimizing the linear model...")   
    optimize!(model_linear)
    
    println("Optimizing the quadratic model...")
    optimize!(model_quadratic)

    println("Optimizing the concave model...")
    optimize!(model_concave)




    println("Solver status for linear model: ", termination_status(model_linear))
    println("Solver status for quadratic model: ", termination_status(model_quadratic))
    println("Solver status for concave model: ", termination_status(model_concave))

    if is_solved_and_feasible(model_linear)
        println("Objective value for linear model: ", objective_value(model_linear))
        println("\nLinear model:")
        println("beta0 = ", value(beta0_L))
        println("beta1 = ", value(beta1_L))
    end

    if is_solved_and_feasible(model_quadratic)
        println("Objective value for quadratic model: ", objective_value(model_quadratic))
        println("\nQuadratic model:")
        println("beta0 = ", value(beta0_Q))
        println("beta1 = ", value(beta1_Q))
        println("beta2 = ", value(beta2_Q))

    end

    if is_solved_and_feasible(model_concave)
        println("Objective value for inverse model: ", objective_value(model_concave))
        println("\nInverse model:")
        # println("beta0 = ", value(beta00_C))
        println("beta1 = ", value(beta01_C))
        println("beta2 = ", value(beta02_C))
        println("beta3 = ", value(beta03_C))
        println("beta4 = ", value(beta04_C))

    end

    println(Pw_scaled)


end

using Plots
function Plot_surrogateModel()

    #points from the test data

    Pw_test = data.Power
    HTO_test = data.HTO

    Pw = range(6000.0, 60000.0, length=100)

    HTO_first(Pw) = 3.66 -5.4782619755281916e-5*Pw + -4.207060459518844e8/(-0.15603331492245331*Pw +1.532258273952951e9 )
    HTO_second(Pw) = 2.4151793460947953e-6*Pw + 18.430160021445616/(0.0005265071104831629*Pw +0.8956009805889193 )
    HTO_third(Pw) = 5.008270813508028 -  -0.0001899823712054724*Pw + 2.0484806766357358e-9*Pw^2

    # plot(Pw, HTO_first.(Pw), label="HTO = b0 + b1Pw + b2/(b3Pw + b4)", xlabel="Power (W)", ylabel="HTO (mol/s)", title="Surrogate Models", legend=:topright)
    plot(Pw_test, HTO_test, seriestype=:scatter, label="Test Data", xlabel="Power (W)",title="Surrogate Models", ylabel="HTO (mol/s)", legend=:topright)
    plot!(Pw, HTO_second.(Pw), label="HTO = b0Pw + b1/(b2Pw + b3), SSE = 0.0037", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
    plot!(Pw, HTO_third.(Pw), label="HTO = b0 + b1Pw + b2Pw^2, SSE = 1.1", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)


end

display(Plot_surrogateModel())

export fit_surrogate, Plot_surrogateModel