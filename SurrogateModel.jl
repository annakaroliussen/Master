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

function scaled(x)
    return (x .- ((maximum(x)+minimum(x))/2))./(maximum(x)-minimum(x))/2
end

data = CSV.read(joinpath(@__DIR__, "surrogate_data.csv"), DataFrame)

data_validation = CSV.read(joinpath("/home/annkar/NewCode/Data", "surrogate_data_Pw_validation.csv"), DataFrame)

function fit_surrogate(;csv_path = joinpath(@__DIR__, "surrogate_data.csv"))
    # gets the data from the CSV file
    data = CSV.read(csv_path, DataFrame)

    ##=====================================================================
    #defining different optimizer models
    ##=====================================================================
        #quadratic

    model_quadratic = Model(Ipopt.Optimizer)
    set_attribute(model_quadratic, "linear_solver", "mumps")

        #hyperbolic
    model_hyperbolic_1 = Model(Ipopt.Optimizer)
    set_attribute(model_hyperbolic_1, "linear_solver", "mumps")
    
    model_hyperbolic_2 = Model(Ipopt.Optimizer)
    set_attribute(model_hyperbolic_2, "linear_solver", "mumps")

    model_hyperbolic_3 = Model(Ipopt.Optimizer)
    set_attribute(model_hyperbolic_3, "linear_solver", "mumps")

    model_hyperbolic_4 = Model(Ipopt.Optimizer)
    set_attribute(model_hyperbolic_4, "linear_solver", "mumps")
    
        #logarithmic
    model_logarithmic = Model(Ipopt.Optimizer)
    set_attribute(model_logarithmic, "linear_solver", "mumps")

        #combination of logarithmic and hyperbolic
    model_log_hyper = Model(Ipopt.Optimizer)
    set_attribute(model_log_hyper, "linear_solver", "mumps")


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
    
        #quadratic regression model
    @variable(model_quadratic, beta0_Q )
    @variable(model_quadratic, beta1_Q )
    @variable(model_quadratic, beta2_Q )
    @variable(model_quadratic, HTO_Q[1:N]) 
    #create the constraints for the polynomial regression model
    @constraint(model_quadratic, [i=1:N], beta0_Q + beta1_Q*Pw_scaled[i] + beta2_Q*Pw_scaled[i]^2 == HTO_Q[i])
    # defining the objective function for the polynomial regression model
    @objective(model_quadratic, Min, sum((HTO_Q[i] - HTO_exp[i])^2 for i in 1:N))

        #hyperbolic regression model
    
        #First version of hyperbolic regression model
    @variable(model_hyperbolic_1, beta00_H1)
    @variable(model_hyperbolic_1, beta01_H1 )
    @variable(model_hyperbolic_1, beta02_H1 )
    @variable(model_hyperbolic_1, HTO_H1[1:N]) 
    #create the constraints for the polynomial regression model
    @constraint(model_hyperbolic_1, [i=1:N], beta00_H1 + beta01_H1*Pw_scaled[i] + beta02_H1*(1/(Pw_scaled[i])) == HTO_H1[i])
    # defining the objective function for the polynomial regression model
    @objective(model_hyperbolic_1, Min, sum((HTO_H1[i] - HTO_exp[i])^2 for i in 1:N))
    
    
        #Second version of hyperbolic regression model
    @variable(model_hyperbolic_2, beta00_H2)
    @variable(model_hyperbolic_2, beta01_H2 )
    @variable(model_hyperbolic_2, beta02_H2 )
    @variable(model_hyperbolic_2, beta03_H2, start=1.0)
    @variable(model_hyperbolic_2, HTO_H2[1:N]) 
    #create the constraints for the polynomial regression model
    @constraint(model_hyperbolic_2, [i=1:N], beta00_H2 + beta01_H2*Pw_scaled[i] + beta02_H2*(1/(Pw_scaled[i]+beta03_H2)) == HTO_H2[i])
    # defining the objective function for the polynomial regression model
    @objective(model_hyperbolic_2, Min, sum((HTO_H2[i] - HTO_exp[i])^2 for i in 1:N))



        #third version of hyperbolic regression model
    @variable(model_hyperbolic_3, beta00_H3)
    @variable(model_hyperbolic_3, beta01_H3 )
    @variable(model_hyperbolic_3, beta02_H3 )
    @variable(model_hyperbolic_3, beta03_H3, start=1.0)
    @variable(model_hyperbolic_3, beta04_H3, start=1.0)
    @variable(model_hyperbolic_3, HTO_H3[1:N]) 
    #create the constraints for the polynomial regression model
    @constraint(model_hyperbolic_3, [i=1:N], beta00_H3 + beta01_H3*Pw_scaled[i] + beta02_H3*(1/(beta03_H3*Pw_scaled[i]+beta04_H3)) == HTO_H3[i])
    # defining the objective function for the polynomial regression model
    @objective(model_hyperbolic_3, Min, sum((HTO_H3[i] - HTO_exp[i])^2 for i in 1:N))

        #fourth version of hyperbolic regression model
    @variable(model_hyperbolic_4, beta01_H4 )
    @variable(model_hyperbolic_4, beta02_H4 )
    @variable(model_hyperbolic_4, beta03_H4, start=1.0)
    @variable(model_hyperbolic_4, beta04_H4, start=1.0)
    @variable(model_hyperbolic_4, HTO_H4[1:N]) 
    #create the constraints for the polynomial regression model
    @constraint(model_hyperbolic_4, [i=1:N], beta01_H4*Pw_scaled[i] + beta02_H4*(1/(beta03_H4*Pw_scaled[i]+beta04_H4)) == HTO_H4[i])
    # defining the objective function for the polynomial regression model
    @objective(model_hyperbolic_4, Min, sum((HTO_H4[i] - HTO_exp[i])^2 for i in 1:N))

        #logarithmic regression model
    @variable(model_logarithmic, beta00_LG)
    @variable(model_logarithmic, beta01_LG)
    @variable(model_logarithmic, beta02_LG)
    @variable(model_logarithmic, HTO_LG[1:N])
    #create the constraints for the logarithmic regression model
    @constraint(model_logarithmic, [i=1:N], beta00_LG  + beta01_LG*log(Pw_scaled[i]) + beta02_LG*log(Pw_scaled[i])^2 == HTO_LG[i])
    # defining the objective function for the logarithmic regression model
    @objective(model_logarithmic, Min, sum((HTO_LG[i] - HTO_exp[i])^2 for i in 1:N))

        #combination of logarithmic and hyperbolic regression model

    @variable(model_log_hyper, beta00_LH)
    @variable(model_log_hyper, beta01_LH)
    @variable(model_log_hyper, beta02_LH)
    @variable(model_log_hyper, beta03_LH)
    @variable(model_log_hyper, beta04_LH)
    # @variable(model_log_hyper, beta05_LH, start=1.0)
    @variable(model_log_hyper, HTO_LH[1:N])
    #create the constraints for the combination model
    @constraint(model_log_hyper, [i=1:N], beta00_LH + beta01_LH*log(Pw_scaled[i]) + beta02_LH*log(Pw_scaled[i])^2 + beta03_LH*(1/(Pw_scaled[i]+beta04_LH)) == HTO_LH[i])
    # @constraint(model_log_hyper, [i=1:N], beta00_LH + beta01_LH*log(Pw_scaled[i]) + beta02_LH*log(Pw_scaled[i])^2 + beta03_LH*(1/(beta05_LH*Pw_scaled[i]+beta04_LH)) == HTO_LH[i])
    # defining the objective function for the combination model
    @objective(model_log_hyper, Min, sum((HTO_LH[i] - HTO_exp[i])^2 for i in 1:N))




    ##=====================================================================
    #solving the optimization probelm
    ##=====================================================================
    
   
    # println("Optimizing the logarithmic model...")
    # optimize!(model_logarithmic)

    # println("Optimizing the hyperbolic model...")
    # optimize!(model_hyperbolic_1)   

    # println("Optimizing the hyperbolic model...")
    # optimize!(model_hyperbolic_2)

    # println("Optimizing the hyperbolic model...")
    # optimize!(model_hyperbolic_4)

    println("Optimizing the combination model...")
    optimize!(model_log_hyper)






    # println("Solver status for linear model: ", termination_status(model_linear))
    # println("Solver status for quadratic model: ", termination_status(model_quadratic))
    # println("Solver status for hyperbolic model: ", termination_status(model_hyperbolic_1))
    # println("Solver status for hyperbolic model: ", termination_status(model_hyperbolic_2))
    # # println("Solver status for hyperbolic model: ", termination_status(model_hyperbolic_3))
    # println("Solver status for hyperbolic model: ", termination_status(model_hyperbolic_4))
    # println("Solver status for logarithmic model: ", termination_status(model_logarithmic))
    println("Solver status for combination model: ", termination_status(model_log_hyper))

    if is_solved_and_feasible(model_hyperbolic_1)
        println("Objective value for hyperbolic 1 model: ", objective_value(model_hyperbolic_1))
        println("\nHyperbolic model 1:")
        println("beta0 = ", value(beta00_H1))
        println("beta1 = ", value(beta01_H1))
        println("beta2 = ", value(beta02_H1))

    end

    if is_solved_and_feasible(model_hyperbolic_2)
        println("Objective value for hyperbolic 2 model: ", objective_value(model_hyperbolic_2))
        println("\nHyperbolic model 2:")
        println("beta0 = ", value(beta00_H2))
        println("beta1 = ", value(beta01_H2))
        println("beta2 = ", value(beta02_H2))
        println("beta3 = ", value(beta03_H2))
    end


    if is_solved_and_feasible(model_hyperbolic_3)
        println("Objective value for hyperbolic 3 model: ", objective_value(model_hyperbolic_3))
        println("\nHyperbolic model 3:")
        println("beta0 = ", value(beta00_H3))
        println("beta1 = ", value(beta01_H3))
        println("beta2 = ", value(beta02_H3))
        println("beta3 = ", value(beta03_H3))
        println("beta4 = ", value(beta04_H3))

    end

    if is_solved_and_feasible(model_hyperbolic_4)
        println("Objective value for hyperbolic 4 model: ", objective_value(model_hyperbolic_4))
        println("\nHyperbolic model 4:")
        println("beta1 = ", value(beta01_H4))
        println("beta2 = ", value(beta02_H4))
        println("beta3 = ", value(beta03_H4))
        println("beta4 = ", value(beta04_H4))

    end

    if is_solved_and_feasible(model_logarithmic)
        println("Objective value for logarithmic model: ", objective_value(model_logarithmic))
        println("\nLogarithmic model:")
        println("beta0 = ", value(beta00_LG))
        println("beta1 = ", value(beta01_LG))
        println("beta2 = ", value(beta02_LG))

    end

    if is_solved_and_feasible(model_log_hyper)
        println("Objective value for combination model: ", objective_value(model_log_hyper))
        println("\nCombination model:")
        println("beta0 = ", value(beta00_LH))
        println("beta1 = ", value(beta01_LH))
        println("beta2 = ", value(beta02_LH))
        println("beta3 = ", value(beta03_LH))
        println("beta4 = ", value(beta04_LH))
        # println("beta5 = ", value(beta05_LH))
    end

    println(Pw_scaled)


end

using Plots
function Plot_surrogateModel()

    #points from the test data

    Pw_test = data.Power
    HTO_test = data.HTO

    Pw_validation = data_validation.Power
    HTO_validation = data_validation.HTO
    Pw = range(6000.0, 60000.0, length=100)

    #hyperbolic regression model
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

    # plot(Pw, HTO_first.(Pw), label="HTO = b0 + b1Pw + b2/(b3Pw + b4)", xlabel="Power (W)", ylabel="HTO (mol/s)", title="Surrogate Models", legend=:topright)
    plot(Pw_test, HTO_test, seriestype=:scatter, label="Test Data", xlabel="Power (W)",title="Surrogate Models", ylabel="HTO (mol/s)", legend=:topright)
    plot!(Pw_validation, HTO_validation, seriestype=:scatter, label="Validation Data", xlabel="Power (W)",title="Surrogate Models", ylabel="HTO (mol/s)", legend=:topright)
    # plot!(Pw, HTO_Hyper3.(Pw), label="HTO = b0 + b1Pw + b2/(b3Pw + b4), SSE = 3.97", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
        #Hyperbolic
    plot!(Pw, HTO_hyper2.(Pw), label="HTO_hyp2 = b0 + b1Pw + b2/(Pw + b3), SSE = 3.7e-6", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
    plot!(Pw, HTO_hyper1.(Pw), label="HTO_hyp1 = b0 + b1Pw + b2/(Pw), SSE = 0.0019 ", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
    plot!(Pw, HTO_hyper4.(Pw), label="HTO_hyp4 = b0 + b1Pw + b2/(b3Pw + b4), SSE = 0.0037", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
    # plot!(Pw, HTO_quadratic.(Pw), label="HTO = b0 + b1Pw + b2Pw^2, SSE = 1.1", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
        #logarithmic
    # plot!(Pw, HTO_logarithmic.(Pw), label="HTO_log = b0 + b1*log(Pw), SSE = 1.0238", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
    # plot!(Pw, HTO_logaritmic2.(Pw), label="HTO_log2 = b0 + b1*log(Pw) + b2*log(Pw)^2, SSE = 0.026", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
        #combination
    # plot!(Pw, HTO_combination_1.(Pw), label="HTO_com_1 = b0 + b1*log(Pw) + b2*log(Pw)^2 + b3/(Pw + b4), SSE = 7.52e-6", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
    # plot!(Pw, HTO_combination_2.(Pw), label="HTO_com_2 = b0 + b1*log(Pw) + b2/(Pw + b3), SSE = 0.00223", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)
    # plot!(Pw, HTO_combination_3.(Pw), label="HTO_com_3 = b0 + b1*log(Pw) + b2*log(Pw)^2 + b3/(b4*Pw + b5), SSE = 8.87e-9", xlabel="Power (W)", ylabel="HTO (mol/s)", legend=:topright)

end

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

function SurrogateModel_SeveralModels()
    data = CSV.read(joinpath("/home/annkar/NewCode/Data", "surrogate_data_Pw_T_training.csv"), DataFrame)

    
    Pw_train = data.Power
    T_train = data.Temperature
    HTO_train = data.HTO

    Pw_train_scaled = scaled(Pw_train)
    T_train_scaled = scaled(T_train)
    HTO_train_scaled = scaled(HTO_train)

    N = length(Pw_train_scaled)

    ##=====================================================================
    #defining different optimizer models
    ##=====================================================================
    model_hyperbolic_Pw_T_01 = Model(Ipopt.Optimizer)
    set_attribute(model_hyperbolic_Pw_T_01, "linear_solver", "mumps")


    ##=====================================================================
    #defining variables for the differet polynomial regression models
    ##=====================================================================
    @variable(model_hyperbolic_Pw_T_01, beta00)
    @variable(model_hyperbolic_Pw_T_01, beta01)
    @variable(model_hyperbolic_Pw_T_01, beta02)
    @variable(model_hyperbolic_Pw_T_01, beta03 >= 0.001)
    @variable(model_hyperbolic_Pw_T_01, beta04 >= 0.001)
    @variable(model_hyperbolic_Pw_T_01, HTO_01[1:N])
    #create the constraints for the polynomial regression model
    @constraint(model_hyperbolic_Pw_T_01, [i=1:N], beta00 + beta01*1/(Pw_train_scaled[i] +beta02) + beta03*(1/(T_train_scaled[i]+beta04)) + beta04*Pw_train_scaled[i]*T_train_scaled[i] == HTO_01[i])
    
    
    # defining the objective function for the polynomial regression model
    @objective(model_hyperbolic_Pw_T_01, Min, sum((HTO_01[i] - HTO_train_scaled[i])^2 for i in 1:length(Pw_train_scaled)))





    ##=====================================================================
    #solving the optimization probelm
    ##=====================================================================

    println("Optimizing the hyperbolic model...")
    optimize!(model_hyperbolic_Pw_T_01)

    ##=====================================================================
    #Cheking the solver status and printing the results
    ##=====================================================================

    println("Solver status for combination model: ", termination_status(model_hyperbolic_Pw_T_01))





    if is_solved_and_feasible(model_hyperbolic_Pw_T_01)
        println("Objective value for hyperbolic 1 model: ", objective_value(model_hyperbolic_Pw_T_01))
        println("\nHyperbolic model 1:")
        println("beta0 = ", value(beta00))
        println("beta1 = ", value(beta01))
        println("beta2 = ", value(beta02))
        println("beta3 = ", value(beta03))
        println("beta4 = ", value(beta04))

    end
    
end




using PlotlyJS, CSV, DataFrames
function Plot_surrogateModel_SeveralModels()
    data = CSV.read(joinpath("/home/annkar/NewCode/Data", "surrogate_data_Pw_T_training.csv"), DataFrame)

    Pw_train = data.Power
    T_train = data.Temperature
    HTO_train = data.HTO

    Pw_train_scaled = scaled(Pw_train)
    T_train_scaled = scaled(T_train)
    HTO_train_scaled = scaled(HTO_train)

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

    plotPolynomial = true

    if plotPolynomial
        Pw_grid = range(-0.3, 0.3, length=50)
        T_grid = range(-0.3, 0.3, length=50)

        # Define the polynomial regression model
    
        HTO_model_01(Pw, T) = -0.11426366984832503 -0.003920002733338371/(Pw -0.03953077268715382) + 0.0027439659308452847/(T + 0.09044948004792186)
        
        #create a meshgrid for HTO

        HTO_grid_01 = [HTO_model_01(Pw, T) for Pw in Pw_grid, T in T_grid]
        
        # Create 3D surface plots
        surface_trace01 = PlotlyJS.surface(x = Pw_grid, y = T_grid, z = HTO_grid_01,
            colorscale = "Viridis", showscale = true, colorbar = attr(title = "HTO"),
            name = "Surrogate Model"
        )

        
        
        # Layout

        trace_train = PlotlyJS.scatter3d(x = Pw_train_scaled, y = T_train_scaled, z = HTO_train_scaled, mode = "markers",
            marker = attr(size = 5, color = HTO_train_scaled, colorscale = "Viridis",  showscale = true, colorbar = attr(title = "HTO")
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

export fit_surrogate, Plot_surrogateModel, SurrogateModel_validationData, Plot_surrogateModel_SeveralModels, SurrogateModel_SeveralModels


