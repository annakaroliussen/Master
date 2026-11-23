using Plots

# data0 = CSV.read(joinpath("/home/annkar/NewCode/Data", "Plot_HTO_vs_varyPc.csv"), DataFrame)

# HTO_343 = data0.HTO[1:217]
# HTO_313 = data0.HTO[218:434]
# HTO_293 = data0.HTO[435:651]

# Pw_343 = data0.Power[1:217]*1e-3
# Pw_313 = data0.Power[218:434]*1e-3
# Pw_293 = data0.Power[435:651]*1e-3

# Plots.plot(Pw_343, HTO_343, label="T=343.15 K, Pw = 30 bar", linestyle=:dashdot)
# Plots.plot!(Pw_313, HTO_313, label="T=313.15 K, Pw = 30 bar", linestyle=:dash)
# Plots.plot!(Pw_293, HTO_293, label="T=293.15 K, Pw = 30 bar", xlabel="Power (kW)", ylabel="HTO (%)", title="HTO vs Power, const Pressure", legend=:topright)

# data1 = CSV.read(joinpath("/home/annkar/NewCode/Data", "Plot_HTO_vs_varyT.csv"), DataFrame)

# HTO_24bar = data1.HTO[1:217]
# HTO_17bar = data1.HTO[218:434] 
# HTO_10bar = data1.HTO[435:651]

# Pw_24bar = data1.Power[1:217]*1e-3
# Pw_17bar = data1.Power[218:434]*1e-3
# Pw_10bar = data1.Power[435:651]*1e-3

# Plots.plot(Pw_24bar, HTO_24bar, label="Pw = 24 bar, T = 333.15 K", linestyle=:dashdot)
# Plots.plot!(Pw_17bar, HTO_17bar, label="Pw = 17 bar, T = 333.15 K", linestyle=:dash)
# Plots.plot!(Pw_10bar, HTO_10bar, label="Pw = 10 bar, T = 333.15 K", xlabel="Power (kW)", ylabel="HTO (%)", title="HTO vs Power, const Temperature", legend=:topright)

# data2 = CSV.read(joinpath("/home/annkar/NewCode/Data", "Plot_HTO_vs_Pc.csv"), DataFrame)

# HTO_9 = data2.HTO[1:291]
# HTO_30 = data2.HTO[292:582]
# HTO_60 = data2.HTO[583:873]

# Pc_9 = data2.P_cathode[1:291]
# Pc_30 = data2.P_cathode[292:582]
# Pc_60 = data2.P_cathode[583:873]

# Plots.plot(Pc_9, HTO_9, label="Pw = 9 kW, T = 333.15 K", linestyle=:dashdot)
# Plots.plot!(Pc_30, HTO_30, label="Pw = 30 kW, T = 333.15 K", linestyle=:dash)
# Plots.plot!(Pc_60, HTO_60, label="Pw = 60 kW, T = 333.15 K", xlabel="Cathode Pressure (bar)", ylabel="HTO (%)", title="HTO vs Pressure, const Temperature", legend=:topleft)



data3 = CSV.read(joinpath("/home/annkar/NewCode/Data", "Plot_HTO_vs_T.csv"), DataFrame)

HTO_9 = data3.HTO[1:701]
HTO_30 = data3.HTO[702:1402]
HTO_60 = data3.HTO[1403:2103]

T_9 = data3.Temperature[1:701]
T_30 = data3.Temperature[702:1402]
T_60 = data3.Temperature[1403:2103]

Plots.plot(T_9, HTO_9, label="Pw = 9 kW, Pc = 30 bar", linestyle=:dashdot)
Plots.plot!(T_30, HTO_30, label="Pw = 30 kW, Pc = 30 bar", linestyle=:dash)
Plots.plot!(T_60, HTO_60, label="Pw = 60 kW, Pc = 30 bar", xlabel="Temperature (K)", ylabel="HTO (%)", title="HTO vs Temperature, const Pressure", legend=:topleft)