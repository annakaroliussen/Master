"""
Dynamic model of a proton exchange membrane (PEM) water electrolyzer stack.

The model is a lumped-parameter DAE describing a single representative cell scaled to
`Nc` series-connected cells. It couples three submodels: electrochemical (cell voltage),
mass transport (species balances, membrane water and gas crossover), and thermal (stack
temperature). Built with ModelingToolkit; connect the inputs and outputs below as
`RealInput`/`RealOutput` ports.

Input ports:
    - Pw        : Power (W)
    - z_c       : Valve opening, cathode (-)
    - z_a       : Valve opening, anode (-)
    - F_in      : Molar flow into anode (mol/s)
    - x_in_H2   : Inlet mole fraction, H2 (-)
    - x_in_O2   : Inlet mole fraction, O2 (-)
    - x_in_H2O  : Inlet mole fraction, H2O (-)
    - T_in      : Inlet temperature, anode (K)
    - P_env_c   : Downstream pressure, cathode (Pa)
    - P_env_a   : Downstream pressure, anode (Pa)

Output ports:
    - V         : Cell voltage (V)
    - F_out_c   : Molar flow out, cathode (mol/s)
    - F_out_a   : Molar flow out, anode (mol/s)
    - x_H2_c    : Outlet mole fraction, H2 cathode (-)
    - x_O2_c    : Outlet mole fraction, O2 cathode (-)
    - x_H2O_c   : Outlet mole fraction, H2O cathode (-)
    - x_H2_a    : Outlet mole fraction, H2 anode (-)
    - x_O2_a    : Outlet mole fraction, O2 anode (-)
    - x_H2O_a   : Outlet mole fraction, H2O anode (-)
    - T         : Stack temperature (K)
    - P_c       : Pressure, cathode (Pa)
    - P_a       : Pressure, anode (Pa)
    - HTO_a     : Hydrogen-to-oxygen content, anode (vol%)

Note: all pressures are in Pa and temperatures in K. Inputs are supplied through
connected source blocks; see `test_Electrolyzer` for an example and test that the 
electrolyzer model is working as expected.
"""

include("Parameters.jl")     # Load parameters for the electrolyzer model

using ModelingToolkit, OrdinaryDiffEq, ModelingToolkitStandardLibrary.Blocks, Plots, Serialization, XLSX, Dates

# Electrolyzer stack model --- par_El contains parameters specifing the electrolyzer stack: 
function Electrolyzer(par_El; name)

    @parameters t
    D = Differential(t)

    # Stack inputs:
    # ----------------------------------------------------
        @named Pw       = RealInput() # [W]     Power
        @named F_in     = RealInput() # [mol/s] Molar flow in anode
        @named x_in_H2  = RealInput() # [-]     Inlet mole fraction H2 (anode)
        @named x_in_O2  = RealInput() # [-]     Inlet mole fraction O2 (anode)
        @named x_in_H2O = RealInput() # [-]     Inlet mole fraction H2O (anode)
        @named T_in     = RealInput() # [K]     Temperature in anode
        @named P_env_c  = RealInput() # [Pa]    Pressure environment cathode
        @named P_env_a  = RealInput() # [Pa]    Pressure environment anode

    # Stack outputs:
    # ----------------------------------------------------
        @named V       = RealOutput() # [V]     Voltage
        @named F_out_c = RealOutput() # [mol/s] Molar flow out cathode
        @named F_out_a = RealOutput() # [mol/s] Molar flow out anode
        @named x_H2_c  = RealOutput() # [-]     Outlet mole fraction H2 (cathode)
        @named x_O2_c  = RealOutput() # [-]     Outlet mole fraction O2 (cathode)
        @named x_H2O_c = RealOutput() # [-]     Outlet mole fraction H2O (cathode)
        @named x_H2_a  = RealOutput() # [-]     Outlet mole fraction H2 (anode)
        @named x_O2_a  = RealOutput() # [-]     Outlet mole fraction O2 (anode)
        @named x_H2O_a = RealOutput() # [-]     Outlet mole fraction H2O (anode)
        @named T       = RealOutput() # [K]     Temperature
        @named P_c     = RealOutput() # [Pa]    Pressure cathode
        @named P_a     = RealOutput() # [Pa]    Pressure anode
        @named HTO_a   = RealOutput() # [vol%]  Hydrogen-to-oxygen content (anode)

    # Stack variables:
    # ----------------------------------------------------
    vars = ModelingToolkit.@variables begin

        # Thermodynamic constants:
            DrH_R(t)           # [J/mol]     Reaction enthalpy at reference T (298 K) 
            DrG_R(t)           # [J/mol]     Reaction Gibbs energy at reference T (298 K)
            Crp(t)             # [J/(mol*K)] Reaction heat-capacity change, dCp assumed independent of temperature
            DrH_0(t)           # [J/mol]     Reaction enthalpy at stack T
            DrG_0(t)           # [J/mol]     Reaction Gibbs energy at stack T

        # Electrochemical submodel:
            i_el(t)            # [A/m^2]     Current density 
            V_rev(t)           # [V]         Reversible cell voltage
            V_ohm(t)           # [V]         Ohmic overpotential 
            R_mem(t)           # [Ohm m^2]   Ionic resistance of membrane
            V_act(t)           # [V]         Activation overpotential
            V_act_a(t)         # [V]         Activation overpotential, anode
            i_a0(t)            # [A/m^2]     Exchange current density, anode
            V_con(t)           # [V]         Concentration overpotential
            V_con_c(t)         # [V]         Concentration overpotential, cathode
            V_con_a(t)         # [V]         Concentration overpotential, anode

        # Molar inventories:
            n_c(t)             # [mol]       Total molar inventory in cathode compartment
            n_a(t)             # [mol]       Total molar inventory in anode compartment
            n_H2_c(t)          # [mol/s]     Molar inventory of H2 (cathode)
            n_O2_c(t)          # [mol/s]     Molar inventory of O2 (cathode)
            n_H2O_c(t)         # [mol/s]     Molar inventory of H2O (cathode)
            n_H2_a(t)          # [mol/s]     Molar inventory of H2 (anode)
            n_O2_a(t)          # [mol/s]     Molar inventory of O2 (anode)
            n_H2O_a(t)         # [mol/s]     Molar inventory of H2O (anode)

        # Production terms:
            n_H2_p(t)          # [mol/s]     Production rate of H2
            n_O2_p(t)          # [mol/s]     Production rate of O2
            n_H2O_p(t)         # [mol/s]     Production rate of H2O
        
        # Water crossover:
            C_H2O_c(t)         # [mol/m^3]   Concentration of water in cathode compartment
            C_H2O_a(t)         # [mol/m^3]   Concentration of water in anode compartment
            C_H2O_c_mem(t)     # [mol/m^3]   Concentration of water at cathode membrane surface
            C_H2O_a_mem(t)     # [mol/m^3]   Concentration of water at anode membrane surface
            n_H2O_diff(t)      # [mol/s]     Water diffusion through membrane
            n_H2O_eod(t)       # [mol/s]     Water electro-osmotic drag through membrane
            n_H2O_dp(t)        # [mol/s]     Water transport by pressure gradient (Darcy) through membrane
            n_H2O_cross(t)     # [mol/s]     Total water crossover through membrane

        # Gas crossover:
            C_H2_c(t)          # [mol/m^3]   Concentration of H2 in cathode compartment 
            C_H2_a(t)          # [mol/m^3]   Concentration of H2 in anode compartment
            C_O2_c(t)          # [mol/m^3]   Concentration of O2 in cathode compartment
            C_O2_a(t)          # [mol/m^3]   Concentration of O2 in anode compartment
            C_H2_c_mem(t)      # [mol/m^3]   Concentration of H2 at cathode membrane surface
            C_H2_a_mem(t)      # [mol/m^3]   Concentration of H2 at anode membrane surface
            C_O2_c_mem(t)      # [mol/m^3]   Concentration of O2 at cathode membrane surface
            C_O2_a_mem(t)      # [mol/m^3]   Concentration of O2 at anode membrane surface
            zeta(t)            # [-]         Empirical crossover factor (temperature-dependent); see E. Afshari et al. 2020
            n_H2_diff(t)       # [mol/s]     H2 crossover by diffusion
            n_H2_eod(t)        # [mol/s]     H2 crossover by current-driven
            n_H2_cross(t)      # [mol/s]     Total H2 crossover
            n_O2_diff(t)       # [mol/s]     O2 crossover by diffusion
            n_O2_eod(t)        # [mol/s]     O2 crossover by current-driven
            n_O2_cross(t)      # [mol/s]     Total O2 crossover

            #oxygen and hydrogen permeability
            #epsilon_H2(t)        # [mol/(m s Pa)] Hydrogen permeability
            #epsilon_O2(t)        # [mol/(m s Pa)] Oxygen permeability

        # Diffusion coefficients:
            D_H(t)             # [m^2/s]     Binary diffusion coefficient H2 - H2O
            D_O(t)             # [m^2/s]     Binary diffusion coefficient O2 - H2O
            D_eff_H(t)         # [m^2/s]     Effective binary diffusion coefficient H2 - H2O
            D_eff_O(t)         # [m^2/s]     Effective binary diffusion coefficient O2 - H2O

        # Pressure:
            P_H2O(t)           # [Pa]        Saturation pressure of water
            n_H2O_c_g(t)       # [mol]       Molar inventory of water vapor in cathode compartment
            n_H2O_a_g(t)       # [mol]       Molar inventory of water vapor in anode compartment
            P_c_H2_mem(t)      # [Pa]        Partial pressure of H2 at cathode membrane surface
            P_c_O2_mem(t)      # [Pa]        Partial pressure of O2 at cathode membrane surface
            P_a_H2_mem(t)      # [Pa]        Partial pressure of H2 at anode membrane surface
            P_a_O2_mem(t)      # [Pa]        Partial pressure of O2 at anode membrane surface

        # Molar inlet and outlet flows:
            n_H2_in_a(t)       # [mol/s]     Molar inlet flow H2 (anode)
            n_O2_in_a(t)       # [mol/s]     Molar inlet flow O2 (anode)
            n_H2O_in_a(t)      # [mol/s]     Molar inlet flow H2O (anode)
            n_H2_out_c(t)      # [mol/s]     Molar outlet flow H2 (cathode)
            n_O2_out_c(t)      # [mol/s]     Molar outlet flow O2 (cathode)
            n_H2O_out_c(t)     # [mol/s]     Molar outlet flow H2O (cathode)
            n_H2_out_a(t)      # [mol/s]     Molar outlet flow H2 (anode)
            n_O2_out_a(t)      # [mol/s]     Molar outlet flow O2 (anode)
            n_H2O_out_a(t)     # [mol/s]     Molar outlet flow H2O (anode)

        # Gas compositions:
            y_H2_c(t)          # [-]         Mole fraction of H2 in cathode compartment
            y_O2_c(t)          # [-]         Mole fraction of O2 in cathode compartment
            y_H2O_c(t)         # [-]         Mole fraction of H2O in cathode compartment
            y_H2_a(t)          # [-]         Mole fraction of H2 in anode compartment
            y_O2_a(t)          # [-]         Mole fraction of O2 in anode compartment
            y_H2O_a(t)         # [-]         Mole fraction of H2O in anode compartment

        # PTL:
            V_PTL_c(t)         # [m^3]       Volume of porous transport layer (cathode)
            V_PTL_a(t)         # [m^3]       Volume of porous transport layer (anode)
            C_H2_PTL_c(t)      # [mol/m^3]   Concentration of H2 in porous transport layer (cathode)
            C_O2_PTL_c(t)      # [mol/m^3]   Concentration of O2 in porous transport layer (cathode)
            C_H2O_PTL_c(t)     # [mol/m^3]   Concentration of H2O in porous transport layer (cathode)
            C_H2_PTL_a(t)      # [mol/m^3]   Concentration of H2 in porous transport layer (anode)
            C_O2_PTL_a(t)      # [mol/m^3]   Concentration of O2 in porous transport layer (anode)
            C_H2O_PTL_a(t)     # [mol/m^3]   Concentration of H2O in porous transport layer (anode)
            n_H2_PTL_c(t)      # [mol]       Molar inventory of H2 in porous transport layer (cathode)
            n_O2_PTL_c(t)      # [mol]       Molar inventory of O2 in porous transport layer (cathode)
            n_H2O_PTL_c(t)     # [mol]       Molar inventory of H2O in porous transport layer (cathode)
            n_H2_PTL_a(t)      # [mol]       Molar inventory of H2 in porous transport layer (anode)
            n_O2_PTL_a(t)      # [mol]       Molar inventory of O2 in porous transport layer (anode)
            n_H2O_PTL_a(t)     # [mol]       Molar inventory of H2O in porous transport layer (anode)

        # Thermal submodel:
            h_H2_in_a(t)       # [J/mol]     Enthalpy of H2 inlet flow (anode)
            h_O2_in_a(t)       # [J/mol]     Enthalpy of O2 inlet flow (anode)
            h_H2O_in_a(t)      # [J/mol]     Enthalpy of H2O inlet flow (anode)
            Cth(t)             # [J/K]       Heat capacity of stack inventories
            H_in(t)            # [W]         Enthalpy of inlet flows
            Q_gen(t)           # [W]         Heat generation

        # Oxygen to hydrogen ratio cathode:
            OTH_c(t)           # [vol%]      Oxygen-to-hydrogen content (cathode)
    end

    # Equations:
    # ----------------------------------------------------
    eqs = [

        # Thermodynamic parameters:
        # ----------------------------------------------------
            # Reaction properties at reference T (298 K) and pressure (10^5 Pa), from species formation values:
                0 ~ - DrG_R - par_H2O[:DG_H2O] + par_H2[:DG_H2] + 1/2*par_O2[:DG_O2] # [J/mol]     Reaction Gibbs free energy 
                0 ~ - DrH_R - par_H2O[:DH_H2O] + par_H2[:DH_H2] + 1/2*par_O2[:DH_O2] # [J/mol]     Reaction enthalpy 
                0 ~ - Crp   - par_H2O[:Cp_H2O] + par_H2[:Cp_H2] + 1/2*par_O2[:Cp_O2] # [J/(mol*K)] Reaction heat-capacity change 

            # Temperature correction enthalpy (constant dCp):
                0 ~ - DrH_0 + DrH_R + Crp*(T.u - 298) # [J/mol] Reaction enthalpy at T and 10^5 Pa 

            # Reaction Gibbs energy at stack T (van't Hoff integrated with constant dCp):
                0 ~ - DrG_0 + T.u/298*(DrG_R - DrH_R) + DrH_R + Crp*(T.u - 298 - T.u*log(T.u/298)) # [J/mol] Reaction Gibbs free energy at T and 10^5 Pa 

        # Electrochemical submodel:
        # ----------------------------------------------------
            # Stack power input: 
                0 ~ - Pw.u + (V.u*par_El[:Nc])*i_el*par_El[:A_cell] # [W]
        
            # Stack voltage:
                0 ~ - V.u + V_rev + V_ohm + V_act + V_con # [V]
        
            # Reversible cell voltage:
                0 ~ - V_rev + DrG_0/(par_El[:z]*par[:F]) + par[:R]*T.u/(par_El[:z]*par[:F])*log(abs((P_c.u*y_H2_c)*(abs(P_a.u*y_O2_a))^(1/2)/(1e5)^1.5)) # [V] (1e5 = P0, standard pressure)

            # Ohmic overpotential:
                0 ~ - V_ohm + (R_mem + par_El[:R_el])*i_el # [V]
                0 ~ - R_mem + par_El[:delta_mem]/((0.5139*par_El[:lambda]-0.326)*exp(1268*(1/303 - 1/T.u))) # [Ohm m^2] Springer conductivity (303 K ref)

            # Activation overpotential:
                0 ~ - V_act + V_act_a # [V]
                0 ~ - i_el + i_a0*exp(par_El[:alpha_a]*par_El[:z]*par[:F]*V_act_a/(par[:R]*T.u)) # [A/m^2] Butler-Volmer (Tafel), anode
                0 ~ - i_a0 + par_El[:i_a0_ref]*exp(- par_El[:E_a]/par[:R]*(1/T.u - 1/328.15)) # [A/m^2] Exhange current density at T (anode)

            # Concentration overpotential:
                0 ~ - V_con + V_con_c + V_con_a # [V]
                0 ~ - V_con_c + par[:R]*T.u/(2*par[:F])*log(abs(C_H2_c_mem/C_H2_c)) # [V] Concentration overpotential (cathode)
                0 ~ - V_con_a + par[:R]*T.u/(4*par[:F])*log(abs(C_O2_a_mem/C_O2_a)) # [V] Concentration overpotential (anode)

        # Mass transfer submodel:
        # ----------------------------------------------------
            # Total molar inventories in cathode and anode compartments:
                0 ~ - n_c + n_H2_c + n_O2_c + n_H2O_c # [mol] Cathode
                0 ~ - n_a + n_H2_a + n_O2_a + n_H2O_a # [mol] Anode

            # Molar balances in cathode and anode compartments:
                D(n_H2_c)  ~ 0          - n_H2_out_c  + n_H2_p  - n_H2_cross  # [mol/s] Hydrogen (cathode)
                D(n_O2_c)  ~ 0          - n_O2_out_c            + n_O2_cross  # [mol/s] Oxygen (cathode)
                D(n_H2O_c) ~ 0          - n_H2O_out_c           + n_H2O_cross # [mol/s] Water (cathode)

                D(n_H2_a)  ~ n_H2_in_a  - n_H2_out_a            + n_H2_cross  # [mol/s] Hydrogen (anode)
                D(n_O2_a)  ~ n_O2_in_a  - n_O2_out_a  + n_O2_p  - n_O2_cross  # [mol/s] Oxygen (anode)
                D(n_H2O_a) ~ n_H2O_in_a - n_H2O_out_a + n_H2O_p - n_H2O_cross # [mol/s] Water (anode)

            # Production rates:
                0 ~ - n_H2_p  + (i_el*par_El[:A_stack])/(par_El[:z]*par[:F])*par_El[:eta_F]   # [mol/s] Hydrogen production rate
                0 ~ - n_O2_p  + (i_el*par_El[:A_stack])/(2*par_El[:z]*par[:F])*par_El[:eta_F] # [mol/s] Oxygen production rate (half of H2 production rate)
                0 ~ - n_H2O_p - (i_el*par_El[:A_stack])/(par_El[:z]*par[:F])*par_El[:eta_F]   # [mol/s] Water consumption rate (negative because water is consumed)

            # Water crossover:
                0 ~ - n_H2O_diff  + par_El[:Dw]/(par_El[:delta_mem])*(par_El[:delta_c]*n_H2O_cross/D_eff_H + par_El[:delta_a]*(n_H2O_cross - n_H2O_p)/D_eff_O) # [mol/s] Water diffusion through membrane (Fick's law)
                0 ~ - n_H2O_eod   + (par_El[:nd]*i_el*par_El[:A_stack])/(par[:F]) # [mol/s] Electro-osmotic drag of water through membrane
                0 ~ - n_H2O_dp    + (par_El[:K_Darcy]*par_El[:A_stack]*par_H2O[:rho_H2O]*((P_c.u)-(P_a.u)))/(par_El[:delta_mem]*par_H2O[:mu_H2O]*par_H2O[:M_H2O]) # [mol/s] Pressure-driven water transport through membrane (Darcy's law)
                0 ~ - n_H2O_cross + n_H2O_diff + n_H2O_eod - n_H2O_dp # [mol/s] Total water crossover through membrane

            # Gas crossover:
              #  0 ~ - n_H2_diff  + (par_El[:e_H2])*(P_c_H2_mem - P_a_H2_mem)*par_El[:A_stack]/par_El[:delta_mem] # [mol/s] H2 crossover by diffusion

                0 ~ - n_H2_diff  + (5.32*1e-14 * exp(-2405*(1/T.u -1/353.15)))  *(P_c_H2_mem - P_a_H2_mem)*par_El[:A_stack]/par_El[:delta_mem] # [mol/s] H2 crossover by diffusion
                0 ~ - n_H2_eod   - (i_el/par[:F])*zeta*(P_c_H2_mem*par_El[:Su_H2]/par_El[:CH2O])*par_El[:A_stack] # [mol/s] H2 crossover by current-driven convection (dissolved and dragged along with the water)
                0 ~ - n_H2_cross + n_H2_diff + n_H2_eod # [mol/s] Total H2 crossover
                
               # 0 ~ - n_O2_diff  + (par_El[:e_O2])*(P_a_O2_mem - P_c_O2_mem)*par_El[:A_stack]/par_El[:delta_mem] # [mol/s] O2 crossover by diffusion
                0 ~ - n_O2_diff  + (2.52*1e-14 * exp(-2034.75*(1/T.u -1/353.15)))*(P_a_O2_mem - P_c_O2_mem)*par_El[:A_stack]/par_El[:delta_mem] # [mol/s] O2 crossover by diffusion
                0 ~ - n_O2_eod   + (i_el/par[:F])*zeta*(P_a_O2_mem*par_El[:Su_O2]/par_El[:CH2O])*par_El[:A_stack] # [mol/s] O2 crossover by current-driven convection (dissolved and dragged along with the water)
                0 ~ - n_O2_cross + n_O2_diff + n_O2_eod # [mol/s] Total O2 crossover

                0 ~ - zeta + 0.0134*T.u + 0.03 # [-] Empirical crossover factor (temperature-dependent); see E. Afshari et al. 2020

            # Binary diffusion coefficients:
                0 ~ - (P_c.u/1.013e5)*(D_H*1e4) + par_El[:a]*(abs(T.u)/(sqrt(par_H2[:T_c_H2]*par_H2O[:T_c_H2O])))^par_El[:b]*(par_H2[:P_c_H2]*par_H2O[:P_c_H2O])^(1/3)*(par_H2[:T_c_H2]*par_H2O[:T_c_H2O])^(5/12)*(1/(par_H2[:M_H2]*1e3) + 1/(par_H2O[:M_H2O]*1e3))^(1/2) # [m^2/s] H2-H2O binary diffusivity
                0 ~ - (P_a.u/1.013e5)*(D_O*1e4) + par_El[:a]*(abs(T.u)/(sqrt(par_O2[:T_c_O2]*par_H2O[:T_c_H2O])))^par_El[:b]*(par_O2[:P_c_O2]*par_H2O[:P_c_H2O])^(1/3)*(par_O2[:T_c_O2]*par_H2O[:T_c_H2O])^(5/12)*(1/(par_O2[:M_O2]*1e3) + 1/(par_H2O[:M_H2O]*1e3))^(1/2) # [m^2/s] O2-H2O binary diffusivity 
            
            # Effective diffusivity: porosity correction
                0 ~ - D_eff_H + (D_H * par_El[:epsilon]*((par_El[:epsilon] - par_El[:epsilon_p])/(1 - par_El[:epsilon_p]))^par_El[:alpha]) # [m^2/s] Effective H2-H2O
                0 ~ - D_eff_O + (D_O * par_El[:epsilon]*((par_El[:epsilon] - par_El[:epsilon_p])/(1 - par_El[:epsilon_p]))^par_El[:alpha]) # [m^2/s] Effective O2-H2O

            # Pressure:
                0 ~ - P_H2O + 610.78*exp(17.2694*((T.u - 273.15)/(T.u - 34.85))) # [Pa] Saturation pressure of water (Tetens equation)

                0 ~ - P_H2O + y_H2O_c*P_c.u # [Pa] Impose saturation: cathode water partial pressure = P_H2O
                0 ~ - P_H2O + y_H2O_a*P_a.u # [Pa] Impose saturation: anode water partial pressure = P_H2O

                0 ~ - P_c.u + (n_H2_c + n_O2_c + n_H2O_c_g)*par[:R]*T.u/(par_El[:V_c] - (n_H2O_c - n_H2O_c_g)*par_H2O[:M_H2O]/par_H2O[:rho_H2O]) # [Pa] Pressure in cathode compartment
                0 ~ - P_a.u + (n_H2_a + n_O2_a + n_H2O_a_g)*par[:R]*T.u/(par_El[:V_a] - (n_H2O_a - n_H2O_a_g)*par_H2O[:M_H2O]/par_H2O[:rho_H2O]) # [Pa] Pressure in anode compartment

            # Pressure at membrane surface as given by Schalenbach et al.:
                0 ~ - P_c_H2_mem + P_c.u + par_El[:A_H2]*i_el - P_H2O
                0 ~ - P_a_H2_mem + P_c_H2_mem*(n_H2_cross/(n_O2_p + n_H2_cross))
                0 ~ - P_c_O2_mem + 0
                0 ~ - P_a_O2_mem + P_a.u + par_El[:A_O2]*i_el - P_H2O

            # Molar inlet and outlet flows:
                0 ~ - n_H2_in_a  + F_in.u*x_in_H2.u # [mol/s] Molar flow H2 in anode
                0 ~ - n_O2_in_a  + F_in.u*x_in_O2.u # [mol/s] Molar flow O2 in anode
                0 ~ - n_H2O_in_a + F_in.u*x_in_H2O.u # [mol/s] Molar flow H2O in anode

                0 ~ - F_out_c.u + par_El[:h_c]*171.55768399498805*sqrt(max(P_c.u - P_env_c.u, 0)) # [mol/s] Cathode outlet valve (0 if P_c <= P_env)
                0 ~ - F_out_a.u + par_El[:h_a]*244.80834801269452*sqrt(max(P_a.u - P_env_a.u, 0)) # [mol/s] Anode outlet valve (0 if P_a <= P_env)

                0 ~ - n_H2_out_c  + F_out_c.u*x_H2_c.u  # [mol/s] Molar flow H2 out cathode
                0 ~ - n_O2_out_c  + F_out_c.u*x_O2_c.u  # [mol/s] Molar flow O2 out cathode
                0 ~ - n_H2O_out_c + F_out_c.u*x_H2O_c.u # [mol/s] Molar flow H2O out cathode

                0 ~ - n_H2_out_a  + F_out_a.u*x_H2_a.u  # [mol/s] Molar flow H2 out anode
                0 ~ - n_O2_out_a  + F_out_a.u*x_O2_a.u  # [mol/s] Molar flow O2 out anode
                0 ~ - n_H2O_out_a + F_out_a.u*x_H2O_a.u # [mol/s] Molar flow H2O out anode

            # Mole fractions and molar gas fractions:
                0 ~ - x_H2_c.u*n_c + n_H2_c   # [-] Mole fraction of H2 in cathode compartment
                0 ~ - x_O2_c.u*n_c + n_O2_c   # [-] Mole fraction of O2 in cathode compartment
                0 ~ - x_H2O_c.u*n_c + n_H2O_c # [-] Mole fraction of H2O in cathode compartment

                0 ~ - x_H2_a.u*n_a + n_H2_a   # [-] Mole fraction of H2 in anode compartment
                0 ~ - x_O2_a.u*n_a + n_O2_a   # [-] Mole fraction of O2 in anode compartment
                0 ~ - x_H2O_a.u*n_a + n_H2O_a # [-] Mole fraction of H2O in anode compartment
                
                0 ~ - y_H2_c*(n_H2_c + n_O2_c + n_H2O_c_g) + n_H2_c     # [-] Molar gas fraction of H2 in cathode compartment
                0 ~ - y_O2_c*(n_H2_c + n_O2_c + n_H2O_c_g) + n_O2_c     # [-] Molar gas fraction of O2 in cathode compartment
                0 ~ - y_H2O_c*(n_H2_c + n_O2_c + n_H2O_c_g) + n_H2O_c_g # [-] Molar gas fraction of H2O in cathode compartment

                0 ~ - y_H2_a*(n_H2_a + n_O2_a + n_H2O_a_g) + n_H2_a     # [-] Molar gas fraction of H2 in anode compartment
                0 ~ - y_O2_a*(n_H2_a + n_O2_a + n_H2O_a_g) + n_O2_a     # [-] Molar gas fraction of O2 in anode compartment
                0 ~ - y_H2O_a*(n_H2_a + n_O2_a + n_H2O_a_g) + n_H2O_a_g # [-] Molar gas fraction of H2O in anode compartment

            # Hydrogen-to-oxygen content and oxygen-to-hydrogen content:
                0 ~ - OTH_c*(y_H2_c + y_O2_c) + y_O2_c*100   # [vol%] Oxygen-to-hydrogen content (cathode)
                0 ~ - HTO_a.u*(y_H2_a + y_O2_a) + y_H2_a*100 # [vol%] Hydrogen-to-oxygen content (anode)

        # Thermal submodel:
        # ----------------------------------------------------
            # Heat capacity of stack inventories: (It is assumed all water in the PTL is liquid)
                0 ~ - Cth + par_H2[:Cp_H2]*(n_H2_c + n_H2_PTL_c + n_H2_a + n_H2_PTL_a) + par_O2[:Cp_O2]*(n_O2_c + n_O2_PTL_c + n_O2_a + n_O2_PTL_a) + par_H2O[:Cp_H2O]*(n_H2O_c - n_H2O_c_g + n_H2O_PTL_c + n_H2O_a - n_H2O_a_g + n_H2O_PTL_a) + par_H2O[:Cp_H2O_g]*(n_H2O_c_g + n_H2O_a_g)

            # Enthalpy of inlet flows:
                0 ~ - h_H2_in_a + par_H2[:Cp_H2]*(T_in.u  - T.u)   # [J/mol] Enthalpy of inlet flows hydrogen (anode)
                0 ~ - h_O2_in_a + par_O2[:Cp_O2]*(T_in.u  - T.u)   # [J/mol] Enthalpy of inlet flows oxygen (anode)
                0 ~ - h_H2O_in_a + par_H2O[:Cp_H2O]*(T_in.u - T.u) # [J/mol] Enthalpy of inlet flows water (anode)

            # Total enthalpy of inlet flows:
                0 ~ - H_in + h_H2_in_a*n_H2_in_a + h_O2_in_a*n_O2_in_a + h_H2O_in_a*n_H2O_in_a # [W]

            # Heat generation term:
                0 ~ - Q_gen + (V.u - (DrH_0)/(par_El[:z]*par[:F]))*(i_el*par_El[:A_stack]) # [W]

            # Energy balance of stack: (Lumped inventories and assumed adiabatic)
                D(T.u) ~ 1/Cth*(H_in + Q_gen) # [K/s]

            # PTL:
                # Volumes:
                    0 ~ - V_PTL_c + par_El[:A_stack]*par_El[:epsilon]*par_El[:delta_c] # [m^3] Volume of porous transport layer (cathode)
                    0 ~ - V_PTL_a + par_El[:A_stack]*par_El[:epsilon]*par_El[:delta_a] # [m^3] Volume of porous transport layer (anode)

                # Chanel consentrations:
                    0 ~ - C_H2O_c + par_H2O[:rho_H2O]/par_H2O[:M_H2O] # [mol/m^3] Water concentration in cathode compartment
                    0 ~ - C_H2O_a + par_H2O[:rho_H2O]/par_H2O[:M_H2O] # [mol/m^3] Water concentration in anode compartment
                                    
                    0 ~ - C_H2_c + (7.8e-6 * exp(530*(1/T.u - 1/298.15)))*P_c.u*y_H2_c # [mol/m^3] H2 concentration in cathode compartment
                    0 ~ - C_O2_c + (1.3e-5 * exp(1500*(1/T.u - 1/298.15)))*P_c.u*y_O2_c # [mol/m^3] O2 concentration in cathode compartment
                    0 ~ - C_H2_a + (7.8e-6 * exp(530*(1/T.u - 1/298.15)))*P_a.u*y_H2_a # [mol/m^3] H2 concentration in anode compartment
                    0 ~ - C_O2_a + (1.3e-5 * exp(1500*(1/T.u - 1/298.15)))*P_a.u*y_O2_a # [mol/m^3] O2 concentration in anode compartment

                # Consentrations at membrane surface:
                    0 ~ - C_H2O_c_mem + C_H2O_c + par_El[:delta_c]/(par_El[:A_stack]*D_eff_H)*(n_H2O_cross) # [mol/m^3] Water concentration at cathode membrane surface
                    0 ~ - C_H2O_a_mem + C_H2O_a - par_El[:delta_a]/(par_El[:A_stack]*D_eff_O)*(n_H2O_cross - n_H2O_p) # [mol/m^3] Water concentration at anode membrane surface

                    0 ~ - C_H2_c_mem + C_H2_c + par_El[:delta_c]/(par_El[:A_stack]*D_eff_H)*(n_H2_p - n_H2_cross) # [mol/m^3] H2 concentration at cathode membrane surface
                    0 ~ - C_O2_c_mem + C_O2_c + par_El[:delta_c]/(par_El[:A_stack]*D_eff_O)*(n_O2_cross) # [mol/m^3] O2 concentration at cathode membrane surface
                    0 ~ - C_H2_a_mem + C_H2_a + par_El[:delta_a]/(par_El[:A_stack]*D_eff_H)*(n_H2_cross) # [mol/m^3] H2 concentration at anode membrane surface
                    0 ~ - C_O2_a_mem + C_O2_a + par_El[:delta_a]/(par_El[:A_stack]*D_eff_O)*(n_O2_p - n_O2_cross) # [mol/m^3] O2 concentration at anode membrane surface

                # Average concentrations:
                    0 ~ - C_H2_PTL_c + (C_H2_c_mem + C_H2_c)/2    # [mol/m^3] Average concentration of H2 in PTL (cathode)
                    0 ~ - C_O2_PTL_c + (C_O2_c_mem + C_O2_c)/2    # [mol/m^3] Average concentration of O2 in PTL (cathode)
                    0 ~ - C_H2O_PTL_c + (C_H2O_c_mem + C_H2O_c)/2 # [mol/m^3] Average concentration of H2O in PTL (cathode)
                    0 ~ - C_H2_PTL_a + (C_H2_a_mem + C_H2_a)/2    # [mol/m^3] Average concentration of H2 in PTL (anode)
                    0 ~ - C_O2_PTL_a + (C_O2_a_mem + C_O2_a)/2    # [mol/m^3] Average concentration of O2 in PTL (anode)
                    0 ~ - C_H2O_PTL_a + (C_H2O_a_mem + C_H2O_a)/2 # [mol/m^3] Average concentration of H2O in PTL (anode)

                # Molar inventory PTL:
                    0 ~ - n_H2_PTL_c + V_PTL_c*C_H2_PTL_c   # [mol] Molar inventory of H2 in porous transport layer (cathode)
                    0 ~ - n_O2_PTL_c + V_PTL_c*C_O2_PTL_c   # [mol] Molar inventory of O2 in porous transport layer (cathode)
                    0 ~ - n_H2O_PTL_c + V_PTL_c*C_H2O_PTL_c # [mol] Molar inventory of H2O in porous transport layer (cathode)
                    0 ~ - n_H2_PTL_a + V_PTL_a*C_H2_PTL_a   # [mol] Molar inventory of H2 in porous transport layer (anode)
                    0 ~ - n_O2_PTL_a + V_PTL_a*C_O2_PTL_a   # [mol] Molar inventory of O2 in porous transport layer (anode)
                    0 ~ - n_H2O_PTL_a + V_PTL_a*C_H2O_PTL_a # [mol] Molar inventory of H2O in porous transport layer (anode)
    ]

    # Finilize and return:
    ports = [
        Pw, F_in, T_in, x_in_H2, x_in_O2, x_in_H2O, # Inlet ports
        V, F_out_c, F_out_a, x_H2_c, x_O2_c, x_H2O_c, x_H2_a, x_O2_a, x_H2O_a, P_c, P_a, T, P_env_c, P_env_a, HTO_a # Outlet ports
    ]
    # ----------------------------------------------------
    compose(ODESystem(eqs, t, vars, []; name), ports)
end


# Function to test the electrolyzer model:
function test_Electrolyzer(save_ss=false, make_plot=false; P_cathode_bar=30.0, P_anode_bar=1.0, T_stack_K=353.15, power=60000.0, use_T_ctrl = true,use_Pc_ctrl = true, use_Pa_ctrl = true, init_state=nothing)

    # Create electrolyzer model:
        @named el       = Electrolyzer(par_El)

    # Create input blocks:
        @named Pw       = Blocks.Constant(k = power) # [w] Input power
        # @named Pw       = Blocks.Constant(k = 60000) # [w] Input power
        @named F_in     = Blocks.Constant(k = 35.09) # [mol/s] Molar inlet flow
        # @named z_a      = Blocks.Ramp(height = 100.0, offset = 244.80834801269452, start_time = 2, duration = 110, smooth = true) 
        @named x_in_H2  = Blocks.Constant(k = 0) # [-] Composition in (hydrogen)
        @named x_in_O2  = Blocks.Constant(k = 0) # [-] Composition in (oxygen)
        @named x_in_H2O = Blocks.Constant(k = 1) # [-] Composition in (water)
        
        
    # For my sensistivity analysis, looking at changes in T, Pa and Pc
    
        #  @named T_ramp = Blocks.Ramp(height = 40.0, offset = 293.15,start_time = 2,duration = 110, smooth = true)

    # Controller blocks:
        # --- Stack temperature loop (setpoint 80 °C) --- 
            use_T_ctrl = use_T_ctrl   # <-- true = PI controller, false = step

            # Simple placeholder tuning (posetive gain: feedwater temperature drops -> stack temperature drops):
            Kc_T     = 50      # [K/K] proportional gain
            tau_i_T  = 8.0         # [s] integral time

            @named T_sp   = Blocks.Constant(k = T_stack_K)   # temperature setpoint [K]
            # @named T_sp = Blocks.Ramp(height = 40.0, offset = 293.15,start_time = 2,duration = 110, smooth = true)   
            @named T_ctrl = Controller(use_controller = use_T_ctrl,
                                       Kc = Kc_T, Ti = tau_i_T, u_min = 0.0, u_max = 400.0,
                                       int_x0 = 332.15/Kc_T,
                                       step_offset = 332.15, step_height = 30.0,
                                       step_start = 10.0)

        # --- Cathode pressure loop (setpoint set by P_cathode_bar) ---
            use_Pc_ctrl = use_Pc_ctrl

            # Simple placeholder tuning (negative gain: valve opens -> pressure drops):
            Kc_Pc  = 1           # [valve unit / Pa]  negative gain
            Ti_Pc  = 20.0         # [s] integral time

            P_env_c_nom = 29.9e5   # nominal cathode valve opening (operating point)

            @named Pc_sp   = Blocks.Constant(k = P_cathode_bar * 1e5) # [Pa]
            @named Pc_ctrl = Controller(use_controller = use_Pc_ctrl,
                                        Kc = Kc_Pc, Ti = Ti_Pc,
                                        u_min = 0.5e5, u_max = 30e5,   # valve opening range
                                        int_x0 = P_env_c_nom / Kc_Pc,      # bumpless start
                                        step_offset = P_env_c_nom, step_height = 0.0,
                                        step_start = 10.0)

        # --- Anode pressure loop (setpoint 1 bar) ---
            use_Pa_ctrl = use_Pa_ctrl

            Kc_Pa  = 1
            Ti_Pa  = 20.0

            P_env_a_nom = 0.9e5   # nominal anode valve opening

            @named Pa_sp   = Blocks.Constant(k = P_anode_bar * 1e5)    # 1 bar setpoint [Pa]
            @named Pa_ctrl = Controller(use_controller = use_Pa_ctrl,
                                        Kc = Kc_Pa, Ti = Ti_Pa,
                                        u_min = 0.5e5, u_max = 4e5,
                                        int_x0 = P_env_a_nom / Kc_Pa,
                                        step_offset = P_env_a_nom, step_height = 0.0,
                                        step_start = 10.0)

    # Connect blocks:
    connections = [
        # Stack inputs: 
            ModelingToolkit.connect(Pw.output,       el.Pw) # Power
            ModelingToolkit.connect(F_in.output,     el.F_in) # Molar inlet flow
            ModelingToolkit.connect(x_in_H2.output,  el.x_in_H2) # Composition (hydrogen)
            ModelingToolkit.connect(x_in_O2.output,  el.x_in_O2) # Composition (oxygen)
            ModelingToolkit.connect(x_in_H2O.output, el.x_in_H2O) # Compostion (water)
            
        # Controller connections:
            # Temperature:
                ModelingToolkit.connect(T_sp.output,     T_ctrl.setpoint)
                ModelingToolkit.connect(el.T,            T_ctrl.measurement)  
                ModelingToolkit.connect(T_ctrl.actuator, el.T_in)          

            # Cathode pressure loop:
                ModelingToolkit.connect(Pc_sp.output,     Pc_ctrl.setpoint)
                ModelingToolkit.connect(el.P_c,           Pc_ctrl.measurement)
                ModelingToolkit.connect(Pc_ctrl.actuator, el.P_env_c)

            # Anode pressure loop:
                ModelingToolkit.connect(Pa_sp.output,     Pa_ctrl.setpoint)
                ModelingToolkit.connect(el.P_a,           Pa_ctrl.measurement)
                ModelingToolkit.connect(Pa_ctrl.actuator, el.P_env_a)

    ]

    # Create and simplify system:
        @parameters t
        @named model = ODESystem(connections, t, systems = [
            el, Pw, F_in, x_in_H2, x_in_O2, x_in_H2O,   
            T_sp, T_ctrl, 
            Pc_sp, Pc_ctrl, 
            Pa_sp, Pa_ctrl
        ])

        model_sim    = structural_simplify(model) # Simplify the system to reduce the number of equations and variables and reduce it to a index 1 DAE or ODE system

    # Solve:
        if init_state === nothing
            init_ss = deserialize(joinpath(@__DIR__, "init_ss_T333_60.jls")) # Gather previus stored solution as initial conditions and guesses
        else
            init_ss = init_state
        end


        prob = ODEProblem(model_sim, init_ss, 1000) # Define the problem to be solved, simulation time = 120 s 
        sol  = solve(prob, Rosenbrock23(), reltol = 1e-6, abstol = 1e-6) # Solve problem using Rosenbrock

        println("Solver status: ", sol.retcode)
        println("Final simulation time: ", sol.t[end])

        #saving the data for the next step, so that we can use it as initial conditions for the next step

        final_state = Dict(ModelingToolkit.states(model_sim) .=> sol.u[end] )

       if make_plot
        # Stack voltage
             p1 = plot(sol.t, sol[el.V.u],
                     title = "Cell voltage", xlabel = "Time [s]", ylabel = "Voltage [V]",
                     label = "V", lw = 2, color = :blue)
         # Stack cathode/hydrogen outflow
              p2 = plot(sol.t, sol[el.F_out_c.u],
                     title = "Cathode outlet flow", xlabel = "Time [s]", ylabel = "Molar flow [mol/s]",
                     label = "F_out_c", lw = 2, color = :green)
         # Stack cathode/hydrogen pressure
              p3 = plot(sol.t, sol[el.P_c.u]*1e-5,
                     title = "Cathode pressure", xlabel = "Time [s]", ylabel = "Pressure [bar]",
                     label = "P_c", lw = 2, color = :red)
         # Stack anode/oxygen pressure
              p4 = plot(sol.t, sol[el.P_a.u]*1e-5,
                     title = "Anode pressure", xlabel = "Time [s]", ylabel = "Pressure [bar]",
                     label = "P_a", lw = 2, color = :orange)
         # Stack temperature
              p5 = plot(sol.t, sol[el.T.u],
                     title = "Temperature", xlabel = "Time [s]", ylabel = "Temperature [K]",
                     label = "T", lw = 2, color = :purple)
         # Stack HTO in anode copartment
             # 
             p6 = plot(sol.t, sol[el.HTO_a.u],
                     title = "Safety - HTO anode", xlabel = "Time [s]", ylabel = "H2/O2 ratio [vol%]",
                     label = "HTO_a", lw = 2, color = :black)
         # Combine into one figure and show plot: 
              display(plot(p1, p2, p3, p4, p5, p6, layout = (3, 2), size = (1200, 900)))
       # --- Polarization curve (voltage vs. current density) ---
           i_exp = [0.12, 0.33, 0.5,  0.75, 1.0,  1.25, 1.5,  1.92] # [A/cm^2]
           V_exp = [1.6,  1.66, 1.71, 1.77, 1.83, 1.89, 1.95, 2.04] # [V]
           p7 = plot(sol[el.i_el]*1e-4, sol[el.V.u],
                   title = "Polarization curve", xlabel = "Current density [A/m^2]", ylabel = "Voltage [V]",
                   label = "Polarization curve", lw = 2, color = :blue)
           p7 = scatter!(i_exp, V_exp, label = "Experimental data", color = :red, marker = :circle, ms = 4)
           
           display(plot(p7, layout = (1, 1), size = (500, 300)))
       # --- HTO (HTO vs current density) ---
        #    1 bar
           i_exp = [0.25, 0.5, 1.0, 1.5, 2.0] # [A/cm^2]
           V_exp = [0.65, 0.6, 0.56, 0.55, 0.54] # [vol%]
         # 30 bar 
            i_exp = [0.25, 0.5, 1.0, 1.5, 2.0] # [A/cm^2]
            V_exp = [7.23, 4.17, 2.41, 1.76, 1.5] # [vol%]

            p8_0 = plot(sol[el.i_el]*1e-4, sol[el.HTO_a.u],
                    title = "Safety - HTO anode", xlabel = "Current density [A/cm^2]", ylabel = "H2 in O2 [vol%]",
                    label = "HTO_a", lw = 2, color = :black,
                    ylims = (0, 8))
            p8_0 = scatter!(i_exp, V_exp, label = "Experimental data", color = :red, marker = :circle, ms = 4)
            # savefig(p8_0, "/home/annkar/NewCode/HTO_plot_CurrentDensity.png")

            p8_01 = plot(sol[el.Pw.u]*1e-4, sol[el.HTO_a.u],
                    title = "Safety - HTO anode", xlabel = "Power [kW]", ylabel = "H2 in O2 [vol%]",
                    label = "HTO_a", lw = 2, color = :black,
                    ylims = (0, 8))


            
            # savefig(p8_01, "/home/annkar/NewCode/HTO_plot_Power.png")


            p8_02 = plot(sol[el.T.u], sol[el.HTO_a.u],
                    title = "Safety - HTO anode", xlabel = "Temperature [K]", ylabel = "H2 in O2 [vol%]",
                    label = "HTO_a", lw = 2, color = :black,
                    ylims = (0, 8))
            # savefig(p8_02, "/home/annkar/NewCode/HTO_plot_Temperature.png")
            
            p8_03 = plot(sol[el.P_a.u]*1e-5, sol[el.HTO_a.u],
                    title = "Safety - HTO anode", xlabel = "Anode Pressure [bar]", ylabel = "H2 in O2 [vol%]",
                    label = "HTO_a", lw = 2, color = :black,
                    ylims = (0, 8))
            # savefig(p8_03, "/home/annkar/NewCode/HTO_plot_PAnode.png")

            after_25s = sol.t .>= 250

            p8_04 = plot(sol[el.P_c.u]*1e-5, sol[el.HTO_a.u],
                title = "Safety - HTO anode", xlabel = "Cathode Pressure [bar]", ylabel = "H2 in O2 [vol%]",
                label = "HTO_a", lw = 2, color = :black,
                ylims = (0, 8))
            # savefig(p8_04, "/home/annkar/NewCode/HTO_plot_PCathode.png")

            p8_05 = plot(sol.t, sol[el.HTO_a.u], xlabel = "Time [s]", ylabel = "HTO [vol%]", label = "HTO_a")
            # savefig(p8_05, "/home/annkar/NewCode/HTO_vs_Time.png")

            display(plot(p8_05, layout = (1, 1), size = (500, 300)))

    # Loop diagnostics: controlled variable + setpoint, actuator + bounds, per loop

        # --- Temperature loop ---
            d1 = plot(sol.t, sol[el.T.u], label="T", color=:purple, lw=2,
                    title="T loop", xlabel="Time [s]", ylabel="T [K]")
            hline!(d1, [333.15], label="setpoint", color=:purple, ls=:dot, lw=1)
            d1t = twinx(d1)
            plot!(d1t, sol.t, sol[el.T_in.u], label="T_in", color=:gray, ls=:dash, ylabel="T_in [K]")

        # --- Cathode pressure loop ---
            d2 = plot(sol.t, sol[el.P_c.u]*1e-5, label="P_c", color=:red, lw=2,
                    title="P_c loop", xlabel="Time [s]", ylabel="P_c [bar]")
            hline!(d2, [P_cathode_bar], label="setpoint", color=:red, ls=:dot, lw=1)
            d2t = twinx(d2)
            plot!(d2t, sol.t, sol[el.P_env_c.u], label="z_c", color=:gray, ls=:dash, ylabel="z_c [-]")
            # savefig(d2, "/home/annkar/NewCode/Pc_loop.png")

        # --- Anode pressure loop ---
            d3 = plot(sol.t, sol[el.P_a.u]*1e-5, label="P_a", color=:orange, lw=2,
                    title="P_a loop", xlabel="Time [s]", ylabel="P_a [bar]")
            hline!(d3, [1.0], label="setpoint", color=:orange, ls=:dot, lw=1)
            d3t = twinx(d3)
            plot!(d3t, sol.t, sol[el.P_env_a.u], label="z_a", color=:gray, ls=:dash, ylabel="z_a [-]")

        display(plot(d1, layout=(1,1), size=(1000, 1000)))
       end 

           
     # Plot results:
         
    # Save steady-state solution for future use:
        if save_ss
            init_dict = Dict(ModelingToolkit.states(model_sim) .=> sol[end])
            serialize(joinpath(@__DIR__, "init_ss_T343.jls"), init_dict)
        end

     return sol, el, final_state
end

using Serialization
using CSV
using DataFrames

function create_initial_conditions()
     Temp_K = collect(353.15:-0.1:290.15) 
     Pressure_bar = collect(0.9:0.001:3.2)
    init_state = nothing  # Reset at the start of each power sweep
     for (i, T) in enumerate(Pressure_bar)
         # println("Running cathode pressure = $(p_bar) bar")
         println("Running pressure = $(T) bar")
         sol, el, final_state = test_Electrolyzer(; T_stack_K=333.15, P_cathode_bar=30, P_anode_bar=T, power=15000.0, init_state=init_state
             )
        init_state = final_state

         
        T = sol[el.T.u][end]

        Pa = sol[el.P_a.u][end]*1e-5

    

        if isapprox(Pa, 1.09, atol=0.05)
            serialize(joinpath(@__DIR__, "init_ss_Pa109_15.jls"), init_state)

        elseif isapprox(Pa, 3.09, atol=0.05)
            serialize(joinpath(@__DIR__, "init_ss_Pa309_15.jls"), init_state)
        #  elseif isapprox(T, 293.15, atol=0.05)
        #      serialize(joinpath(@__DIR__, "init_ss_T293_15.jls"), init_state)
        # elseif isapprox(T, 323.15, atol=0.05)
            # serialize(joinpath(@__DIR__, "init_ss_T323_10.jls"), init_state)
        # elseif isapprox(T, 333.15, atol=0.05)
        #     serialize(joinpath(@__DIR__, "init_ss_T333_10.jls"), init_state)
        # elseif isapprox(T, 343.15, atol=0.05)
            # serialize(joinpath(@__DIR__, "init_ss_T343_10.jls"), init_state)   
        end

        # if isapprox(Pa, 1.09, atol=0.01)
            # serialize(joinpath(@__DIR__, "init_ss_Pa109.jls"), init_state)

        # elseif isapprox(Pa, 2.0, atol=0.05)
        #     serialize(joinpath(@__DIR__, "init_ss_Pa2.jls"), init_state)
        # elseif isapprox(Pa, 3.0, atol=0.05)
        #     serialize(joinpath(@__DIR__, "init_ss_Pa3.jls"), init_state)
        # elseif isapprox(Pa, 4.0, atol=0.05)
        #     serialize(joinpath(@__DIR__, "init_ss_Pa4.jls"), init_state)
        # end

     end

    return init_state

end


function generate_surrogate_data()

    # Definer effektverdiene
    Power_values = collect(6000.0:6000.0:60000.0)

    # Konstante betingelser
    T_stack = 333.15
    P_anode = 1.0 
    P_cathode = 30.0

    # Les inn initialbetingelser
    initial_state = deserialize(
        joinpath(@__DIR__, "init_ss_Pa1_60.jls")
    )

    # Opprett DataFrame
    data = DataFrame(
        Power = Float64[],
        Temperature = Float64[],
        P_anode = Float64[],
        P_cathode = Float64[],
        HTO = Float64[]
    )

    # Kjør simuleringer
    for P in Power_values

        println("Simulating power = $P W")

        sol, el, final_state = test_Electrolyzer(;
            T_stack_K = T_stack,
            P_cathode_bar = P_cathode,
            P_anode_bar = P_anode,
            power = P,
            init_state = initial_state
        )

        # Hent HTO fra modellen
        # ERSTATT med riktig variabel fra din modell
        HTO = sol[el.HTO_a.u][end]
        
        println("HTO at power $P W: $HTO")

        # Lagre resultat
        push!(data, (P, T_stack,P_anode, P_cathode, HTO
        ))

        # Oppdater initialbetingelser
        initial_state = final_state

    end

    # Lagre resultatene
    # CSV.write(joinpath(@__DIR__, "surrogate_data.csv"), data)

    println("Simulation data saved!")

    return data

end


 

###===============================================================================
#Plot the pressure and temperature against the HTO to see if there is a correlation between them and the HTO.
function sweep_HTO(;save_excel=false)
    
    pressure_bar = collect(29.0:-1.0:1.0)
    Power=collect(60000.0:500:60000.0)
    Temp_K = collect(353.15:-2:283.15) #plot this again, begin from high to low
    
    Power_i = [6000.0, 15000.0, 30000.0, 60000.0]


    
    power_results = Float64[]
    temperature_setpoints = Float64[]
    temperature_results = Float64[]
    cathode_pressures = Float64[]
    cathode_pressure_setpoints = Float64[]
    anode_pressures = Float64[]
    anode_pressure_setpoints = Float64[]
    hto_results = Float64[]

    #steady state values
    # T_ss = zeros(length(Temp_K), length(Power_i))
    
  
    # Power_ss = zeros(length(Power))
    # P_anode_ss = zeros(length(pressure_bar), length(Power_i))
    P_cathode_ss = zeros(length(pressure_bar), length(Power_i))

    HTO_ss = zeros(length(pressure_bar), length(Power_i))



    # A new folder for each sweep preserves previous results.
    output_dir = joinpath(@__DIR__, "results", "sweep_HTO_" * Dates.format(now(), "yyyymmdd_HHMMSS_sss"))

    excel_path = joinpath(output_dir, "sweep_HTO_PressureAnode_results.xlsx")

    if save_excel
        mkpath(output_dir)
    end



    for (j, power) in  enumerate(Power_i)
        println("Running power = $(power) W")
        init_state = nothing  # Reset at the start of each power sweep

        for (i, P) in enumerate(pressure_bar)
            # println("Running cathode pressure = $(p_bar) bar")
            println("Running pressure = $(P) bar")

            sol, el, final_state = test_Electrolyzer(false, false; T_stack_K=353.15, P_cathode_bar=P, P_anode_bar=1, power=power, init_state=init_state
            )

            # Last solver result is taken as steady state
            P_cathode_ss[i, j] = sol[el.P_c.u][end] * 1e-5
            # P_anode_ss[i, j] = sol[el.P_a.u][end] * 1e-5
            HTO_ss[i, j] = sol[el.HTO_a.u][end]
            # T_ss[i, j] = sol[el.T.u][end]

            #These needs to be changed when we change what is constant, T, P_c or P_a
            push!(power_results, power)
            push!(temperature_setpoints, 353.15)
            push!(temperature_results, sol[el.T.u][end])
            push!(cathode_pressure_setpoints, pressure_bar[i])
            push!(anode_pressure_setpoints,1.0)
            push!(cathode_pressures, sol[el.P_c.u][end] * 1e-5)
            push!(anode_pressures, sol[el.P_a.u][end] * 1e-5)
            push!(hto_results, HTO_ss[i, j])

            init_state = final_state

            if save_excel
                # Update the workbook after every simulation.S
                # All previously completed runs remain in the columns.
                XLSX.writetable(
                    excel_path,
                    [   power_results, temperature_setpoints, temperature_results, cathode_pressure_setpoints,anode_pressure_setpoints,     cathode_pressures,  
                        anode_pressures,  hto_results
                    ],
                    [
                        "Power_W","Temperature_setpoint_K","Final_temperature_K", "Final_cathode_pressure_setpoint_bar", "Final_anode_pressure_setpoint_bar",
                        "Final_cathode_pressure_bar",  "Final_anode_pressure_bar","Final_HTO_vol_percent"
                    ];
                    sheetname = "Sweep_results",
                    overwrite = true
                )
            end
        end

        # p = plot(Temp_K, HTO_ss;
        #     marker = :circle,
        #     linewidth = 2,
        #     xlabel = "Cathode pressure [bar]",
        #     ylabel = "HTO [vol%]",
        #     title = "Steady-state HTO vs cathode pressure at power = $(power) W",
        #     label = "HTO steady state"
        # )

        # display(p)

   
        
    end

    # Sort temperature from low to high for the surface
    
    # order = sortperm(Temp_K)
    
    # p3d = surface( Power_i ./ 1000,  [order],  HTO_ss[order, :];
    #     xlabel = "Power [kW]", ylabel = "Pressure Cathode [bar]",  zlabel = "HTO [vol%]",
    #     title = "HTO versus Pressure cathode and Power",

    #     size = (1200, 900),   guidefontsize = 10,   tickfontsize = 20,    titlefontsize = 30,

    #     xticks = Power_i ./ 1000,    
    #     yticks = [1, 5, 10, 15, 20, 25, 29],
    #     # yticks = 290:60:350,

    #     color = :plasma,  colorbar = true,  colorbar_title = "HTO [vol%]",  camera = (45, 30),
    #     legend = false
    # )
    
    # p3d = surface( Power_i ./ 1000, pressure_bar[order],  HTO_ss[order, :];
    #     xlabel = "Power [kW]", ylabel = "Pressure Cathode [bar]",  zlabel = "HTO [vol%]",
    #     title = "HTO versus Pressure cathode and Power",

    #     size = (1200, 900),   guidefontsize = 10,   tickfontsize = 20,    titlefontsize = 30,

    #     xticks = Power_i ./ 1000,    
    #     yticks = [1, 5, 10, 15, 20, 25, 29],
    #     # yticks = 290:60:350,

    #     color = :plasma,  colorbar = true,  colorbar_title = "HTO [vol%]",  camera = (45, 30),
    #     legend = false
    # )

    order  = sortperm(pressure_bar)   # ascending y-axis; Temp_K is listed low→high
    porder = sortperm(Power_i)   # ascending x-axis; Power_i is listed high→low     
    
    plotlyjs()   # interactive backend    
    p3d = surface( Power_i[porder] ./ 1000, pressure_bar[order],  HTO_ss[order, porder];
    xlabel = "Power [kW]", ylabel = "Cathode Pressure [bar]",  zlabel = "HTO [vol%]",
    title = "HTO versus Cathode pressure and Power",

    size = (800, 600),
    guidefontsize = 16,
    tickfontsize = 10,
    colorbar_titlefontsize = 14,
    colorbar_tickfontsize = 12,
    titlefontsize = 18,

    xticks = Power_i ./ 1000,    
    yticks = [1, 5, 10, 15, 20, 25, 29],
    # ytick = [293.15, 303.15, 318.15, 333.15, 353.15],

    color = :plasma,  colorbar = true,  colorbar_title = "HTO [vol%]",  camera = (45, 30),
    legend = false    
    )

    

    savefig(p3d, joinpath(@__DIR__, "HTO_CathodePressure_Power_3D.html"))    
    gr()   # switch back so later plots use the normal backend
    # savefig(p3d, joinpath(@__DIR__, "HTO_PressureCathode_Power_3D.png"))
    # savefig(p3d, joinpath(@__DIR__, "HTO_TemperatureCathode_Power_3D.png"))
    # savefig(p3d, joinpath(@__DIR__, "HTO_AnodeCathode_Power_3D.png"))

    if save_excel
        # savefig(p, joinpath(output_dir, "HTO_vs_Temperature.png"))
        println("Results saved to: ", excel_path)
    end
    return temperature_results, HTO_ss, Power_i
end

###===============================================================================
#Create a G matrix explaining how sensitive the system is for different changes


function run_Electrolyzer(;Pw=30000, 
    F_in=35.09, 
    z_c=171.55768399498805,
    z_a = 244.80834801269452,
    P_env_c =26.0,    
    P_env_a=0.9,
    T=343.15,
    init_state=nothing)
    
    
    # Create electrolyzer model:
    @named el       = Electrolyzer(par_El)

    #inputs
    @named Pw       = Blocks.Constant(k = Pw) # [w] Input power
    @named F_in     = Blocks.Constant(k = F_in) # [mol/s] Molar inlet flow
    @named z_c      = Blocks.Constant(k = z_c) # [-] Valve opening (unscaled) 
    @named z_a      = Blocks.Constant(k = z_a) # [-] Valve opening (unscaled)
    @named x_in_H2  = Blocks.Constant(k = 0) # [-] Composition in (hydrogen)
    @named x_in_O2  = Blocks.Constant(k = 0) # [-] Composition in (oxygen)
    @named x_in_H2O = Blocks.Constant(k = 1) # [-] Composition in (water)
    @named P_env_c  = Blocks.Constant(k = P_env_c*1e5) # [Pa] Pressure cathode seperator
    @named P_env_a  = Blocks.Constant(k = P_env_a*1e5) # [Pa] Pressure anode seperator
    @named T_in     = Blocks.Constant(k=T) #[K]


    connections = [
        connect(Pw.output, el.Pw)
        connect(F_in.output, el.F_in)

        connect(z_c.output, el.z_c)
        connect(z_a.output, el.z_a)

        connect(x_in_H2.output, el.x_in_H2)
        connect(x_in_O2.output, el.x_in_O2)
        connect(x_in_H2O.output, el.x_in_H2O)

        connect(P_env_c.output, el.P_env_c)
        connect(P_env_a.output, el.P_env_a)

        connect(T_in.output, el.T_in)
    ]

        # Create and simplify system:
        @parameters t
        @named model = ODESystem(connections, t, systems = [
            el, Pw, F_in,
            z_c, z_a,
            x_in_H2, x_in_O2, x_in_H2O,
            P_env_c, P_env_a,
            T_in
        ])

        model_sim    = structural_simplify(model) # Simplify the system to reduce the number of equations and variables and reduce it to a index 1 DAE or ODE system

        # Solve:
        if init_state === nothing
            init_ss = deserialize(joinpath(@__DIR__, "init_ss_26_bar.jls")) # Gather previus stored solution as initial conditions and guesses
        else
            init_ss = init_state
        end


        prob = ODEProblem(model_sim, init_ss, 1000) # Define the problem to be solved, simulation time = 120 s 
        sol  = solve(prob, Rosenbrock23(), reltol = 1e-6, abstol = 1e-6) # Solve problem using Rosenbrock

        #saving the data for the next step, so that we can use it as initial conditions for the next step

        final_state = Dict(ModelingToolkit.states(model_sim) .=> sol.u[end]
                    )

        #Here i need to solve wuth the init_ss_30_bar the first step and then 



        p1 = plot(sol.t, sol[el.HTO_a.u], xlabel = "Time [s]", ylabel = "HTO [vol%]", label = "HTO_a")
        p2 = plot(sol.t, sol[el.P_a.u]*1e-5, xlabel = "Time [s]", ylabel = "P anode [bar]", label = "Pressure anode")
        p3 = plot(sol.t, sol[el.P_c.u]*1e-5, xlabel = "Time [s]", ylabel = "P cathode [bar]", label = "Pressure cathode")
        p4 = plot(sol.t, sol[el.T.u], xlabel = "Time [s]", ylabel = "T [K]", label = "Temperature")
        
        display(plot(p1, p2, p3, p4, layout=(4,1), size=(500, 500)))   

        HTO = sol[el.HTO_a.u][end]
        Pc = sol[el.P_c.u][end]
        Pa = sol[el.P_a.u][end]
        T = sol[el.T.u][end]

        @show HTO Pc Pa T


        return HTO, Pc, Pa, T, final_state
end

function G_matrix(Pw0 = 60000.0,     # W
    P_env_c0 = 27.0, # bar
    P_env_a0 = 1.0,    # bar
    T_in0 = 333.15,     # K
    F_in0 = 35.09
    )
###############################################################################
    # Nominal input values.
 ###############################################################################  
    
 state0 = deserialize(joinpath(@__DIR__, "init_ss_T333_60.jls"))

    sol, el, final_state = test_Electrolyzer(;
        P_cathode_bar=P_env_c0, P_anode_bar=P_env_a0, T_stack_K=T_in0, power=Pw0, init_state=state0

    )

    HTO0 = sol[el.HTO_a.u][end]
    Pc0 = sol[el.P_c.u][end]*1e-5
    Pa0 = sol[el.P_a.u][end]*1e-5
    T0 = sol[el.T.u][end]


    @show HTO0 Pc0 Pa0 T0

    #SAFETY: cheking if the HTO0 is zero such that it will not divide on zero


    if iszero(HTO0)
        error("Cannot normalise gains because nominal HTO is zero.")
    end
###############################################################################
# Change in power
###############################################################################
    # delta_Pw = -3000
    # sol_Pw, el_Pw, final_state_Pw =test_Electrolyzer(;
    #  P_cathode_bar=P_env_c0, P_anode_bar=P_env_a0, T_stack_K=T_in0, power=Pw0+delta_Pw, init_state=state0
 
    # )

###############################################################################
    # Normalizing the gain in order to compare them to eachother
##############################################################################

    # HTO_Pw = sol_Pw[el_Pw.HTO_a.u][end]
    # G_Pw = ((HTO_Pw-HTO0)/HTO0)/(delta_Pw/Pw0)

###############################################################################
# Change in pressure at cathode
###############################################################################
    delta_Pc = 2.6
    sol_Pc, el_Pc, final_state_Pc =test_Electrolyzer(;
     P_cathode_bar=P_env_c0 + 2.6 , P_anode_bar=P_env_a0, T_stack_K=T_in0, power=Pw0, init_state=state0
 
    )

    HTO_Pc = sol_Pc[el_Pc.HTO_a.u][end]
    Pc_Pc = sol_Pc[el_Pc.P_c.u][end]*1e-5
    Pa_Pc = sol_Pc[el_Pc.P_a.u][end]*1e-5
 
    G_Pc = (HTO_Pc-HTO0)/(delta_Pc)

    s_Pc = G_Pc * (P_env_c0 / HTO0)

    @show G_Pc s_Pc

###############################################################################
#change in pressure at anode
##############################################################################
    
    # delta_Pa = 0.09
    # state_Pa = deserialize(joinpath(@__DIR__, "init_ss_Pa309_15.jls"))
    # sol_Pa, el_Pa, final_state_Pa =test_Electrolyzer(;
    #  P_cathode_bar=P_env_c0, P_anode_bar=P_env_a0 + 0.09, T_stack_K=T_in0, power=Pw0, init_state=state_Pa
 
    # )



    # HTO_Pa = sol_Pa[el_Pa.HTO_a.u][end]
    # Pc_Pa = sol_Pa[el_Pa.P_c.u][end]*1e-5
    # Pa_Pa = sol_Pa[el_Pa.P_a.u][end]*1e-5
    # G_Pa = (HTO_Pa-HTO0)/(delta_Pa)

    # s_Pa = G_Pa * (P_env_a0 / HTO0)

    # @show G_Pa s_Pa Pa_Pa

###############################################################################   
    #change in temperature
###############################################################################

    # state_T = deserialize(joinpath(@__DIR__, "init_ss_T293_60.jls"))
    # delta_T = 10
    # sol_T, el_T, final_state_T =test_Electrolyzer(;
    #  P_cathode_bar=P_env_c0, P_anode_bar=P_env_a0, T_stack_K=T_in0 + 10.0, power=Pw0, init_state=state_T
    # )


    # HTO_T = sol_T[el_T.HTO_a.u][end]
    # Pc_T = sol_T[el_T.P_c.u][end]*1e-5
    # Pa_T = sol_T[el_T.P_a.u][end]*1e-5
    # T_T = sol_T[el_T.T.u][end]

    # # Unnormalized temperature gain
    # G_T = (HTO_T - HTO0) / delta_T

    # # Normalized temperature sensitivity
    # s_T = G_T * (T_in0 / HTO0)

    # @show HTO_T Pc_T Pa_T T_T G_T s_T


    println("\nG matrix:")
    println("\nNominal HTO: $HTO0")

    return s_Pa
end




#### ==========================================================================================================================
# For convinience I have included a PI-controller from PI.jl here: 

using ModelingToolkit, ModelingToolkitStandardLibrary.Blocks

"""
Controller(; name, use_controller, Kc, Ti, u_min, u_max, int_x0,
                step_offset, step_height=0.0, step_start=100.0)

Switchable single-loop input source, chosen at build time.

- `use_controller = true`  : a `LimPID` drives `actuator` from `setpoint` and
`measurement` (closed loop). Tuning via `Kc`, `Ti`, `u_min`, `u_max`, `int_x0`.
- `use_controller = false` : a `Step` drives `actuator` (open loop). The
`measurement` port is present but unused. Signal set by `step_offset`,
`step_height`, `step_start`.

Ports (both modes):
    setpoint    : RealInput   reference / target value
    measurement : RealInput   process measurement (used only in closed loop)
    actuator    : RealOutput  signal to the actuator (e.g. valve opening)
"""
function Controller(; name, use_controller::Bool,
                      Kc = 1.0, Ti = 1.0, u_min = 0.0, u_max = 1e3, int_x0 = 0.0,
                      step_offset = 0.0, step_height = 0.0, step_start = 100.0)

    @named setpoint    = RealInput()
    @named measurement = RealInput()
    @named actuator    = RealOutput()

    if use_controller
        @named pid = LimPID(; k = Kc, Ti = Ti, Td = false,
                              u_max = u_max, u_min = u_min, int__x = int_x0)
        eqs = [
            connect(setpoint,       pid.reference)
            connect(measurement,    pid.measurement)
            connect(pid.ctr_output, actuator)
        ]
        systems = [setpoint, measurement, actuator, pid]
    else
        @named step = Blocks.Step(height = step_height, offset = step_offset,
                                   start_time = step_start, duration = Inf, smooth = true)
        eqs = [ connect(step.output, actuator) ]
        systems = [setpoint, measurement, actuator, step]
    end

    # Use the caller-supplied name, NOT @named (which hardcodes :sys):
    return ODESystem(eqs, ModelingToolkit.get_iv(setpoint);
                     systems = systems, name = name)
end






export Electrolyzer, test_Electrolyzer, sweep_HTO, run_Electrolyzer, G_matrix, create_initial_conditions, iterate, generate_surrogate_data