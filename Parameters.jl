par = Dict{Symbol,Float64}(                 # General:
     :F         => 96485.336521             # s A/mol        Farady constant
    ,:R         => 8.314462175              # J/(K mol)      Universal gas constant
    ,:g         => 9.81                     # m/s^2          Gravitational acceleration
)

par_El = Dict{Symbol,Float64}(              # Electrolyzer:
    # Dimesions stack:
     :Nc        => 65                       # -              Number of stacks
    ,:A_cell    => 214*1e-4                 # m^2            Area of membrane
    ,:A_stack   => 214*65*1e-4              # m^2            Area of membrane per stack
    ,:delta_mem => 0.000322                 # m              Tichnes of membrane
    ,:delta_c   => 1e-3                     # m              Thichnes cathode
    ,:delta_a   => 1e-3                     # m              Thichnes anode
    ,:V_c       => 0.0006955                # m^3            Volume cathode
    ,:V_a       => 0.0006955                # m^3            Volume anode

    # Electrochemical parameters:
    ,:z         => 2                        # -              Number of electrones
    ,:lambda    => 22                       # -              Water content membrane
    ,:R_el      => 3.193346538545269e-7     # Ohm            Electrical resistance
    ,:alpha_a   => 0.5                      # -              Charge transfer coefficient anode
    ,:i_a0_ref  => 0.014013822460373167     # A/m^2          Reference exhange current density anode
    ,:E_a       => 76*1e3                   # J/mol          Activation energy anode

    # Mass transfer parameters:
    ,:eta_F     => 0.99                     # -              Carrnot efficiency

    ,:a         => 3.640e-4                 # -              
    ,:b         => 2.334                    # -
    ,:epsilon   => 0.3                      # -
    ,:epsilon_p => 0.11                     # -
    ,:alpha     => 0.785                    # -

    ,:Dw        => 1.28e-10                 # m^2/s           Diffusion coefficient of water in memebrane

    ,:nd        => 2.5036                   #                 Electrosmotic drag coefficient (usully not constant)
    ,:K_Darcy   => 1.58*10^(-18)            # m^2             Membrane permability
    ,:e_H2      => 5.31*1e-14               # mol /(m s Pa)
    ,:A_H2      => 80                       # Pa/(A/m^2)      Pressure drop coefficient hydrogen
    ,:A_O2      => 94                       # Pa/(A/m^2)      Pressure drop coefficient oxygen
    ,:e_O2      => 2.26*1e-14               # mol /(m s Pa)
    ,:CH2O      => 2.3*1e4                  # mol/(m^3)
    ,:Su_H2     => 0.72*1e-6                # mol/(m^3*Pa)    Solubility hydrogen in membrane
    ,:Su_O2     => 0.8*1e-6                 # mol/(m^3*Pa)    Solubility oxygen in membrane

    #diffusion
    ,:E_e_H2    =>  3.32*1e-20              # J              Activation energy hydrogen diffusion 
    ,:E_e_O2    =>  2.81*1e-20              # J              Activation energy oxygen diffusion
    # Div:
    ,:h_a       => 0.013657283227322643     # (Calibrate to fit valve opening 1)
    ,:h_c       => 0.00014279660440722695   # (Calibrate to fit valve opening 1)
)

par_H2 = Dict{Symbol,Float64}(              # Hydrogen (H2):
     :M_H2      => 2.016*1e-3               # kg/mol          Molar mass hydrogen
    ,:rho_H2    => 0.1                      # g/cm^3          Density hydrogen (liquid)
    ,:Cp_H2     => 29                       # J/(K mol)       Heat capacity hydrogen
    ,:DH_H2     => 0                        # J/mol           Enthalpy of formation hydrogen
    ,:DG_H2     => 0                        # J/mol           Gibbs free energy hydrogen
    ,:S_H2      => 131                      # J/(K mol)       Entropy hydrogen
    ,:P_c_H2    => 33.3                     # Pa              Critical pressure hydrogen
    ,:T_c_H2    => 12.8                     # K               Critical temperature hydrogen

)

par_O2 = Dict{Symbol,Float64}(              # Oxygen (O2):
     :M_O2      => 32.00*1e-3               # kg/mol          Molar mass oxygen
    ,:rho_O2    => 1.1                      # g/cm^3          Density oxygen (liquid)
    ,:Cp_O2     => 29                       # J/(K mol)       Heat capacity oxygen
    ,:DH_O2     => 0                        # J/mol           Enthalpy of formation oxygen
    ,:DG_O2     => 0                        # J/mol           Gibbs free energy oxygen
    ,:S_O2      => 205                      # J/(K mol)       Entropy oxygen
    ,:P_c_O2    => 154.4                    # Pa              Critical pressure oxygen
    ,:T_c_O2    => 49.7                     # K               Critical temperature oxygen
)

par_H2O = Dict{Symbol,Float64}(             # Water (H2O):
     :M_H2O     => 18.0*1e-3                # kg/mol          Molar mass water
    ,:mu_H2O    => 1.1e-2                   # Pa*s            Dynamic viscosity water
    ,:rho_H2O   => 1000                     # kg/m3           Density water
    ,:Cp_H2O    => 75                       # J/(K mol)       Heat capacity water
    ,:DH_H2O    => -286*1e3                 # J/mol           Enthalpy of formation water
    ,:DG_H2O    => -237*1e3                 # J/mol           Gibbs free energy water
    ,:S_H2O     => 70                       # J/(K mol)       Entropy water
    ,:P_c_H2O   => 647.3                    # Pa              Critical pressure water
    ,:T_c_H2O   => 218.3                    # K               Critical temperature water
    ,:Cp_H2O_g  => 34                       # J/(K mol)       Heat capacity water gas
)
