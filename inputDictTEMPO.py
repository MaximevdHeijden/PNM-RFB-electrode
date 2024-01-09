input_dict = {
    # Dimensions (x-dimension = 0, y-dimension = 1, z-dimension = 2)
    # Specifying the dimensions is required for assigning the right boundary
    # conditions to the right labeled pores. 
    "height_dimension": 1,              # Dimension corresponding with the height of the electrode (sides of the electrode)
    "length_dimension": 2,              # Dimension corresponding with the length of the electrode (flow dimension)
    "width_dimension": 0,               # Dimension corresponding with the thickness of the electrode (current collector <-> membrane)
    
    # Phases
    "anolyte_conductivity": 3.1 ,       # Anolyte electrical conductivity [S m-1]   
    "anolyte_density": 852,             # Anolyte density [kg m-3]
    "anolyte_viscosity": 3.43e-4,       # Anolyte viscosity [Pa s] 
    "anolyte_inlet_velocity": 2e-2,     # Anolyte inlet velocity [m s-1] (used if Neumann inlet boundary condition is used in the StokesFlow algorithm)
    "anode_pressure_drop": 1e3,         # Pressure drop over the anode [Pa] (used if Dirichlet inlet boundary condition is used in the StokesFlow algorithm)
    
    "catholyte_conductivity": 3.1 ,     # Anolyte electrical conductivity [S m-1] 
    "catholyte_density": 852,           # Anolyte density [kg m-3]
    "catholyte_viscosity": 3.43e-4,     # Anolyte viscosity [Pa s] 
    "catholyte_inlet_velocity": 2e-2,   # Anolyte inlet velocity [m s-1] (used if Neumann inlet boundary condition is used in the StokesFlow algorithm)
    "cathode_pressure_drop": 1e3,       # Pressure drop over the cathode [Pa] (used if Dirichlet inlet boundary condition is used in the StokesFlow algorithm)
    
    # Active species / kinetic parameters
    "D_a": 1.3e-9,                      # Anolyte active species diffusivity (Fe2+) [m2 s-1]
    "conc_in_a": 100,                   # Anolyte active species inlet concentration [mol m-3]
    "conc_ref_a": 100,                  #  Anolyte active species reference concentration [mol m-3] (concentration at which j0 was measured)
    "j0_a": 1750,                       # Anolyte exchange current density [A m-2]
    "val_a": 1,                         # Amount of valence electrons participating in the half-reaction [-]
    "alpha_a_a": 0.5,                   # Anodic transfer coefficient for the half-reaction in the anode [-]
    "alpha_c_a": 0.5,                   # Cathodic transfer coefficient for the half-reaction in the anode [-]
    
    "D_c": 1.3e-9,                      # Catholyte active species diffusivity (Fe3+) [m2 s-1]
    "conc_in_c": 100,                   # Catholyte active species inlet concentration [mol m-3]
    "conc_ref_c": 100,                  # Catholyte active species reference concentration [mol m-3] (concentration at which j0 was measured)
    "j0_c": 1750,                       # Catholyte exchange current density [A m-2]
    "val_c": 1,                         # Amount of valence electrons participating in the half-reaction [-]
    "alpha_a_c": 0.5,                   # Anodic transfer coefficient for the half-reaction in the cathode [-]
    "alpha_c_c": 0.5,                   # Cathodic transfer coefficient for the half-reaction in the cathode [-]

    # Cell potential parameters
    "E_red_a": 0.77,                    # Standard reduction potential of the anodic half reaction [V]
    "E_red_c": 0.77,                    # Standard reduction potential of the cathodic half reaction [V]
    "V_step": -0.1,                     # Step change in the cell voltage [V] 
    "E_cell_final": -1.0,                 # Final value of the cell voltage range [V] # (E_cell < 0 --> Charging, E_cell > 0 --> Discharging)
    
    # Membrane properties
    "membrane_porosity": 0.58,          # Membrane porosity [m3 void/ m3 total]
    "membrane_thickness": 0.000175,     # Membrane thickness [m]
    "membrane_resistivity":  1.278E-04, # Membrane ionic resistivity [Ohm m2] 
    
    # Universal constants
    "R_g": 8.314,                       # Universal gas constant [J K-1 mol-1]
    "F": 96485.333,                     # Faradaic constant [C mol-1]
    "T": 298,                           # Operating temperature [K]
    
    # fullCellRFBSeries network parameter
    "total_electrode_length": 1.7e-2,   # Electrode length that is approximated [m] (using the ceil function, so the calculated network is slightly larger)
    
    # Numerical settings
    "rel_tol": 5e-5,                    # Relative tolerance between the solutions for two sequential iterations  
    "abs_tol": 6e-4,                    # Absolute tolerance between the estimated external current density between anodic and cathodic side
    "max_iter": 50000,                  # Maximum amount of iterations
    "omega": 0.1,                       # Damping factor
    
    # Mass transfer coefficient - fit km 
    # The Mass transfer equation is: km = C1_old/a * vel^(C2) = C1 * vel^(C2) 
    "MT_coeff_C1": 0.00016 ,            # Altered network 0.00016     Original network 0.000143        
    "MT_coeff_C2": 0.80,                # Altered network 0.82        Original network 0.80       
    }