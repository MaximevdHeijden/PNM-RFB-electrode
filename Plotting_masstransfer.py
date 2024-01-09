""" ---------Importing Packages and set Load/Read directories --------------"""
import openpnm as op
import numpy as np
import inputDictTEMPO as inputDict  
import os
import Custom_functions_km as cf
import Custom_network as cn
import Custom_functions_transport as cf_trans
import Custom_functions_phase as cf_phase
import matplotlib.pyplot as plt
import Custom_functions_pressure_fitting as cf_pres_fit  # Pressure
param = inputDict.input_dict                        # Load in the input dictionary
import scipy.optimize as optimize

# Input velocity
v_in = 1.5e-2 

# If Original_network == 0 -> We use the adapted network (i.e. downsized throat diameters).
# If Original_network == 1 -> We use the original extracted network.
# **************** Do not forget to switch the local MT correlation exponents 
#                  C1 and C2 (inputDict) before running the script.************
Original_network = 1

# Creating the normalized km plot in Origin:
# 0. Run the script at the desired input velocity
# 1. Open "bins_vel_ind" (the x-axis) in the variable explorer and copy paste to Origin
# 2. Open the NORMALIZED velocity dependent distribution (n_Vel_dep_nrom) 
#    and velocity independent distribution (n_Vel_ind_norm) and copy paste to Origin

"""----------Loading in the Network, phases and fluid fields ---------------"""
network_name = 'Freudenberg_check_CSV'
# Setting import network directory
cwd = os.getcwd()                                   # Curent Working directory

# Load in the *original* extracted network. Moreover, initialize the phase, and stokesflow
# and compute the mass transfer coefficients for the velocity dependent and independent case. 
if Original_network == 1:
    network_path = '\\input\\'                      # Folder containing network
    file_name_net_c = 'Freudenberg_check_CSV'       # Network name
    path_net_c = cwd + network_path + file_name_net_c
    file_name = path_net_c + '.pnm'
    project_cat = op.Workspace().load_project(filename = file_name)
    net_c = project_cat.network   
               
    # Checking for isolated and disconnected pores:
    health_c = op.utils.check_network_health(network = net_c)
    print("The health of the cathodic and anodic electrodes are as",
          "follows (note that they are the same mirrored network):\n")
    print(health_c); 

    # Update the face labels that are defined in the SNOW algorithm
    cf.validate_face_labels(net_c) 
    
    """-----------Define and compute network properties-------------"""
    # Define electrode dimensions [x=0, y=1, z=2]
    # NOTE: Make sure you know which dimension correlates with which face in your image.
    H_dim = param['height_dimension']               # The dimension corresponding to the height of the electrodes (sides of the electrode)
    L_dim = param['length_dimension']               # The dimension corresponding to the length of the electrodes (flow field direction)
    W_dim = param['width_dimension']                # The dimension corresponding to the width of the electrodes (thickness: current collector -> membrane)            
    
    L_c = np.amax(net_c['pore.coords'][net_c.pores(), L_dim])  # Length [m]
    H_c = np.amax(net_c['pore.coords'][net_c.pores(), H_dim])  # Height [m]
    W_c = np.amax(net_c['pore.coords'][net_c.pores(), W_dim])  # Width [m]         
    
    A_ext_c = L_c * H_c                             # Cathode area to compute external current density [m2]
    A_in_c = W_c * H_c                              # Cathode inlet area [m2]
    mem_area = A_ext_c                              # Fictitious membrane area [m2]            
    SA_network = net_c['pore.surface_area']         # Network surface area [m2]    
    Total_network_area = SA_network.sum()           # Total network surface area [m2]

    # Assign the labels 'membrane', 'current collector' 'flow inlet' and 'flow outlet'
    # to the correct boundary pores.
    cf.assign_boundary_pores(net_c, W_dim, L_dim)
    
    """-----------Define and compute phase properties---------------"""    
    # Phase constants are assigned via models, as during phase
    # regeneration the water object inherent functions would overwrite 
    # assigned phase constants.
    Conductivity_value = param['anolyte_conductivity']
    catholyte = op.phase.Water(network=net_c, name='catholyte')  
    catholyte.add_model(propname = 'pore.electrical_conductivity',              # Catholyte electrical conductivity [S m-1]
                      model = cf_phase.custom_electrical_conductivity_fit,
                      conductivity_fitted = Conductivity_value) 
    catholyte.add_model(propname = 'pore.diffusivity',                          # Catholyte active species diffusivity in electrolyte [m2 s-1]
                      model = cf_phase.custom_diffusivity,
                      parameter_script = param, prop = 'D_c')           
    catholyte.add_model(propname = 'pore.viscosity',                            # Catholyte viscosity [Pa s]
                      model = cf_phase.custom_viscosity,
                      parameter_script = param, prop = 'catholyte_viscosity')                 
    catholyte.add_model(propname = 'pore.density',                              # Catholyte density [kg m-3]  
                      model = cf_phase.custom_density,
                      parameter_script = param, prop = 'catholyte_density')
    catholyte['pore.molecular_weight'] = 0.01802
    
    #------------------------ Adding phase physics ---------------------------#
    # Custom functions (copy of those used in OpenPNM V2.2.0) are used
    # due to discrepencies w.r.t. the diameter of the throat exceeding
    # that of the pores from the network extraction.
    
    f_hyd = cf_trans.Flow_shape_factors_ball_and_stick
    catholyte.add_model(propname = 'throat.flow_shape_factors', model = f_hyd)    
    f_Hyd_cond = cf_trans.Hydraulic_conductance_Hagen_Poiseuille
    catholyte.add_model(propname = 'throat.hydraulic_conductance', model = f_Hyd_cond)
    f_poisson = cf_trans.Poisson_shape_factors_ball_and_stick
    catholyte.add_model(propname = 'throat.poisson_shape_factors', model = f_poisson)
    f_diff_cond = cf_trans.Diffusive_conductance_mixed_diffusion
    catholyte.add_model(propname = 'throat.diffusive_conductance', model = f_diff_cond)
    f_ad_dif = cf_trans.Advection_diffusion
    catholyte.add_model(propname = 'throat.ad_dif_conductance', model = f_ad_dif)
    f_elec_con = cf_trans.Electrical_conductance_series_resistors
    catholyte.add_model(propname = 'throat.electrical_conductance', model = f_elec_con)            
    catholyte.regenerate_models()   
    
    """-----------------------Stokes algorithm----------------------"""
    # Create and initialize stokes flow algorithm. 
    sf_c_odd = op.algorithms.StokesFlow(network=net_c, phase=catholyte)
     
    Q_in_c = v_in * A_in_c                      # total flow rate entering the network [m3/s]

    # Finding the inlet and outlet throats of the network.
    inlet_throats_c_odd = net_c.find_neighbor_throats(pores = net_c.pores('flow_inlet'))        
    inlet_throats_c_even = net_c.find_neighbor_throats(pores=net_c.pores('flow_outlet'))
    
    # The flow rate entering a pore has to be scaled with the area of the inlet pore.
    cross_area_c = np.pi / 4 * net_c['pore.diameter'] ** 2
    total_area_inlet_c = np.sum(cross_area_c[net_c.pores('flow_inlet')])
    total_area_outlet_c = np.sum(cross_area_c[net_c.pores('flow_outlet')])
    
    # Imposing the Neumann inlet boundary condition 
    for pore in net_c.pores('flow_inlet'):
        sf_c_odd.set_rate_BC(rates=Q_in_c * cross_area_c[pore] / total_area_inlet_c, pores=pore)    # Neumann inlet boundary condition
    
    # Dirichlet outlet boundary condition
    Pout = 0    # Pressure outlet boundary condition [Pa]
    sf_c_odd.set_value_BC(values=Pout, pores=net_c.pores('flow_outlet'))                            # Dirichlet outlet boundary condition
    
    # Run electrolyte transport algorithms
    sf_c_odd.run()
    
    # Now we iteratively solve for an inlet pressure boundary condition
    # that adheres to the set inlet flow rate. This is done via fsolve 
    # (solve for zero) of cf.inlet_pressure, which outputs the following 
    # ratio: (Q_tot - Q_desired)/Q_tot where Q is the inlet (odd net) or
    # outflow (even net) volumetric flowrate. 
    # For this function we need to supply an initial pressure
    # guess, which is the found mean pressures from the stokesflow alg.
    # However, each electrolyte phase can only hold the results of one 
    # network section. Hence now, and throughout the iterative algorithm,
    # we will update these phases results depending on the network number
    # that is iterated over.

    # Updating phase with sf_c_odd to extract initial guess 
    catholyte.update(sf_c_odd.soln)
    x0_c_odd = catholyte['pore.pressure'][net_c.pores('flow_inlet')].mean()
    
    # Delete the initial guess Neumann (flow rate) boundary condition
    sf_c_odd.clear_rate_BCs()

    def inlet_pressure(P_in, alg, Q_desired, inlet_throats, inlet_pores, net, phase, outlet_pores):
        r''' inlet_pressure is used as a target function to find the required inlet pressure
        in every pore to match the desired inlet flow rate Q_desired'''
        alg.set_value_BC(pores = net.pores('internal', mode = 'nor'), mode='remove')
        alg.set_value_BC(values = 0.0, pores=outlet_pores, mode='overwrite')    # Dirichlet boundary inlet condition
        alg.set_value_BC(values=P_in, pores=inlet_pores, mode='overwrite')      # Dirichlet boundary inlet condition
        alg.run()
        
        phase.update(alg.soln)
        
        pore_1 = net['throat.conns'][inlet_throats][:, 0]
        pore_2 = net['throat.conns'][inlet_throats][:, 1]
        
        delta_P = (phase['pore.pressure'][pore_2] - phase['pore.pressure'][pore_1])
        Q_tot = (phase['throat.hydraulic_conductance'][inlet_throats] * delta_P).sum()
    
        return (Q_tot - Q_desired)/Q_tot

    # Find the pressure at the inlet at which the total flow rate matches the desired flow rate.
    optimize.fsolve(func=inlet_pressure, x0=x0_c_odd, args=(sf_c_odd, Q_in_c, inlet_throats_c_odd, net_c.pores('flow_inlet'), net_c, catholyte, net_c.pores('flow_outlet')))
        
    # Compute and extract pore and throat velocity.
    cf_pres_fit.Throat_pore_velocity_extraction_OpenPNM(net_c, 1, sf_c_odd)
    
    # Retrieve the StokesFlow data for the odd network.
    catholyte.update(sf_c_odd.soln)
    cf_pres_fit.Update_phase_with_SF_OpenPNM(sf_c_odd, catholyte, 1)

    '''------------------ Mass transfer coefficient ----------------------------'''
    # The local mass transfer equation is: km = C1 * vel^(C2) 
    
    # Update phases with Odd networks
    cf_pres_fit.Update_phase_with_SF_OpenPNM(sf_c_odd, catholyte, 1)    
    
    # Initializing mass transfer coefficient
    catholyte['pore.km_exp'] = np.nan
    
    # Velocity-independent
    cf_pres_fit.Mass_transfer_coefficient_OpenPNM(net_c, catholyte, 1, 0, param)
    MT_OpenPNM_cat_odd_vel_independent = catholyte['pore.km_exp'].copy()
    
    cf_pres_fit.Mass_transfer_coefficient_OpenPNM(net_c, catholyte, 1, 1, param)
    MT_OpenPNM_cat_odd_vel_dependent = catholyte['pore.km_exp'].copy()

    '''---------------------- Extract fluid field data -------------------------'''
    # Obtain permeability in the flow direction
    inlet_throats_c_odd = net_c.find_neighbor_throats(pores = net_c.pores('flow_inlet'))
    Inlet_surface_pores_c_odd = net_c.find_connected_pores(throats = inlet_throats_c_odd)[:,0]
    outlet_throats_c_odd = net_c.find_neighbor_throats(pores = net_c.pores('flow_outlet'))
    Outlet_surface_pores_c_odd = net_c.find_connected_pores(throats = outlet_throats_c_odd)[:,0]
    # Extracting pressure
    mean_pressure_inlet_c_odd = sf_c_odd['pore.pressure'][Inlet_surface_pores_c_odd].mean()
    mean_pressure_outlet_c_odd = sf_c_odd['pore.pressure'][Outlet_surface_pores_c_odd].mean()
    Pressure_drop_network_c_odd = mean_pressure_inlet_c_odd - mean_pressure_outlet_c_odd
    Q_throats_new_HC_c_odd, Abs_pres_diff_new_HC_c_odd = cf_pres_fit.Throat_flowrate_total_hydraulic_conductance(net_c, 
                                                                                                                 throat_area='throat.area',
                                                                                                                 pore_diameter='pore.diameter',
                                                                                                                 throat_diameter='throat.diameter',
                                                                                                                 Hydraulic_conductance='throat.hydraulic_conductance')
    # Check if there is conservation of mass
    Q_throats_new_HC_c_odd_inlet = Q_throats_new_HC_c_odd[inlet_throats_c_odd].sum()
    Q_throats_new_HC_c_odd_outlet = Q_throats_new_HC_c_odd[outlet_throats_c_odd].sum()
    
    # Compute the permeability following Darcy's law
    mu_w = param['catholyte_viscosity']
    A_in_c = W_c * H_c
    K = Q_throats_new_HC_c_odd_inlet*mu_w*L_c/(A_in_c*Pressure_drop_network_c_odd)

    '''-------------------------------- Plotting -------------------------------'''
    # Selecting internal pores
    pores_internal_c = net_c.pores('internal')
    
    # Histogram of absolute flowrate 
    plt.figure(1)
    bins = np.linspace(0, 5e-11, 200)
    plt.hist(x=catholyte['throat.absolute_flow_rate_odd_network'][net_c.throats('internal')], bins=bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. odd.')
    plt.title(label = 'Throat absolute flowrate histogram (Internal)')
    plt.xlabel('Throat absolute flowrate')
    plt.ylabel('# Throats')    
    plt.legend()
    
    # Histogram of hydraulic conductance 
    plt.figure(2)
    bins = np.linspace(0, 1e-11, 200)
    plt.hist(x=catholyte['throat.hydraulic_conductance'][net_c.throats('internal')], bins=bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. odd.')
    plt.title(label = 'Throat hydraulic conductance Histogram (Internal)')
    plt.xlabel('Throat hydraulic conductance')
    plt.ylabel('# Throats')    
    plt.legend()
    
    # Mass transfer coefficients (logarithmic plot)
    fig, ax = plt.subplots(figsize=(8, 4))
    bins = np.logspace(np.log10(7e-7),np.log10(1.0e-3),200)
    n_Vel_ind, bins_vel_ind, patches = ax.hist(MT_OpenPNM_cat_odd_vel_independent[pores_internal_c], bins, label='Velocity independent cat. odd. (TFM)')
    n_Vel_dep, bins_vel_dep, patches = ax.hist(MT_OpenPNM_cat_odd_vel_dependent[pores_internal_c], bins, label='Velocity dependent cat. odd. (OpenPNM)')
    
    # Creating the normalized origin plot. 
    # 1. Select the pores with a km between 1e-6 to 1e-3 
    # 2. Normalize these pore to the total number of pores
    # 3. You can copy paste bins_vel_ind (the x-axis) and the NORMALIZED velocity dependent
    # distribution (n_Vel_dep_nrom) and velocity independent distribution (n_Vel_ind_norm) to Origin.
    norm_number_vel_ind = len(np.where(MT_OpenPNM_cat_odd_vel_independent[pores_internal_c]<1e-3)[0])
    norm_number_vel_dep = len(np.where(MT_OpenPNM_cat_odd_vel_dependent[pores_internal_c]>1e-6)[0])
    n_Vel_ind_norm = n_Vel_ind / norm_number_vel_ind
    n_Vel_dep_nrom = n_Vel_dep / norm_number_vel_dep
    ax.set_xscale('log')
    plt.title(label = 'Mass transfer coefficients of internal pores')
    plt.xlabel('Mass transfer coefficient [m/s]')
    plt.ylabel('No. pores')    
    plt.legend()

# Load in the *altered* networks. Moreover, initialize the phase, and stokesflow
# and compute the mass transfer coefficients for the velocity dependent and independent case. 
elif Original_network == 0:
    Folder = '\\Input_networks\\'
    
    # Load in the altered network
    file_name_net_c = network_name + '_catholyte_network_throat_diameter_scaled.pnm'
    file_name_net_a = network_name + '_anolyte_network_throat_diameter_scaled.pnm'
    path_net_c = cwd + Folder + file_name_net_c
    path_net_a = cwd + Folder + file_name_net_a
    project_cat_net = op.Workspace().load_project(filename = path_net_c)
    project_ano_net = op.Workspace().load_project(filename = path_net_a)
    net_c = project_cat_net.network
    net_a = project_ano_net.network
    
    ####### Create catholyte and anolyte phase for the odd and even networks ######
    Conductivity_value = param['anolyte_conductivity']
    anolyte = op.phase.Water(network=net_a, name='anolyte')
    anolyte.add_model(propname = 'pore.electrical_conductivity',                # Anolyte electrical conductivity [S m-1]
                      model = cf_phase.custom_electrical_conductivity_fit,
                      conductivity_fitted = Conductivity_value)    
    anolyte.add_model(propname = 'pore.viscosity',                              # Anolyte viscosity [Pa s]   
                      model = cf_phase.custom_viscosity,
                      parameter_script = param, prop = 'anolyte_viscosity' )
    anolyte.add_model(propname = 'pore.density',                                # Anolyte density [kg m-3]  
                      model = cf_phase.custom_density,
                      parameter_script = param, prop = 'anolyte_density')
    anolyte.add_model(propname = 'pore.diffusivity',                            # Anolyte active species diffusivity in electrolyte [m2 s-1]
                      model = cf_phase.custom_diffusivity,
                      parameter_script = param, prop = 'D_a')            
    anolyte['pore.molecular_weight'] = 0.01802
    
    catholyte = op.phase.Water(network=net_c, name='catholyte')   
    catholyte.add_model(propname = 'pore.electrical_conductivity',              # Catholyte electrical conductivity [S m-1]
                      model = cf_phase.custom_electrical_conductivity_fit,
                      conductivity_fitted = Conductivity_value) 
    catholyte.add_model(propname = 'pore.diffusivity',                          # Catholyte active species diffusivity in electrolyte [m2 s-1]
                      model = cf_phase.custom_diffusivity,
                      parameter_script = param, prop = 'D_c')           
    catholyte.add_model(propname = 'pore.viscosity',                            # Catholyte viscosity [Pa s]
                      model = cf_phase.custom_viscosity,
                      parameter_script = param, prop = 'catholyte_viscosity' )                 
    catholyte.add_model(propname = 'pore.density',                              # Catholyte density [kg m-3]  
                      model = cf_phase.custom_density,
                      parameter_script = param, prop = 'catholyte_density')
    catholyte['pore.molecular_weight'] = 0.01802
    anolyte.regenerate_models()
    catholyte.regenerate_models()  
        
    ################### Define and compute network properties #####################
    # Define electrode dimensions [x=0, y=1, z=2]
    # NOTE: Make sure you know which dimension correlates with which face in your image.
    H_dim = param['height_dimension']   # The dimension corresponding to the height of the electrodes (sides of the electrode)
    L_dim = param['length_dimension']   # The dimension corresponding to the length of the electrodes (flow field direction)
    W_dim = param['width_dimension']    # The dimension corresponding to the width of the electrodes (thickness: current collector -> membrane)            
    
    L_c = np.amax(net_c['pore.coords'][net_c.pores(), L_dim])  # Length [m]
    H_c = np.amax(net_c['pore.coords'][net_c.pores(), H_dim])  # Height [m]
    W_c = np.amax(net_c['pore.coords'][net_c.pores(), W_dim])  # Width [m]
    L_a = np.amax(net_a['pore.coords'][net_a.pores(), L_dim])
    H_a = np.amax(net_a['pore.coords'][net_a.pores(), H_dim])
    W_a = np.amax(net_a['pore.coords'][net_a.pores(), W_dim])

    ########################### Importing Stokesflow ##############################
    # Creating stokes flow algorithm
    sf_c_odd = op.algorithms.StokesFlow(network=net_c, phase=catholyte)
    sf_c_even = op.algorithms.StokesFlow(network=net_c, phase=catholyte)
    sf_a_odd = op.algorithms.StokesFlow(network=net_a, phase=anolyte)
    sf_a_even = op.algorithms.StokesFlow(network=net_a, phase=anolyte)
    
    # Setting input paths
    Folder_SF = '\\Input_Fluid_field\\'
    file_name_sf_c_odd = 'SF_' + str(v_in*100).replace(".", "_") + 'cms_catholyte_odd_network'
    file_name_sf_a_odd = 'SF_' + str(v_in*100).replace(".", "_") + 'cms_anolyte_odd_network'
    file_name_sf_c_even = 'SF_' + str(v_in*100).replace(".", "_") + 'cms_catholyte_even_network'
    file_name_sf_a_even = 'SF_' + str(v_in*100).replace(".", "_") + 'cms_anolyte_even_network'
    
    path_sf_c_odd = cwd + Folder_SF + file_name_sf_c_odd + '.pnm'
    path_sf_a_odd = cwd + Folder_SF + file_name_sf_a_odd + '.pnm'
    path_sf_c_even = cwd + Folder_SF + file_name_sf_c_even + '.pnm'
    path_sf_a_even = cwd + Folder_SF + file_name_sf_a_even + '.pnm'
       
    # Importing stokes flow simulation data
    cf_pres_fit.Import_stokesflow_algo(sf_c_odd, 1, path_sf_c_odd)
    cf_pres_fit.Import_stokesflow_algo(sf_a_odd, 1, path_sf_a_odd)
    cf_pres_fit.Import_stokesflow_algo(sf_c_even, 0, path_sf_c_even)
    cf_pres_fit.Import_stokesflow_algo(sf_a_even, 0, path_sf_a_even)
    
    '''------------------ Mass transfer coefficient ----------------------------'''
    # The Mass transfer equation is: km = C1 * vel^(C2) 
    
    # Update phases with Odd networks
    cf_pres_fit.Update_phase_with_SF(sf_c_odd, catholyte, 1, 1)
    cf_pres_fit.Update_phase_with_SF(sf_a_odd, anolyte, 1, 1)
    
    # Initializing mass transfer coefficient
    catholyte['pore.km_exp'] = np.nan
    
    # Velocity-independent. Note that only the network properties influence km, 
    # it shows no depedency on the phyics used.
    cf_pres_fit.Mass_transfer_coefficient(net_c, catholyte, 1, 0, 0, param)
    MT_OpenPNM_cat_odd_vel_independent = catholyte['pore.km_exp'].copy()
    
    cf_pres_fit.Mass_transfer_coefficient(net_a, anolyte, 1, 0, 0, param)
    MT_OpenPNM_ano_odd_vel_independent = anolyte['pore.km_exp'].copy()
    
    ############################# ODD NETWORKS ####################################
    # OpenPNM physics velocity-dependent
    cf_pres_fit.Mass_transfer_coefficient(net_c, catholyte, 1, 1, 0, param)
    MT_OpenPNM_cat_odd_vel_dependent = catholyte['pore.km_exp'].copy()
    
    cf_pres_fit.Mass_transfer_coefficient(net_a, anolyte, 1, 1, 0, param)
    MT_OpenPNM_ano_odd_vel_dependent = anolyte['pore.km_exp'].copy()
    
    # Larachi physics velocity-dependent
    cf_pres_fit.Mass_transfer_coefficient(net_c, catholyte, 1, 1, 1, param)
    MT_Larachi_cat_odd_vel_dependent = catholyte['pore.km_exp'].copy()
    
    cf_pres_fit.Mass_transfer_coefficient(net_a, anolyte, 1, 1, 1, param)
    MT_Larachi_ano_odd_vel_dependent = anolyte['pore.km_exp'].copy()
    
    ############################# EVEN NETWORKS ###################################
    # Update phases with even networks
    cf_pres_fit.Update_phase_with_SF(sf_c_even, catholyte, 0, 1)
    cf_pres_fit.Update_phase_with_SF(sf_a_even, anolyte, 0, 1)
    
    # OpenPNM physics velocity-dependent
    cf_pres_fit.Mass_transfer_coefficient(net_c, catholyte, 0, 1, 0, param)
    MT_OpenPNM_cat_even_vel_dependent = catholyte['pore.km_exp'].copy()
    
    cf_pres_fit.Mass_transfer_coefficient(net_a, anolyte, 0, 1, 0, param)
    MT_OpenPNM_ano_even_vel_dependent = anolyte['pore.km_exp'].copy()
    
    
    # Larachi physics velocity-dependent
    cf_pres_fit.Mass_transfer_coefficient(net_c, catholyte, 0, 1, 1, param)
    MT_Larachi_cat_even_vel_dependent = catholyte['pore.km_exp'].copy()
    
    cf_pres_fit.Mass_transfer_coefficient(net_a, anolyte, 0, 1, 1, param)
    MT_Larachi_ano_even_vel_dependent = anolyte['pore.km_exp'].copy()

    '''---------------------- Extract fluid field data -------------------------'''
    # Obtain permeability in the flow direction
    inlet_throats_c_odd = net_c.find_neighbor_throats(pores = net_c.pores('flow_inlet'))
    Inlet_surface_pores_c_odd = net_c.find_connected_pores(throats = inlet_throats_c_odd)[:,0]
    outlet_throats_c_odd = net_c.find_neighbor_throats(pores = net_c.pores('flow_outlet'))
    Outlet_surface_pores_c_odd = net_c.find_connected_pores(throats = outlet_throats_c_odd)[:,0]
    # Extracting pressure
    mean_pressure_inlet_c_odd = sf_c_odd['pore.pressure'][Inlet_surface_pores_c_odd].mean()
    mean_pressure_outlet_c_odd = sf_c_odd['pore.pressure'][Outlet_surface_pores_c_odd].mean()
    Pressure_drop_network_c_odd = mean_pressure_inlet_c_odd - mean_pressure_outlet_c_odd
    Q_throats_new_HC_c_odd, Abs_pres_diff_new_HC_c_odd = cf_pres_fit.Throat_flowrate_total_hydraulic_conductance(net_c, 
                                                                                                                 throat_area='throat.area',
                                                                                                                 pore_diameter='pore.diameter',
                                                                                                                 throat_diameter='throat.diameter',
                                                                                                                 Hydraulic_conductance='throat.hydraulic_conductance')
    # Check if there is conservation of mass
    Q_throats_new_HC_c_odd_inlet = Q_throats_new_HC_c_odd[inlet_throats_c_odd].sum()
    Q_throats_new_HC_c_odd_outlet = Q_throats_new_HC_c_odd[outlet_throats_c_odd].sum()
    
    # Compute the permeability following Darcy's law
    mu_w = param['catholyte_viscosity']
    A_in_c = W_c * H_c
    K = Q_throats_new_HC_c_odd_inlet*mu_w*L_c/(A_in_c*Pressure_drop_network_c_odd)
    
    # Obtain Euler number 
    pores_internal_c = net_c.pores('internal')
    n = net_c['pore.Fitting_parameter_n'][0]
    m = net_c['pore.Fitting_parameter_m'][0]
    Gamma = 1
    Eu_inv = cf_pres_fit.Euler_number(net_c, Gamma)
    
    '''-------------------------------- Plotting -------------------------------'''
    # Selecting internal pores
    pores_internal_c = net_c.pores('internal')
    
    # Histogram of absolute flowrate 
    plt.figure(1)
    bins = np.linspace(0, 5e-11, 200)
    plt.hist(x=catholyte['throat.absolute_flow_rate_OpenPNM_odd_network'][net_c.throats('internal')], bins=bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. odd.')
    plt.hist(x=catholyte['throat.absolute_flow_rate_odd_network'][net_c.throats('internal')], bins=bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. odd.')
    plt.title(label = 'Throat absolute flowrate histogram (Internal)')
    plt.xlabel('Throat absolute flowrate')
    plt.ylabel('# Throats')    
    plt.legend()
    
    # Euler number (cumulative logarithmic plot)
    n_bins = np.logspace(np.log10(1e-6),np.log10(1.0),200)
    fig, ax = plt.subplots(figsize=(8, 4))
    n, bins, patches = ax.hist(Eu_inv, n_bins, density=True, histtype='step',
                                cumulative=-1, label='Empirical')
    ax.set_xscale('log')
    plt.xlabel('Inverse Euler number')
    plt.ylabel('Normalized no. throats')   
    
    # Histogram of hydraulic conductance 
    plt.figure(3)
    bins = np.linspace(0, 1e-11, 200)
    plt.hist(x=catholyte['throat.hydraulic_conductance'][net_c.throats('internal')], bins=bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. odd.')
    plt.hist(x=catholyte['throat.hydraulic_conductance'][net_c.throats('internal')], bins=bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. odd.')
    plt.title(label = 'Throat hydraulic conductance Histogram (Internal)')
    plt.xlabel('Throat hydraulic conductance')
    plt.ylabel('# Throats')    
    plt.legend()
    
    # Mass transfer coefficients (logarithmic plot)
    fig, ax = plt.subplots(figsize=(8, 4))
    bins = np.logspace(np.log10(7e-7),np.log10(1.0e-3),200)
    n_Vel_ind, bins_vel_ind, patches = ax.hist(MT_OpenPNM_cat_odd_vel_independent[pores_internal_c], bins, label='Velocity independent cat. odd. (TFM)')
    n_Vel_dep, bins_vel_dep, patches = ax.hist(MT_Larachi_cat_odd_vel_dependent[pores_internal_c], bins, label='Velocity dependent cat. odd. (Larachi)')
    
    # Creating the normalized origin plot. 
    # 1. Select the pores with a km between 1e-6 to 1e-3 
    # 2. Normalize these pore to the total number of pores
    # 3. You can copy paste bins_vel_ind (the x-axis) and the NORMALIZED velocity dependent
    # distribution (norm_number_vel_dep) and velocity independent distribution (norm_number_vel_ind) to Origin.
    norm_number_vel_ind = len(np.where(MT_OpenPNM_cat_odd_vel_independent[pores_internal_c]<1e-3)[0])
    norm_number_vel_dep = len(np.where(MT_Larachi_cat_odd_vel_dependent[pores_internal_c]>1e-6)[0])
    n_Vel_ind_norm = n_Vel_ind / norm_number_vel_ind
    n_Vel_dep_nrom = n_Vel_dep / norm_number_vel_dep
    ax.set_xscale('log')
    plt.title(label = 'Mass transfer coefficients of internal pores')
    plt.xlabel('Mass transfer coefficient [m/s]')
    plt.ylabel('No. pores')    
    plt.legend()