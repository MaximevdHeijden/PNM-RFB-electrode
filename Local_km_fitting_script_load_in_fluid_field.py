"""   
Script fits the pre-exponential (alpha) and exponential (Beta) factor of a 
*LOCAL* mass transfer correlation (the form of which given below) to a global mass 
transfer correlation (obtained from e.g. literature or limiting current density experiments).
    km_{local} = alpha * velocity_{pore} ^ {Beta}.
Here velocity_{pore} is the velocity in the *pores*, which is derived following the 
methodology presented by Larachi et al. (http://dx.doi.org/10.1016/j.cej.2013.11.077 )

The methodology of the local km fitting is adapted from:
    Chapter 4: "Structure-performance relationship in multi-layer VRFB electrodes"
    of the PhD Thesis "Multiphysics pore network simulation of electrochemical
    devices as a design tool" by N. Misaghian ( http://hdl.handle.net/10012/18998 ).
    AND     
    "Mass transfer in fibrous media with varying anisotropy for flow battery 
    electrodes: Direct numerical simulations with 3D X-ray computed tomography""
    by M. Kok et al. (https://doi.org/10.1016/j.ces.2018.10.049 ).
     
     In summary, the script numerically repeats a limiting current density
     experiment. For a set range of inlet velocities it computes global mass transfer
     values, and compares that to value obtained experimentally (Ilim experiments)
     or from that obtained from a mass transfer correlation published in literature. 
     
     It should be noted that the used approach assumes the bulk concentration to 
     remain constant throughout the electrode (i.e. equal to the inlet concentration).
     To minimize depletion, only 6 to 8 % of all pores are selected for a reaction 
     to occur. This corresponds to about 11-14 % of the total surface area.
     This process is repeated "Itermax" times for randomly selected group of 
     pores, and the final global km is obtained through averaging the results of these iterations.
     
Working principle:
    1. Loading in the network 
    
    2. Initializing the phase and stokesflow object. Note that you can choose
        between employing the conduit conductance model of OpenPNM (shape factor)
        that by Larachi et al. (mechanical energy balance based on fitting pressure drop)
        with the option *Larachi_HC*, further detailed below. All stokesflow 
        algorithms have been created and ran in the function script "Flow_field_main", 
        placed in the folder "Hydrodynamic_fitting". 
    
    3. Creating the advection-diffusion-reaction algorithm. A reaction source
       term is ascribed to the function:
           R_i [mol/s] = km_{local} * A_i * (C_i_bulk - C_i_surface).
           At I_Lim: C_i_surface is assumed zero. Hence rxn rate becomes:
           R_i [mol/s] = km_{local} * A_i * (C_i_bulk - 0).
           
     4. Fitting the km_{local} pre-exponential (alpha) and exponential (Beta). 
     
         a. Load in the global mass transfer values at defined inlet velocities, 
            and transfer the network (and associated algorithms) to the fitting function. 
         
         b. Run the advection-diffusion-reaction algorithm. Extract the total 
            reaction rate at each inlet velocity and numerically compute the
            global mass transfer value.
       
Options:
    1. Larachi_HC
        If Larachi_HC == 0: OpenPNM physics are used (i.e. shape factor)
        If Larachi_HC == 1: Physics based on mechanical energy balance are used (Larachi et al.)
        
    2. Flow field type
        If Flow_field == 0: The FTFF is simulated.
        If Flow_field == 1: The IDFF is simulated. 
"""

Larachi_HC = 1
Itermax = 10                                    # No. iterations for averaging of global km per applied velocity
Flow_field = 0

""" ---------Importing Packages and set Load/Read directories --------------"""
import openpnm as op
import numpy as np
import inputDictTEMPO as inputDict  
import os
import Custom_functions_km as cf
import Costum_functions_transport as cf_trans
import Costum_functions_phase as cf_phase
import Custom_functions_pressure_fitting as cf_pres_fit
import time
import openpyxl
import random

# Setting import network directory
cwd = os.getcwd()                               # Curent Working directory
network_path = '\\Input_networks\\'             # Folder containing network
network_name = 'Freudenberg_check_CSV'
file_name_net_c = network_name + '_catholyte_network_throat_diameter_scaled.pnm'
path_net_c = cwd + network_path + file_name_net_c
file_vec = [path_net_c]
param = inputDict.input_dict                    # Load in the input dictionary
v_in = 1.5e-2

"""----------Loading in the Network and simulation constants----------------"""
for NrNetwork, file in enumerate(file_vec):
    
    op.Workspace().clear()  
    
    # Loading in network
    file_name = file + '.pnm'
    project_cat = op.Workspace().load_project(filename = path_net_c)
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
    H_dim = param['height_dimension']           # The dimension corresponding to the height of the electrodes (sides of the electrode)
    L_dim = param['length_dimension']           # The dimension corresponding to the length of the electrodes (flow field direction)
    W_dim = param['width_dimension']            # The dimension corresponding to the width of the electrodes (thickness: current collector -> membrane)            
    
    L_c = np.amax(net_c['pore.coords'][net_c.pores(), L_dim])  # Length [m]
    H_c = np.amax(net_c['pore.coords'][net_c.pores(), H_dim])  # Height [m]
    W_c = np.amax(net_c['pore.coords'][net_c.pores(), W_dim])  # Width [m]         
    
    A_ext_c = L_c * H_c                          # Cathode area to compute external current density [m2]
    if Flow_field == 0:
        A_in_c = W_c * H_c                       # Cathode inlet area [m2]
    if Flow_field == 1:
        A_in_c = W_c * L_c                       # Cathode inlet area [m2] Note for IDFF: this is the thickness of the electrodetimes its length
    mem_area = A_ext_c                           # Fictitious membrane area [m2]            
    SA_network = net_c['pore.surface_area']      # Network surface area [m2]    
    Total_network_area = SA_network.sum()        # Total network surface area [m2]
    
    """-----------Define and compute phase properties---------------"""    
    # Phase constants are assigned via models, as during phase
    # regeneration the water object inherent functions would overwrite 
    # assigned phase constants.
    Conductivity_value = param['anolyte_conductivity']
    catholyte = op.phase.Water(network=net_c, name='catholyte')  
    catholyte.add_model(propname = 'pore.electrical_conductivity',              # Catholyte electrical conductivity [S m-1]
                      model = cf_phase.costum_electrical_conductivity_fit,
                      conductivity_fitted = Conductivity_value) 
    catholyte.add_model(propname = 'pore.diffusivity',                          # Catholyte active species diffusivity in electrolyte [m2 s-1]
                      model = cf_phase.costum_diffusivity,
                      parameter_script = param, prop = 'D_c')           
    catholyte.add_model(propname = 'pore.viscosity',                            # Catholyte viscosity [Pa s]
                      model = cf_phase.costum_viscosity,
                      parameter_script = param, prop = 'catholyte_viscosity')                 
    catholyte.add_model(propname = 'pore.density',                              # Catholyte density [kg m-3]  
                      model = cf_phase.costum_density,
                      parameter_script = param, prop = 'catholyte_density')
    catholyte['pore.molecular_weight'] = 0.01802
    
    """-----------------------Stokes algorithm----------------------"""
    # Loading the pre-run stokes algorithm. Note that some phase models require
    # certain phase properties (i.e advection-diffusion required hydraulic conductance). 
    # Hence, the SF-algorithm should be initialized prior to assigning these models.
    # In the iterative loop we will update the stokesflow algorithm based on the flow velocity.
    # The StokesFlow loaded in here will  thus not be used for the other velocities.
    
    # Creating stokes flow algorithm
    sf_c_odd = op.algorithms.StokesFlow(network=net_c, phase=catholyte)
    
    # Simply to initialize the ad_dif_conductance of catholyte. 
    Folder_SF = '\\Input_Fluid_field\\'
    file_name_sf_c_odd = 'SF_' + str(v_in*100).replace(".", "_") + 'cms_catholyte_odd_network'
    path_sf_c_odd = cwd + Folder_SF + file_name_sf_c_odd + '.pnm'
    cf_pres_fit.Import_stokesflow_algo(sf_c_odd, 1, path_sf_c_odd)
    
    if Larachi_HC == 0:
        cf_pres_fit.Update_phase_with_SF(sf_c_odd, catholyte, 1, 0)
        
    elif Larachi_HC == 1:
        cf_pres_fit.Update_phase_with_SF(sf_c_odd, catholyte, 1, 1) 
    
    #------------------------ Adding phase physics ---------------------------#
    f_poisson = cf_trans.Poisson_shape_factors_ball_and_stick
    catholyte.add_model(propname = 'throat.poisson_shape_factors', model = f_poisson)
    f_diff_cond = cf_trans.Diffusive_conductance_mixed_diffusion
    catholyte.add_model(propname = 'throat.diffusive_conductance', model = f_diff_cond)
    f_ad_dif = cf_trans.Advection_diffusion
    catholyte.add_model(propname = 'throat.ad_dif_conductance', model = f_ad_dif)
    f_elec_con = cf_trans.Electrical_conductance_series_resistors
    catholyte.add_model(propname = 'throat.electrical_conductance', model = f_elec_con)            
    catholyte.regenerate_models()            
    
    """----------------- Mass and charge transport models ------------------""" 
    # Set-up advection-Diffusion algorithm
    ad_c_odd = op.algorithms.AdvectionDiffusion(network=net_c, phase=catholyte)
            
    # Setting boundary conditions:
    ad_c_odd.set_outflow_BC(pores=net_c.pores('flow_outlet'))
    
    # Source term: Limiting current density experiment. Pore rxn rate equates
    # to mass transfer to the pore:
    # R_i [mol/s] = km_{local} * A_i * (C_i_bulk - C_i_surface).
    # At I_Lim: C_i_surface is assumed zero. Hence rxn rate becomes:
    # R_i [mol/s] = km_{local} * A_i * (C_i_bulk - 0).
    # The  linear source term is defined as: r = A_{1} X + A_{2}. The inter-
    # cept  A_{2} = 0. A_{1} will be defined in the iterative loop. The 
    # independent variable X is the concentration inside the pores. 
    
    source_term = op.models.physics.source_terms.linear
    catholyte['pore.rxn_A2'] = 0.0
    catholyte['pore.rxn_A1'] = 0.0
    catholyte.add_model(propname='pore.reaction_term', 
                        model=source_term,
                        A1='pore.rxn_A1', 
                        A2="pore.rxn_A2", 
                        X='pore.concentration',
                        regen_mode="deferred") 
    ad_c_odd.set_source(propname='pore.reaction_term', pores=net_c.pores('internal'))

    """--------------- Optimization function mass transport ----------------"""
    def iterative_solver(params, catholyte, net_c, data,
                         ad_c_odd, SA_network, Total_network_area, sf_c_odd, number_of_networks, param, Larachi_HC, v_in_vec, Itermax):                         
            
        alpha = params['alpha_factor'].value                # Pre-exponential factor of *Local* MT correlation
        beta = params['beta_factor'].value                  # Exponential factor of *Local* MT correlation
        print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\nPre-exponential factor (alpha) is:', alpha,'\nExponential factor (Beta) is:', beta ,'\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
        
        # Initializing matrices used for storing global km
        Global_km_vec = np.zeros([1, len(v_in_vec)])        # Saving the averaged global km per applied velocity 
        Global_vec_int = np.zeros([1, Itermax])             # Saving the global km per network at each inlet velocity       
        Global_km_vec_stdev = np.zeros([1, len(v_in_vec)])  # Saving the Standard deviation of computed global km per applied velocity

        # Matrixes to store inlet concentration for the network
        conc0_in_c_odd = np.zeros((len(net_c.pores('flow_inlet')), 1))
        conc_in_c0 = param['conc_in_c']                     # Catholyte active species inlet concentration in network [mol m-3]
        conc0_in_c_odd[:, 0] = conc_in_c0                   # Set the inlet concentration of the 0th network to the inlet concentration of the electrode    
        
        """----------- Start of fitting local km coefficients  -------------"""
        # In this algorithm we loop over each network at a specified inlet
        # velocity and numerically compute the global km.
        
        vel_idx = 0                                         # Index for velocity
        
        for v_in in v_in_vec:
            
            # Setting input paths to load in Stokesflow algorithm for specified 
            # inlet velocity.
            Folder_SF = '\\Input_Fluid_field\\'
            file_name_sf_c_odd = 'SF_' + str(v_in*100).replace(".", "_") + 'cms_catholyte_odd_network'
            path_sf_c_odd = cwd + Folder_SF + file_name_sf_c_odd + '.pnm'
            
            # Importing StokesFlow algorithm simulation data
            cf_pres_fit.Import_stokesflow_algo(sf_c_odd, 1, path_sf_c_odd)
            
            network = 0
            for Iter in range(Itermax):
                # Reset network concentration profile and assign inlet concentrations. 
                conc0_in_c_odd[:, 0] = conc_in_c0                    
                conc_in_c = conc0_in_c_odd[:, network]
                    
                # Selection of reactive pockets. We select 1% of random pores 
                # and find their neighboring pores (mean conductivity Freudenberg ~ 6).
                # All these are subsequently labeled as reactive pores. About 
                # 11 - 14 % of the total internal surface area is labeled as reactive 
                # with this approach. 
                Fraction = 0.01
                No_reactive_pores = int(len(net_c.pores('internal')) * Fraction)
                Reactive_pores = random.sample(range(len(net_c.pores('internal'))), No_reactive_pores)
                # Find the internal neighbor pores of the random pores.
                Neighbor_pores = net_c.find_neighbor_pores(pores = Reactive_pores)
                net_c.set_label(label='pore.neighbor', pores=Neighbor_pores)
                Neighbor_pores_internal = net_c.pores(['neighbor', 'internal'], mode='and')
                # Concetatenate randomly selected and neighbor pores into a single array
                Reactive_pores_all = np.concatenate((Reactive_pores, Neighbor_pores_internal))
                # Remove the label of neighboring pores
                net_c.set_label(label='pore.neighbor', mode='purge')
                # Label reactive pores
                net_c.set_label(label='pore.reactive', pores=Reactive_pores_all)
                
                # Reactive surface area fraction is computed as: 
                # SA_network[net_c.pores('reactive')].sum() / SA_network[net_c.pores('internal')].sum()
                
                ####### Solving advection-diffusion-reaction algorithm ########
                # Retrieve the StokesFlow data
                catholyte.update(sf_c_odd.soln)
                
                if Larachi_HC == 0:                     # OpenPNM conduit conductance model
                    cf_pres_fit.Update_phase_with_SF(sf_c_odd, catholyte, 1, 0)
                
                    # Local MT coefficient
                    cf_pres_fit.Local_mass_transfer_coefficient_fitting(net_c, catholyte, 1, alpha, beta)
                    
                elif Larachi_HC == 1:                   # Larachi et al. conduit conductance model
                    cf_pres_fit.Update_phase_with_SF(sf_c_odd, catholyte, 1, 1) 
                
                    # Local MT coefficient
                    cf_pres_fit.Local_mass_transfer_coefficient_fitting(net_c, catholyte, 1, alpha, beta)  

                # Update the Linear reaction source term : R_i = km_local * A_i * (C_b_i - 0)
                # Note that only the reactive pores are assigned the term A1.
                catholyte['pore.rxn_A1'][net_c.pores('reactive')] = - catholyte['pore.km_loc'][net_c.pores('reactive')] * SA_network[net_c.pores('reactive')]
                catholyte['pore.rxn_A2'] = 0
                
                # Regenerate (the ad_dif_conductance) models and set boundary conditions
                catholyte.regenerate_models('throat.ad_dif_conductance')
                ad_c_odd.set_value_BC(pores = net_c.pores('internal', mode = 'nor'), mode='remove')
                ad_c_odd.set_value_BC(values=conc_in_c, pores=net_c.pores('flow_inlet'), mode='overwrite')
                ad_c_odd.set_outflow_BC(pores = net_c.pores('internal', mode = 'nor'), mode='remove')
                ad_c_odd.set_outflow_BC(pores=net_c.pores('flow_outlet'), mode='overwrite')
                
                # Solve the advection-diffusion-reaction equation
                ad_c_odd.run()      
                    
                # Update the concentration
                conc_c = ad_c_odd['pore.concentration']
                    
                # Update the value for the concentration
                catholyte['pore.concentration'] = conc_c     
        
                ######### Compute global km and update inlet concentration new network ##########
                # Total reaction rate is computed as:
                # km_{global} = Sum(R_i)/(Sum(A_i) * C_in)
                # with R_i [mol/s] = km_{local} * A_i * (C_i_bulk - 0).
                Reaction_rate = - catholyte['pore.reaction_term.rate']
                Total_reaction_rate = Reaction_rate.sum()
                Surface_area_reactive_pores = SA_network[net_c.pores('reactive')].sum()
                Global_vec_int[:,Iter] = Total_reaction_rate / (Surface_area_reactive_pores * param['conc_in_c'])
                
                # Uncomment for visualization of network (with reactive pores indicated) in Paraview:
                # pn_cat_vtk = op.network.Network()
                # pn_cat_vtk.update(net_c)
                # catholyte_vtk = op.phase.Phase(network=pn_cat_vtk, name='catholyte_vtk')
                # catholyte_vtk.update(catholyte)
                # op.io.project_to_vtk(project = pn_cat_vtk.project, filename = '.\\output\\' + 'Cathode' )
                                
                # Reset the reactive pores and their reactive term
                catholyte['pore.rxn_A1'][net_c.pores('reactive')]  = 0
                net_c.set_label(label='pore.reactive', mode='purge')
                
            # Compute mean of global km of each network for each inlet velocity
            Global_km_vec[:,vel_idx] = Global_vec_int.mean()
            Global_km_vec_stdev[:,vel_idx] = np.std(Global_vec_int)   

            vel_idx = vel_idx + 1
            
        model = Global_km_vec          
        return model - data     
            
    import pandas as pd
    import lmfit
    
    """ --------------------- Fitting script -------------------------------"""
    starting_time = time.time()                 # Start recording fitting time
    num_net = 1                                 # Only fit for a single network

    # Import experimental global km at specified inlet velocities
    df = pd.read_excel('.//input//Experiment_fitting_masstransfer.xlsx')
    v_in_vec = df['velocity'].to_numpy() 
    Global_km_vec = df['km_exp'].to_numpy() 
    data = Global_km_vec
    method = 'leastsq'        
        
    # Intializing fitting parameters
    params = lmfit.Parameters()
    params.add('alpha_factor', 0.00016, min= 1.42789011e-5, max=1.42789011e-3) 
    params.add('beta_factor', 0.82 , min = 0.6, max = 1.0) 

    out = lmfit.minimize(iterative_solver, params, method=method, args=(catholyte, net_c, data,
                                                                        ad_c_odd, SA_network, Total_network_area, sf_c_odd, num_net, param, Larachi_HC, v_in_vec, Itermax))
    lmfit.report_fit(out)    
    final_time = time.time()
    print(f'Network simulations finished in {np.floor((final_time-starting_time)/60):.0f} minutes and {(final_time-starting_time)%60:.1f} seconds.')