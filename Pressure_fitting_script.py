"""
Script fits the pressure drop with the conduit conductance model of 
F. Larachi et al.   - http://dx.doi.org/10.1016/j.cej.2013.11.077 
X. Liu et al.       - https://doi.org/10.1002/aic.16258

This is done by fitting three parameters:
    1.	The contraction curvature parameter N
    2.	The expansion curvature parameter M
    3.	The throat effective aspect ratio F

Note, however, that for the IDFF we only fit the first factor due to the linear
pressure drop regime in which the IDFF resides. 

Working principle:
    1. Initialization of network. For each conduit containing a 
       throat diameter > pore diameter:
           Assign the smallest pore diameter as throat diameter.
    
    2. For each velocity - load in network in pressure optimization function. 
        a. Create a new network and a new phase. A new network is created for 
           fitting of the throat effective aspect ratio F, as it is defined w.r.t.
           the ORIGINAL throat radius. 
        
        b. Create an initial pressure field guess based on the conduit conductance
            model of OpenPNM.

        c. Run Stokes-flow algo with "new" conduit conductance. Note 
           that the transport.py has been altered to include its computation. After 
           every pressure computation the conduit conductance is updated and
           assigned to the throats.
           
        d. Iteratively solve the system of equations by a direct substitution 
           approach. Extract the throat flow rate and compare with throat flow rate
          of previous iteration. If the mean and maximum deviation is within
          the threshold the simulation is stopped. If not, the iterations continue. 
    
    3. Script re-runs the optimization funtion for all velocities and obtains the 
       final fitting parameters.
    
Options script
    1. SNOW
        If SNOW == 1 -> Use a network extracted with SNOW 1
        If SNOW == 2 -> Use a network extracted with SNOW 2
        
    2. Inlet_channel
        If Inlet_channel == 0 -> Assign inlet, rib and outlet region with boundary conditions
        If Inlet_channel == 1 -> Assign inlet, rib and outlet region with flow field pores 
        
    3. Flow field type
        If Flow_field == 0: The FTFF is simulated.
        If Flow_field == 1: The IDFF is simulated.
"""

SNOW = 1
Inlet_channel = 0
Flow_field = 0

"""------------ Importing Packages and simulation constants-----------------"""
import openpnm as op
import numpy as np
import time
import inputDictTEMPO as inputDict  
import os
import Custom_functions_km as cf
import Costum_network as cn
import Costum_functions_transport as cf_trans
import Costum_functions_phase as cf_phase
import Custom_functions_pressure_fitting as cf_pres_fit
import pandas as pd
import lmfit

Gamma = 1                                       # Flow pattern constant (1 = flat velocity profile, 2 = Parabolic velocity profile)
C0 = 26                                         # Laminar constant for contraction
E0 = 27                                         # Laminar constant for expansion
param = inputDict.input_dict                    # Load in the input dictionary

############ Setting import network directory and loading network##############
cwd = os.getcwd()                               # Curent Working directory
network_path = '\\input\\'                      # Folder containing network

# Network selection
network_name = 'Freudenberg_check_CSV'
# network_name = 'Freudenberg_MC_pore_inscribed_throat_inscribed'

file_vec = cwd + network_path + network_name    # Directory of network
file_name = file_vec + '.pnm'
project_cat = op.Workspace().load_project(filename = file_name)
net = project_cat.network 

########################## Setting input parameters ###########################
"""----------------------- Adjusting throat diameter -----------------------"""
# Throats with Dt> Dp are reset ot Dt = Dp (smallest of each connecting pore)
throats_internal = net.throats('internal')
throats_boundary = net.throats('boundary')
Dt_altered = cn.Shrink_throat_radius_pressure_fitting(net, 
                                                      pore_diameter='pore.diameter',
                                                      throat_diameter='throat.diameter')
net['throat.diameter'][throats_internal] = Dt_altered

################### Define and compute network properties #####################
# Define electrode dimensions [x=0, y=1, z=2]
# NOTE: Make sure you know which dimension correlates with which face in your image.
H_dim = param['height_dimension']               # The dimension corresponding to the height of the electrodes (sides of the electrode)
L_dim = param['length_dimension']               # The dimension corresponding to the length of the electrodes (flow field direction)
W_dim = param['width_dimension']                # The dimension corresponding to the width of the electrodes (thickness: current collector -> membrane)            

L_c = np.amax(net['pore.coords'][net.pores(), L_dim])       # Length [m]
H_c = np.amax(net['pore.coords'][net.pores(), H_dim])       # Height [m]
W_c = np.amax(net['pore.coords'][net.pores(), W_dim])       # Width [m]
H_c_min = np.amin(net['pore.coords'][net.pores(), H_dim])   # Minimum Height [m] (offset from y = 0)

# Assign the labels 'membrane', 'current collector' 'flow inlet' and 'flow outlet'
# To the correct boundary pores.
if Flow_field == 0:
    cf.assign_boundary_pores(net, W_dim, L_dim)
if Flow_field == 1:
    cf.assign_boundary_pores_IDFF(net, W_dim, L_dim, H_dim, H_c, H_c_min, Inlet_channel)

A_ext_c = L_c * H_c                             # Cathode area to compute external current density [m2]

if Flow_field == 0:
    A_in_c = W_c * H_c                          # Cathode inlet area [m2]
if Flow_field == 1:
    A_in_c = W_c * L_c                          # Cathode inlet area [m2]. Note for IDFF: this is the thickness of the electrodetimes its length

######################## Setting throat boundary labels #######################
# Labeling boundary throats. Used in check for mass conservation
membrane_throats = net.find_neighbor_throats(pores = net.pores('membrane'))
current_collector_throats = net.find_neighbor_throats(pores = net.pores('current_collector'))
right_throats_boundary = net.find_neighbor_throats(pores = net.pores('right'))
left_throats_boundary = net.find_neighbor_throats(pores = net.pores('left'))
net.set_label(label='membrane', throats = membrane_throats)
net.set_label(label='current_collector', throats = current_collector_throats)
net.set_label(label='right', throats = right_throats_boundary)
net.set_label(label='left', throats = left_throats_boundary)

'''----------------------- Start optimization function ---------------------'''
def Fitting_pressure_drop(params, net, v_in_vec, P_drop_vec_exp, param, mean_pressure_inlet, 
                          mean_pressure_outlet, Pressure_drop_network, del_P_internal, num_net, Pressure_drop_electrode, Ori_thr_diam, Ori_thr_area, A_ext_c, A_in_c, SNOW, Gamma):
                          
    """-------------------------Fitting Factors---------------------------------"""
    # Copy original throat diameter and area (Required for fitting Throat_diameter_factor)
    net['throat.diameter'] = np.copy(Ori_thr_diam)
    net['throat.area'] = np.copy(Ori_thr_area)
    
    Throat_diameter_factor = params['Throat_diameter_factor'].value   # Factor = Effective_throat_radius / Real_throat_radius
    n_factor_fit = params['Contraction_factor'].value                 # Contraction curvature fit factor
    m_factor_fit = params['Expansion_factor'].value                   # Expansion curvature fit factor 
    
    print('Throat diameter fit factor:', Throat_diameter_factor)
    print('Contraction curvature fit factor:', n_factor_fit)
    print('Expansion curvature fit factor :', m_factor_fit)
    
    """------------- Initializing network and electrolyte ----------------------"""
    ####################### Initializing fitting network ##########################
    net_fit = op.network.Network()
    net_fit.update(net)
    
    # Adjusting internal throat diameter
    net_fit['throat.diameter'][net_fit.throats('internal')] = net_fit['throat.diameter'][net_fit.throats('internal')] * Throat_diameter_factor
    
    ########################## Adjusting throat diameter ######################
    # Throats with Dt> Dp are reset ot Dt = Dp (smallest of each connecting pore)
    throats_internal = net_fit.throats('internal')
    throats_boundary = net_fit.throats('boundary')
    Dt_altered = cn.Shrink_throat_radius_pressure_fitting(net_fit, 
                                                          pore_diameter='pore.diameter',
                                                          throat_diameter='throat.diameter')
    net_fit['throat.diameter'][throats_internal] = Dt_altered

    ################### Define and compute phase properties #######################
    catholyte = op.phase.Water(network=net_fit, name='catholyte')
    catholyte.add_model(propname = 'pore.electrical_conductivity',              # Catholyte electrical conductivity [S m-1]
                      model = cf_phase.costum_electrical_conductivity,
                      parameter_script = param, prop = 'catholyte_conductivity')    
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
    
    ############################ Adding phase physics #############################
    if SNOW == 2:
        net_fit['throat.area'] = net_fit['throat.cross_sectional_area']
        net_fit['pore.area'] = net_fit['pore.surface_area']
    
    f_ctc = cn.Center_to_center_length
    net_fit.add_model(propname = 'throat.ctc_length', model = f_ctc)    
    f_endpoints = cn.end_points_spherical_pores
    net_fit.add_model(propname = 'throat.endpoints', model = f_endpoints)    
    f_conduit_length = cn.spherical_pores_conduit_lengths
    net_fit.add_model(propname = 'throat.conduit_lengths', model = f_conduit_length)   
    f_hyd = cf_trans.Flow_shape_factors_ball_and_stick
    catholyte.add_model(propname = 'throat.flow_shape_factors', model = f_hyd)   
    f_Hyd_cond = cf_trans.Hydraulic_conductance_Hagen_Poiseuille
    catholyte.add_model(propname = 'throat.hydraulic_conductance', model = f_Hyd_cond)
    # Make a hard copy of the OpenPNM hydraulic conductance values
    Hydraulic_conductance_OpenPNM = catholyte['throat.hydraulic_conductance'].copy()
    catholyte['throat.hydraulic_conductance_OP_copy'] = Hydraulic_conductance_OpenPNM.copy()
    
    #################### Setting up Stokes-flow algorithm #########################
    # Run Stokes Flow algorithm. The pressure drop over all networks in the NIS
    # approach are the same, so we only need to fit the data for one network.
    sf = op.algorithms.StokesFlow(network=net_fit, phase=catholyte)
    
    # Assigning the inlet flow rate. Done based on the fractional area of the 
    # internal surface pores - instead of the equally sized boundary pores.
    # Finding the inlet throats of even and odd networks. 
    inlet_throats= net_fit.find_neighbor_throats(pores = net_fit.pores('flow_inlet'))
    
    # Find the inlet internal surface pores
    Inlet_surface_pores = net_fit.find_connected_pores(throats = inlet_throats)[:,0]
    
    # The flow rate entering a pore has to be scaled with the area of the inlet pore.
    cross_area_c = np.pi / 4 * net_fit['pore.diameter'] ** 2
    
    # total_area_inlet_c = np.sum(cross_area_c[net.pores('flow_inlet')])
    total_area_internal_inlet_surface_pores = np.sum(cross_area_c[Inlet_surface_pores])
    
    # Loop over the different velocities and extract the pressure drop
    idx = 0
    print(('Vel [m/s]\t P_model [bar]\t P_exp [bar] \tP_model/P_exp [-]'))
    
    for v_in_c in v_in_vec:
        
        # Regenerate OpenPNM hydraulic conductance 
        catholyte.regenerate_models()
        # set inlet flowrate
        Q_in_c = v_in_c * A_in_c # total flow rate entering the network [m3/s]
        
        # Imposing the Neumann inlet boundary condition 
        for pore in net_fit.pores('flow_inlet'):
                    # Find the throat connection between flow inlet boundary pore and internal surface pore
                    Throat_inlet_Q = net_fit.find_neighbor_throats(pores = pore)
                    # Find the connected internal surface pore
                    Connected_internal_surface_pore = net_fit.find_connected_pores(throats = Throat_inlet_Q)[:,0]
                    # Extract the pore area of the corresponding internal surface pore
                    Area_internal_surface_pore = cross_area_c[Connected_internal_surface_pore]
                    # Assign the corresponding flow rate (Neumann inlet boundary condition)
                    sf.set_rate_BC(rates=Q_in_c * Area_internal_surface_pore/total_area_internal_inlet_surface_pores, pores=pore, mode = 'overwrite') 
        
        # Dirichlet outlet boundary condition
        Pout = 0    # Pressure outlet boundary condition [Pa]
        sf.set_value_BC(values=Pout, pores=net_fit.pores('flow_outlet'), mode = 'overwrite')    
        
        ############## Compute stokesflow for OpenPNM physics #################
        # Assign value of n, m, Gamma and Init to the network, and run the Stokesflow algorithm
        net_fit['pore.Fitting_parameter_n'] = n_factor_fit
        net_fit['pore.Fitting_parameter_m'] = m_factor_fit
        net_fit['pore.Fitting_parameter_Gamma'] = Gamma
        net_fit['pore.parameter_Init'] = 0   # Indicate we work with hydraulic conductance of OpenPNM     
        sf.run()
        catholyte.update(sf.soln)
        
        # Find the throat flowrates and mean flowrate
        Q_throats_OpenPNM_physics, Abs_pres_diff_OpenPNM_physics = cf_pres_fit.Throat_flowrate_total_hydraulic_conductance(net_fit, 
                                                                                                                           throat_area='throat.area',
                                                                                                                           pore_diameter='pore.diameter',
                                                                                                                           throat_diameter='throat.diameter',
                                                                                                                           Hydraulic_conductance='throat.hydraulic_conductance')
        
        Q_throats_OpenPNM_physics_mean = Q_throats_OpenPNM_physics.mean()
        
        ###### Compute pressure drop over network with OpenPNM physics ########
        Inlet_surface_pores = net_fit.find_connected_pores(throats = inlet_throats)[:,0]
        outlet_throats = net_fit.find_neighbor_throats(pores = net_fit.pores('flow_outlet'))
        Outlet_surface_pores = net_fit.find_connected_pores(throats = outlet_throats)[:,0]
        
        # Mean pressure inlet and outlet
        mean_pressure_inlet = sf['pore.pressure'][Inlet_surface_pores].mean()
        mean_pressure_outlet = sf['pore.pressure'][Outlet_surface_pores].mean()
        Pressure_drop_network = mean_pressure_inlet - mean_pressure_outlet
        maximum_pressure_network = sf['pore.pressure'].max()
        
        Init_guess = sf['pore.pressure'].copy() + 1e-3
        Init_guess[np.where(Init_guess == 1e-3)[0]] = 0
        
        """--- Setting and solving Stokes-flow algorithm with new Hydraulic model---"""    
        # Run stokesflow algorithm with the new hydraulic model
        net_fit['pore.parameter_Init'] = 1          # Indicate we work with hydraulc conductance of Larachi et al.
        # sf.run(x0=Init_guess)                     # Compute with an intitial pressure guess based on the converged results from OpenPNM physics stokesflow
        sf.run()                                    # Compute with an initial pressure guess of zero
        # Find the throat flowrates and mean flowrate
        Q_throats_new_HC, Abs_pres_diff_new_HC = cf_pres_fit.Throat_flowrate_total_hydraulic_conductance(net_fit, 
                                                                                                         throat_area='throat.area',
                                                                                                         pore_diameter='pore.diameter',
                                                                                                         throat_diameter='throat.diameter',
                                                                                                         Hydraulic_conductance='throat.hydraulic_conductance')
        
        Q_throats_new_HC_mean = Q_throats_new_HC.mean()
        
        ############ Pressure drop over network Larachi Physics ###############
        Inlet_surface_pores = net_fit.find_connected_pores(throats = inlet_throats)[:,0]
        outlet_throats = net_fit.find_neighbor_throats(pores = net_fit.pores('flow_outlet'))
        Outlet_surface_pores = net_fit.find_connected_pores(throats = outlet_throats)[:,0]
        
        # Mean pressure inlet and outlet
        mean_pressure_inlet_new_HC = sf['pore.pressure'][Inlet_surface_pores].mean()
        mean_pressure_outlet_new_HC= sf['pore.pressure'][Outlet_surface_pores].mean()
        Pressure_drop_network_new_HC = mean_pressure_inlet_new_HC - mean_pressure_outlet_new_HC
        maximum_pressure_network_new_HC = sf['pore.pressure'].max()
        
        # Save the pressure drop over the entire electrode per inlet velocity
        Pressure_drop_electrode[:,idx] = num_net * Pressure_drop_network_new_HC

        # Print the results to the console        
        print(f'{v_in_c:.2f}\t\t {Pressure_drop_electrode[0,idx]/1e5:.3f}\t\t\t {P_drop_vec_exp[idx]/1e5:.3f}\t\t\t {Pressure_drop_electrode[0,idx]/P_drop_vec_exp[idx]:.3f}\t')
        
        # Index for new velocity
        idx = idx + 1
        
    # Delete network and catholyte
    del net_fit, catholyte, sf
    
    model = Pressure_drop_electrode
    return model - P_drop_vec_exp

""" ------------------Fitting script------------------------------------"""
starting_time = time.time()

# Import experimental velocity and pressure drop data
df = pd.read_excel('.//input//Experiment_fitting_pressure.xlsx.xlsx') #
v_in_vec = df['Velocity'].to_numpy()                    # Experimental Velocity [m/s]
P_drop_vec_exp = df['Pressure_drop'].to_numpy()         # Experimental pressure drop [Pa]

# Number of networks
if Flow_field == 0:
    num_net = int(np.ceil(param['total_electrode_length']/L_c)) 
if Flow_field == 1:
    num_net = 1
    
# Initialization of saving the data
mean_pressure_inlet = np.zeros([1, len(v_in_vec)])
mean_pressure_outlet = np.zeros([1, len(v_in_vec)])
Pressure_drop_network = np.zeros([1, len(v_in_vec)])
Pressure_drop_electrode = np.zeros([1, len(v_in_vec)])
maximum_pressuree_network = np.zeros([1, len(v_in_vec)])
del_P_internal = np.zeros([len(v_in_vec), 3])

method = 'leastsq'           

params_fit = lmfit.Parameters()
params_fit.add('Throat_diameter_factor', 1.2618734059716306, min=0.5, max=2.5)      # Factor = Effective_throat_radius / Real_throat_radius
params_fit.add('Contraction_factor', 0.6769651781769254, min=0.1, max=2.5)          # Contraction curvature fit factor
params_fit.add('Expansion_factor', 1.3592818186419697, min=0.1, max=2.5)            # Expansion curvature fit factor

# Extract the original (though adapted for Dt>Dp) throat and pore diameter and area
Ori_thr_diam_ = net['throat.diameter'].copy()
Ori_thr_area_ = net['throat.area'].copy()
Ori_thr_diam = np.copy(Ori_thr_diam_)
Ori_thr_area = np.copy(Ori_thr_area_)

out = lmfit.minimize(Fitting_pressure_drop, params = params_fit, method=method, args=(net, v_in_vec, P_drop_vec_exp, param, mean_pressure_inlet, mean_pressure_outlet, Pressure_drop_network, del_P_internal, num_net, Pressure_drop_electrode, Ori_thr_diam, Ori_thr_area, A_ext_c, A_in_c, SNOW, Gamma))
                                               
lmfit.report_fit(out)    
final_time = time.time()
print(f'Network simulations finished in {np.floor((final_time-starting_time)/60):.0f} minutes and {(final_time-starting_time)%60:.1f} seconds.')