"""
Script outputs the fluid field (i.e. results from the stokesflow algorithm) to 
be loaded into e.g. scripts used for polarization simulations or for obtaining 
the local mass transfer coefficient. It does this for:
    1. Catholyte even networks
    2. Catholyte odd networks
    3. Anolyte even networks
    4. Anlyte odd networks

The script computes the fluid fields using the hydraulic conductance models of both:
    1. OpenPNM physics (i.e. with shape factors)
    2. The hydraulic model presented by 
        F. Larachi et al.   - http://dx.doi.org/10.1016/j.cej.2013.11.077 
        X. Liu et al.       - https://doi.org/10.1002/aic.16258
Inside the respective scripts where the fluid fiels are loaded in, you can assign
which physics you want to use. 

Working principle:
    1. Set a velocity for which you want the fluid field. Further, set the fitted
        contraction curvature parameter N, expansion curvature parameter M, and
        throat effective aspect ratio F obtained in the "Pressure fitting script".
    
    2. Initialization of network and phase. The throat diameters of the network
       are altered based on the assigned throat effective aspect ratio F. Hereafter
       for each conduit containing a throat diameter > pore diameter:
           The smallest pore diameter is assgined as throat diameter.
         
    3. Compute the fluid fields of the *odd* networks for the anolyte and catholyte. 
        This is done for both physics, and throat/pore velocity are extracted.
        
    4. Compute the fluid fields of the *even* networks for the anolyte and catholyte. 
        This is done for both physics, and throat/pore velocity are extracted.
    
    5. Export the phases (i.e. results of the stokesflow algorithms) and the
       *altered* network (altered throat diameter).
    
The fluid fields are saved in the folder "Input_Fluid_field"
The altered networks are saved in the folder "Input_networks"
    
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
import inputDictTEMPO as inputDict  
import os
import Custom_functions_km as cf
import Costum_network as cn
import Costum_functions_transport as cf_trans
import Costum_functions_phase as cf_phase
import matplotlib.pyplot as plt
import Custom_functions_pressure_fitting as cf_pres_fit

Gamma = 2                                           # Flow pattern constant (1 = flat velocity profile, 2 = Parabolic velocity profile)
C0 = 26                                             # Laminar constant for contraction
E0 = 27                                             # Laminar constant for expansion
param = inputDict.input_dict                        # Load in the input dictionary
v_in_vec = [1.5e-2]                                 # Set inlet velocity [m/s]

# Set fitting factors
Throat_diameter_factor = 1.2618734059716306         # Factor = Effective_throat_radius / Real_throat_radius
n_factor_fit = 0.6769651781769254                   # Contraction curvature fit factor
m_factor_fit =  1.3592818186419697                  # Expansion curvature fit factor   
print('Throat diameter fit factor:', Throat_diameter_factor)
print('Contraction curvature fit factor:', n_factor_fit)
print('Expansion curvature fit factor :', m_factor_fit)

############ Setting import network directory and loading network##############
cwd = os.getcwd()                                   # Curent Working directory
network_path = '\\input\\'                          # Folder containing network

# Network selection
network_name = 'Freudenberg_check_CSV'

########################## Setting input parameters ###########################
"""----------Loading in the Network and simulation constants----------------"""

file_vec = cwd + network_path + network_name    # Directory of network
file_name = file_vec + '.pnm'

project_cat = op.Workspace().load_project(filename = file_name)
net_c = project_cat.network 

project_ano = op.Workspace().load_project(filename = file_name)
net_a = project_ano.network  

# There are three pores which cause the stokesflow algorithm not to converge 
# for the new hydraulic computation in the Freudenberg_check_CSV network. These
# pores are trimmed. 
if Flow_field == 0:
    op.topotools.trim(net_c, pores = [2609, 3458, 5080])
    op.topotools.trim(net_a, pores = [2609, 3458, 5080])
if Flow_field == 1:
    op.topotools.trim(net_c, pores = [5507, 751, 9217])
    op.topotools.trim(net_a, pores = [5507, 751, 9217])
    
# Checking for isolated and disconnected pores:
health_c = op.utils.check_network_health(network = net_c)
print("The health of the cathodic and anodic electrodes are as",
      "follows (note that they are the same mirrored network):\n")
print(health_c); 

# Update the face labels that are defined in the SNOW algorithm
cf.validate_face_labels(net_c) 
cf.validate_face_labels(net_a)

"""----------------------- Adjusting throat diameter -----------------------"""
# Throats with Dt> Dp are reset ot Dt = Dp (smallest of each connecting pore)
throats_internal_c = net_c.throats('internal')
throats_boundary_c = net_c.throats('boundary')
Original_throat_diameter_cat = net_c['throat.diameter'][throats_internal_c].copy()
Dt_altered = cn.Shrink_throat_radius(net_c, 0, Original_throat_diameter_cat,
                                      pore_diameter='pore.diameter',
                                      throat_diameter='throat.diameter')
net_c['throat.diameter'][throats_internal_c] = Dt_altered
Original_throat_diameter_cat = net_c['throat.diameter'][throats_internal_c].copy()

throats_internal_a = net_a.throats('internal')
throats_boundary_a = net_a.throats('boundary')
Original_throat_diameter_ano = net_a['throat.diameter'][throats_internal_c].copy()
Dt_altered = cn.Shrink_throat_radius(net_a, 0,Original_throat_diameter_ano,
                                      pore_diameter='pore.diameter',
                                      throat_diameter='throat.diameter')
net_a['throat.diameter'][throats_internal_a] = Dt_altered
Original_throat_diameter_ano = net_a['throat.diameter'][throats_internal_c].copy()

################### Define and compute network properties #####################
# Define electrode dimensions [x=0, y=1, z=2]
# NOTE: Make sure you know which dimension correlates with which face in your image.
H_dim = param['height_dimension']   # The dimension corresponding to the height of the electrodes (sides of the electrode)
L_dim = param['length_dimension']   # The dimension corresponding to the length of the electrodes (flow field direction)
W_dim = param['width_dimension']    # The dimension corresponding to the width of the electrodes (thickness: current collector -> membrane)            

L_c = np.amax(net_c['pore.coords'][net_c.pores(), L_dim])       # Maximum Length [m]
H_c = np.amax(net_c['pore.coords'][net_c.pores(), H_dim])       # Maximum Height [m]
H_c_min = np.amin(net_c['pore.coords'][net_c.pores(), H_dim])   # Minimum Height [m] (offset from y = 0)
W_c = np.amax(net_c['pore.coords'][net_c.pores(), W_dim])       # Maximum Width [m]
L_a = np.amax(net_a['pore.coords'][net_a.pores(), L_dim])
H_a = np.amax(net_a['pore.coords'][net_a.pores(), H_dim])
H_a_min = np.amin(net_a['pore.coords'][net_a.pores(), H_dim])
W_a = np.amax(net_a['pore.coords'][net_a.pores(), W_dim])

# Assign the labels 'membrane', 'current collector' 'flow inlet' and 'flow outlet'
# To the correct boundary pores.
if Flow_field == 0:
    cf.assign_boundary_pores(net_c, W_dim, L_dim)
    cf.assign_boundary_pores(net_a, W_dim, L_dim)  
if Flow_field == 1:
    cf.assign_boundary_pores_IDFF(net_c, W_dim, L_dim, H_dim, H_c, H_c_min, Inlet_channel)
    cf.assign_boundary_pores_IDFF(net_a, W_dim, L_dim, H_dim, H_a, H_a_min, Inlet_channel)

# Position networks next to each other (left: anode, middle: fictitious membrane, right: cathode)
net_a['pore.coords'] -= [W_a * 1.5, 0, 0]

A_ext_c = L_c * H_c                                 # Cathode area to compute external current density [m2]
A_ext_a = L_a * H_a                                 # Anode area [m2]
mem_area = A_ext_c                                  # Fictitious membrane area [m2]  

# Invert anodic network in the width in order to obtain the same pores at the membrane side for both electrodes
inversion_factor = net_a['pore.coords'][:, W_dim] - min(net_a['pore.coords'][:, W_dim])
net_a['pore.coords'][:, W_dim]  = max(net_a['pore.coords'][:, W_dim]) - inversion_factor

if Flow_field == 0:
    A_in_c = W_c * H_c                              # Cathode inlet area [m2]
    A_in_a = W_a * H_a                              # Anode inlet area [m2]
    
    # The boundaries of the network are also reversed.
    pores_left = net_a['pore.right']
    pores_right = net_a['pore.left']
    net_a['pore.right'] = pores_right
    net_a['pore.left'] = pores_left
    
elif Flow_field == 1:
    A_in_c = W_c * L_c                              # Cathode inlet area [m2]
    A_in_a = W_a * L_a                              # Anode inlet area [m2]

"""------------- Initializing network and electrolyte ----------------------"""
# Adjusting internal throat diameter
net_c['throat.diameter'][net_c.throats('internal')] = net_c['throat.diameter'][net_c.throats('internal')] * Throat_diameter_factor
net_a['throat.diameter'][net_a.throats('internal')] = net_a['throat.diameter'][net_a.throats('internal')] * Throat_diameter_factor

########################## Adjusting throat diameter ######################
# Throats with Dt> Dp are reset ot Dt = Dp (smallest of each connecting pore)
throats_internal_c = net_c.throats('internal')
throats_boundary_c = net_c.throats('boundary')
Dt_altered = cn.Shrink_throat_radius(net_c, 1, Original_throat_diameter_cat,
                                      pore_diameter='pore.diameter',
                                      throat_diameter='throat.diameter')

net_c['throat.diameter'][throats_internal_c] = Dt_altered

throats_internal_a = net_a.throats('internal')
throats_boundary_a = net_a.throats('boundary')
Dt_altered = cn.Shrink_throat_radius(net_a, 1, Original_throat_diameter_ano,
                                      pore_diameter='pore.diameter',
                                      throat_diameter='throat.diameter')

net_a['throat.diameter'][throats_internal_a] = Dt_altered

################### Define and compute phase properties #######################
anolyte = op.phase.Water(network=net_a, name='anolyte')
anolyte.add_model(propname = 'pore.electrical_conductivity',                # Anolyte electrical conductivity [S m-1]
                  model = cf_phase.costum_electrical_conductivity,
                  parameter_script = param, prop = 'anolyte_conductivity')            
anolyte.add_model(propname = 'pore.viscosity',                              # Anolyte viscosity [Pa s]   
                  model = cf_phase.costum_viscosity,
                  parameter_script = param, prop = 'anolyte_viscosity')
anolyte.add_model(propname = 'pore.density',                                # Anolyte density [kg m-3]  
                  model = cf_phase.costum_density,
                  parameter_script = param, prop = 'anolyte_density')
anolyte.add_model(propname = 'pore.diffusivity',                            # Anolyte active species diffusivity in electrolyte [m2 s-1]
                  model = cf_phase.costum_diffusivity,
                  parameter_script = param, prop = 'D_a')            
anolyte['pore.molecular_weight'] = 0.01802

catholyte = op.phase.Water(network=net_c, name='catholyte')
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
    net_c['throat.area'] = net_c['throat.cross_sectional_area']
    net_c['pore.area'] = net_c['pore.surface_area']
    net_a['throat.area'] = net_a['throat.cross_sectional_area']
    net_a['pore.area'] = net_a['pore.surface_area']

f_ctc = cn.Center_to_center_length
net_c.add_model(propname = 'throat.ctc_length', model = f_ctc)   
net_a.add_model(propname = 'throat.ctc_length', model = f_ctc)  

f_endpoints = cn.end_points_spherical_pores
net_c.add_model(propname = 'throat.endpoints', model = f_endpoints)    
net_a.add_model(propname = 'throat.endpoints', model = f_endpoints)    

f_conduit_length = cn.spherical_pores_conduit_lengths
net_c.add_model(propname = 'throat.conduit_lengths', model = f_conduit_length)   
net_a.add_model(propname = 'throat.conduit_lengths', model = f_conduit_length)   

f_hyd = cf_trans.Flow_shape_factors_ball_and_stick
catholyte.add_model(propname = 'throat.flow_shape_factors', model = f_hyd)   
anolyte.add_model(propname = 'throat.flow_shape_factors', model = f_hyd)   

f_Hyd_cond = cf_trans.Hydraulic_conductance_Hagen_Poiseuille
catholyte.add_model(propname = 'throat.hydraulic_conductance', model = f_Hyd_cond)
anolyte.add_model(propname = 'throat.hydraulic_conductance', model = f_Hyd_cond)

# Make a hard copy of the OpenPNM hydraulic conductance values
Hydraulic_conductance_OpenPNM_cat = catholyte['throat.hydraulic_conductance'].copy()
Hydraulic_conductance_OpenPNM_ano = anolyte['throat.hydraulic_conductance'].copy()

#################### Setting up Stokes-flow algorithm #########################
# Run Stokes Flow algorithm. The pressure drop over all networks in the NIS
# approach are the same, so we only need to obtain the data for one network.

r''' Run Stokes Flow algorithm in both half cells for both the odd
and the even networks. Boundary conditions for the even networks are 
reversed so that the outlet pores of e.g. network 0 can be directly
mapped to the inlet pores of network 1 (it are the same pores).'''            
sf_c_odd = op.algorithms.StokesFlow(network=net_c, phase=catholyte)
sf_c_even = op.algorithms.StokesFlow(network=net_c, phase=catholyte)
sf_a_odd = op.algorithms.StokesFlow(network=net_a, phase=anolyte)
sf_a_even = op.algorithms.StokesFlow(network=net_a, phase=anolyte)

# Assigning the inlet flow rate. Done based on the fractional area of the 
# internal surface pores - instead of the equally sized boundary pores.
# Finding the inlet throats of even and odd networks. 
if Flow_field == 0:
    inlet_throats_c_odd = net_c.find_neighbor_throats(pores = net_c.pores('flow_inlet'))
    inlet_throats_a_odd = net_a.find_neighbor_throats(pores = net_a.pores('flow_inlet'))
    inlet_throats_c_even = net_c.find_neighbor_throats(pores=net_c.pores('flow_outlet'))
    inlet_throats_a_even = net_a.find_neighbor_throats(pores=net_a.pores('flow_outlet'))  
if Flow_field == 1:
    inlet_throats_c_odd = net_c.find_neighbor_throats(pores = net_c.pores('flow_inlet'))
    inlet_throats_a_odd = net_a.find_neighbor_throats(pores = net_a.pores('flow_inlet'))
    inlet_throats_c_even = net_c.find_neighbor_throats(pores=net_c.pores('flow_inlet'))
    inlet_throats_a_even = net_a.find_neighbor_throats(pores=net_a.pores('flow_inlet')) 

# Find the inlet internal surface pores
Inlet_surface_pores_c_odd = net_c.find_connected_pores(throats = inlet_throats_c_odd)[:,0]
Inlet_surface_pores_a_odd = net_a.find_connected_pores(throats = inlet_throats_a_odd)[:,0]
Inlet_surface_pores_c_even = net_c.find_connected_pores(throats = inlet_throats_c_even)[:,0]
Inlet_surface_pores_a_even = net_a.find_connected_pores(throats = inlet_throats_a_even)[:,0]

# The flow rate entering a pore has to be scaled with the area of the inlet pore.
cross_area_c = np.pi / 4 * net_c['pore.diameter'] ** 2
cross_area_a = np.pi / 4 * net_a['pore.diameter'] ** 2

# total_area_inlet_c = np.sum(cross_area_c[net.pores('flow_inlet')])
total_area_internal_inlet_surface_pores_c_odd = np.sum(cross_area_c[Inlet_surface_pores_c_odd])
total_area_internal_inlet_surface_pores_a_odd = np.sum(cross_area_a[Inlet_surface_pores_a_odd])
total_area_internal_inlet_surface_pores_c_even = np.sum(cross_area_c[Inlet_surface_pores_c_even])
total_area_internal_inlet_surface_pores_a_even = np.sum(cross_area_a[Inlet_surface_pores_a_even])

for v_in_c in v_in_vec:
    
    # Regenerate OpenPNM hydraulic conductance 
    catholyte.regenerate_models()
    anolyte.regenerate_models()

    # set inlet flowrate
    Q_in_c = v_in_c * A_in_c # total flow rate entering the network [m3/s]
    Q_in_a = v_in_c * A_in_a # total flow rate entering the network [m3/s] 
    
    # Imposing the Neumann inlet boundary condition 
    # Catholyte odd
    for pore in net_c.pores('flow_inlet'):
               # Find the throat connection between flow inlet boundary pore and internal surface pore
               Throat_inlet_Q_c_odd = net_c.find_neighbor_throats(pores = pore)
               # Find the connected internal surface pore
               Connected_internal_surface_pore_c_odd = net_c.find_connected_pores(throats = Throat_inlet_Q_c_odd)[:,0]
               # Extract the pore area of the corresponding internal surface pore
               Area_internal_surface_pore_c_odd = cross_area_c[Connected_internal_surface_pore_c_odd]
               # Assign the corresponding flow rate (Neumann inlet boundary condition)
               sf_c_odd.set_rate_BC(rates=Q_in_c * Area_internal_surface_pore_c_odd/total_area_internal_inlet_surface_pores_c_odd, pores=pore, mode = 'overwrite') 
   
    # Catholyte even
    if Flow_field == 0:
        for pore in net_c.pores('flow_outlet'):
                   # Find the throat connection between flow inlet boundary pore and internal surface pore
                   Throat_inlet_Q_c_even = net_c.find_neighbor_throats(pores = pore)
                   # Find the connected internal surface pore
                   Connected_internal_surface_pore_c_even = net_c.find_connected_pores(throats = Throat_inlet_Q_c_even)[:,0]
                   # Extract the pore area of the corresponding internal surface pore
                   Area_internal_surface_pore_c_even = cross_area_c[Connected_internal_surface_pore_c_even]
                   # Assign the corresponding flow rate (Neumann inlet boundary condition)
                   sf_c_even.set_rate_BC(rates=Q_in_c * Area_internal_surface_pore_c_even/total_area_internal_inlet_surface_pores_c_even, pores=pore, mode = 'overwrite') 
    if Flow_field == 1:
        for pore in net_c.pores('flow_inlet'):
                   # Find the throat connection between flow inlet boundary pore and internal surface pore
                   Throat_inlet_Q_c_even = net_c.find_neighbor_throats(pores = pore)
                   # Find the connected internal surface pore
                   Connected_internal_surface_pore_c_even = net_c.find_connected_pores(throats = Throat_inlet_Q_c_even)[:,0]
                   # Extract the pore area of the corresponding internal surface pore
                   Area_internal_surface_pore_c_even = cross_area_c[Connected_internal_surface_pore_c_even]
                   # Assign the corresponding flow rate (Neumann inlet boundary condition)
                   sf_c_even.set_rate_BC(rates=Q_in_c * Area_internal_surface_pore_c_even/total_area_internal_inlet_surface_pores_c_even, pores=pore, mode = 'overwrite') 

    # Anolyte even
    for pore in net_a.pores('flow_inlet'):
               # Find the throat connection between flow inlet boundary pore and internal surface pore
               Throat_inlet_Q_a_odd = net_a.find_neighbor_throats(pores = pore)
               # Find the connected internal surface pore
               Connected_internal_surface_pore_a_odd = net_a.find_connected_pores(throats = Throat_inlet_Q_a_odd)[:,0]
               # Extract the pore area of the corresponding internal surface pore
               Area_internal_surface_pore_a_odd = cross_area_a[Connected_internal_surface_pore_a_odd]
               # Assign the corresponding flow rate (Neumann inlet boundary condition)
               sf_a_odd.set_rate_BC(rates=Q_in_a * Area_internal_surface_pore_a_odd/total_area_internal_inlet_surface_pores_a_odd, pores=pore, mode = 'overwrite') 
 
    # Anolyte even
    if Flow_field == 0:
        for pore in net_a.pores('flow_outlet'):
                   # Find the throat connection between flow inlet boundary pore and internal surface pore
                   Throat_inlet_Q_a_even = net_a.find_neighbor_throats(pores = pore)
                   # Find the connected internal surface pore
                   Connected_internal_surface_pore_a_even = net_a.find_connected_pores(throats = Throat_inlet_Q_a_even)[:,0]
                   # Extract the pore area of the corresponding internal surface pore
                   Area_internal_surface_pore_a_even = cross_area_a[Connected_internal_surface_pore_a_even]
                   # Assign the corresponding flow rate (Neumann inlet boundary condition)
                   sf_a_even.set_rate_BC(rates=Q_in_a * Area_internal_surface_pore_a_even/total_area_internal_inlet_surface_pores_a_even, pores=pore, mode = 'overwrite') 
    if Flow_field == 1:
        for pore in net_a.pores('flow_inlet'):
                   # Find the throat connection between flow inlet boundary pore and internal surface pore
                   Throat_inlet_Q_a_even = net_a.find_neighbor_throats(pores = pore)
                   # Find the connected internal surface pore
                   Connected_internal_surface_pore_a_even = net_a.find_connected_pores(throats = Throat_inlet_Q_a_even)[:,0]
                   # Extract the pore area of the corresponding internal surface pore
                   Area_internal_surface_pore_a_even = cross_area_a[Connected_internal_surface_pore_a_even]
                   # Assign the corresponding flow rate (Neumann inlet boundary condition)
                   sf_a_even.set_rate_BC(rates=Q_in_a * Area_internal_surface_pore_a_even/total_area_internal_inlet_surface_pores_a_even, pores=pore, mode = 'overwrite') 

    # Dirichlet outlet boundary condition
    Pout = 0    # Pressure outlet boundary condition [Pa]
    if Flow_field == 0:
        sf_c_odd.set_value_BC(values=Pout, pores=net_c.pores('flow_outlet'), mode = 'overwrite')    # Dirichlet outlet boundary condition
        sf_c_even.set_value_BC(values=Pout, pores=net_c.pores('flow_inlet'), mode = 'overwrite')    # Dirichlet outlet boundary condition
        sf_a_odd.set_value_BC(values=Pout, pores=net_a.pores('flow_outlet'), mode = 'overwrite')    # Dirichlet outlet boundary condition
        sf_a_even.set_value_BC(values=Pout, pores=net_a.pores('flow_inlet'), mode = 'overwrite')    # Dirichlet outlet boundary condition
    if Flow_field == 1:
        sf_c_odd.set_value_BC(values=Pout, pores=net_c.pores('flow_outlet'), mode = 'overwrite')    # Dirichlet outlet boundary condition
        sf_c_even.set_value_BC(values=Pout, pores=net_c.pores('flow_outlet'), mode = 'overwrite')   # Dirichlet outlet boundary condition
        sf_a_odd.set_value_BC(values=Pout, pores=net_a.pores('flow_outlet'), mode = 'overwrite')    # Dirichlet outlet boundary condition
        sf_a_even.set_value_BC(values=Pout, pores=net_a.pores('flow_outlet'), mode = 'overwrite')   # Dirichlet outlet boundary condition

    '''------------------------ Running odd networks -----------------------'''
    ################# Compute stokesflow for OpenPNM physics ##################
    # Assign value of n, m, Gamma and Init to the network 
    net_c['pore.Fitting_parameter_n'] = n_factor_fit
    net_c['pore.Fitting_parameter_m'] = m_factor_fit
    net_c['pore.Fitting_parameter_Gamma'] = Gamma
    net_c['pore.parameter_Init'] = 0            # Indicate we work with hydraulic conductance of OpenPNM    
    
    net_a['pore.Fitting_parameter_n'] = n_factor_fit
    net_a['pore.Fitting_parameter_m'] = m_factor_fit
    net_a['pore.Fitting_parameter_Gamma'] = Gamma
    net_a['pore.parameter_Init'] = 0            # Indicate we work with hydraulic conductance of OpenPNM    
    
    # Run electrolyte transport algorithms for odd networks
    sf_c_odd.run()
    sf_a_odd.run()
       
    # Find the throat flowrates and velocity, and the pore velocity of odd flow networks and save them to phase
    catholyte.update(sf_c_odd.soln)
    anolyte.update(sf_a_odd.soln)
    Q_throats_OpenPNM_physics_c_odd, Abs_pres_diff_OpenPNM_physics_c_odd = cf_pres_fit.Throat_flowrate_total_hydraulic_conductance(net_c, 
                                                                                                                                   throat_area='throat.area',
                                                                                                                                   pore_diameter='pore.diameter',
                                                                                                                                   throat_diameter='throat.diameter',
                                                                                                                                   Hydraulic_conductance='throat.hydraulic_conductance')
    Q_throats_OpenPNM_physics_c_odd_mean = Q_throats_OpenPNM_physics_c_odd.mean()
    cf_pres_fit.Throat_pore_velocity_extraction(net_c, Q_throats_OpenPNM_physics_c_odd, 0, 1, sf_c_odd)
    Q_throats_OpenPNM_physics_a_odd, Abs_pres_diff_OpenPNM_physics_a_odd = cf_pres_fit.Throat_flowrate_total_hydraulic_conductance(net_a, 
                                                                                                                                   throat_area='throat.area',
                                                                                                                                   pore_diameter='pore.diameter',
                                                                                                                                   throat_diameter='throat.diameter',
                                                                                                                                   Hydraulic_conductance='throat.hydraulic_conductance')
    Q_throats_OpenPNM_physics_a_odd_mean = Q_throats_OpenPNM_physics_a_odd.mean()
    cf_pres_fit.Throat_pore_velocity_extraction(net_a, Q_throats_OpenPNM_physics_a_odd, 0, 1, sf_a_odd)
    
    sf_c_odd['throat.Hydraulic_conductance_OpenPNM'] = Hydraulic_conductance_OpenPNM_cat.copy() 
    sf_a_odd['throat.Hydraulic_conductance_OpenPNM'] = Hydraulic_conductance_OpenPNM_ano.copy()  
    
    ######## Compute pressure drop over network with OpenPNM physics ##########
    # Catholyte odd network
    Inlet_surface_pores_c_odd = net_c.find_connected_pores(throats = inlet_throats_c_odd)[:,0]
    outlet_throats_c_odd = net_c.find_neighbor_throats(pores = net_c.pores('flow_outlet'))
    Outlet_surface_pores_c_odd = net_c.find_connected_pores(throats = outlet_throats_c_odd)[:,0]
    # Mean pressure inlet and outlet
    mean_pressure_inlet_c_odd = sf_c_odd['pore.pressure'][Inlet_surface_pores_c_odd].mean()
    mean_pressure_outlet_c_odd = sf_c_odd['pore.pressure'][Outlet_surface_pores_c_odd].mean()
    Pressure_drop_network_c_odd = mean_pressure_inlet_c_odd - mean_pressure_outlet_c_odd
    maximum_pressure_network_c_odd = sf_c_odd['pore.pressure'].max()
    
    # Anolyte odd network
    Inlet_surface_pores_a_odd = net_a.find_connected_pores(throats = inlet_throats_a_odd)[:,0]
    outlet_throats_a_odd = net_a.find_neighbor_throats(pores = net_a.pores('flow_outlet'))
    Outlet_surface_pores_a_odd = net_a.find_connected_pores(throats = outlet_throats_a_odd)[:,0]
    # Mean pressure inlet and outlet
    mean_pressure_inlet_a_odd = sf_a_odd['pore.pressure'][Inlet_surface_pores_c_odd].mean()
    mean_pressure_outlet_a_odd = sf_a_odd['pore.pressure'][Outlet_surface_pores_c_odd].mean()
    Pressure_drop_network_a_odd = mean_pressure_inlet_a_odd - mean_pressure_outlet_a_odd
    maximum_pressure_network_a_odd = sf_a_odd['pore.pressure'].max()
    
    print('\n ---------------------- ODD NETWORK OpenPNM physics---------------------------\n')
    ################### Outputting Stokes flow algorithm data #####################
    print(f'Pressure_drop_single_network_net_c_odd: {Pressure_drop_network_c_odd*1e-5 :.5f} bar')
    print(f'Pressure_drop_total_network_net_c_odd: {Pressure_drop_network_c_odd*1e-5*17 :.5f} bar')
    print(f'maximum_pressure_network_net_c_odd: {maximum_pressure_network_c_odd*1e-5 :.5f} bar')
    print(f'mean_pressure_inlet_net_c_odd: {mean_pressure_inlet_c_odd*1e-5 :.5f} bar')
    print(f'mean_pressure_outlet_net_c_odd: {mean_pressure_outlet_c_odd*1e-5 :.5f} bar')
    Init_guess_c_odd = sf_c_odd['pore.pressure'].copy() + 1e-3
    Init_guess_c_odd[np.where(Init_guess_c_odd == 1e-3)[0]] = 0    
    
    print(f'Pressure_drop_single_network_net_a_odd: {Pressure_drop_network_a_odd*1e-5 :.5f} bar')
    print(f'Pressure_drop_total_network_net_a_odd: {Pressure_drop_network_a_odd*1e-5*17 :.5f} bar')
    print(f'maximum_pressure_network_net_a_odd: {maximum_pressure_network_a_odd*1e-5 :.5f} bar')
    print(f'mean_pressure_inlet_net_a_odd: {mean_pressure_inlet_a_odd*1e-5 :.5f} bar')
    print(f'mean_pressure_outlet_net_a_odd: {mean_pressure_outlet_a_odd*1e-5 :.5f} bar')
    Init_guess_a_odd = sf_a_odd['pore.pressure'].copy() + 1e-3
    Init_guess_a_odd[np.where(Init_guess_a_odd == 1e-3)[0]] = 0   
    
    """--- Setting and solving Stokes-flow algorithm with new Hydraulic model---"""
    # Run stokesflow algorithm with the new hydraulic model
    net_c['pore.parameter_Init'] = 1            # Indicate we work with hydraulc conductance of Larachi et al.
    # sf_c_odd.run(x0=Init_guess_c_odd)         # Compute with an intitial pressure guess based on the converged results from OpenPNM physics stokesflow
    sf_c_odd.run()                              # Compute with an initial pressure guess of zero

    net_a['pore.parameter_Init'] = 1            # Indicate we work with hydraulc conductance of Larachi et al.
    # sf_a_odd.run(x0=Init_guess_a_odd)         # Compute with an intitial pressure guess based on the converged results from OpenPNM physics stokesflow
    sf_a_odd.run()                              # Compute with an initial pressure guess of zero

    # Find the throat flowrates and velocity, and the pore velocity of odd flow networks and save them to phase
    catholyte.update(sf_c_odd.soln)
    anolyte.update(sf_a_odd.soln)
    Q_throats_new_HC_c_odd, Abs_pres_diff_new_HC_c_odd = cf_pres_fit.Throat_flowrate_total_hydraulic_conductance(net_c, 
                                                                                     throat_area='throat.area',
                                                                                     pore_diameter='pore.diameter',
                                                                                     throat_diameter='throat.diameter',
                                                                                     Hydraulic_conductance='throat.hydraulic_conductance')
    Q_throats_new_HC_c_odd_mean = Q_throats_new_HC_c_odd.mean()
    cf_pres_fit.Throat_pore_velocity_extraction(net_c, Q_throats_new_HC_c_odd, 1, 1, sf_c_odd)
    
    Q_throats_new_HC_a_odd, Abs_pres_diff_new_HC_a_odd = cf_pres_fit.Throat_flowrate_total_hydraulic_conductance(net_a, 
                                                                                     throat_area='throat.area',
                                                                                     pore_diameter='pore.diameter',
                                                                                     throat_diameter='throat.diameter',
                                                                                     Hydraulic_conductance='throat.hydraulic_conductance')
    
    Q_throats_new_HC_a_odd_mean = Q_throats_new_HC_a_odd.mean()
    cf_pres_fit.Throat_pore_velocity_extraction(net_a, Q_throats_new_HC_a_odd, 1, 1, sf_a_odd)
    
    sf_c_odd['throat.Hydraulic_conductance_Larachi'] = catholyte['throat.hydraulic_conductance'].copy()
    sf_a_odd['throat.Hydraulic_conductance_Larachi'] = anolyte['throat.hydraulic_conductance'].copy()  
    
    ############## Pressure drop over network Larachi Physics #################
    # Catholyte odd network
    Inlet_surface_pores_c_odd = net_c.find_connected_pores(throats = inlet_throats_c_odd)[:,0]
    outlet_throats_c_odd = net_c.find_neighbor_throats(pores = net_c.pores('flow_outlet'))
    Outlet_surface_pores_c_odd = net_c.find_connected_pores(throats = outlet_throats_c_odd)[:,0]
    # Mean pressure inlet and outlet
    mean_pressure_inlet_new_HC_c_odd = sf_c_odd['pore.pressure'][Inlet_surface_pores_c_odd].mean()
    mean_pressure_outlet_new_HC_c_odd = sf_c_odd['pore.pressure'][Outlet_surface_pores_c_odd].mean()
    Pressure_drop_network_new_HC_c_odd = mean_pressure_inlet_new_HC_c_odd - mean_pressure_outlet_new_HC_c_odd
    maximum_pressure_network_new_HC_c_odd = sf_c_odd['pore.pressure'].max()
    
    # Anolyte odd network
    Inlet_surface_pores_a_odd = net_a.find_connected_pores(throats = inlet_throats_a_odd)[:,0]
    outlet_throats_a_odd = net_a.find_neighbor_throats(pores = net_a.pores('flow_outlet'))
    Outlet_surface_pores_a_odd = net_a.find_connected_pores(throats = outlet_throats_a_odd)[:,0]
    # Mean pressure inlet and outlet
    mean_pressure_inlet_new_HC_a_odd = sf_a_odd['pore.pressure'][Inlet_surface_pores_c_odd].mean()
    mean_pressure_outlet_new_HC_a_odd = sf_a_odd['pore.pressure'][Outlet_surface_pores_c_odd].mean()
    Pressure_drop_network_new_HC_a_odd = mean_pressure_inlet_new_HC_a_odd - mean_pressure_outlet_new_HC_a_odd
    maximum_pressure_network_new_HC_a_odd = sf_a_odd['pore.pressure'].max()
    
    print('\n ---------------------- ODD NETWORK Larachi physics---------------------------\n')
    # ################### Outputting Stokes flow algorithm data #####################
    print(f'Pressure_drop_single_network new hyd. cond. net_c_odd: {Pressure_drop_network_new_HC_c_odd*1e-5 :.5f} bar')
    print(f'Pressure_drop_total_network new hyd. cond. net_c_odd: {Pressure_drop_network_new_HC_c_odd*1e-5*17 :.5f} bar')
    print(f'maximum_pressure_network new hyd. cond. net_c_odd: {maximum_pressure_network_new_HC_c_odd*1e-5 :.5f} bar')
    print(f'mean_pressure_inlet new hyd. cond. net_c_odd: {mean_pressure_inlet_new_HC_c_odd*1e-5 :.5f} bar')
    print(f'mean_pressure_outlet new hyd. cond. net_c_odd: {mean_pressure_outlet_new_HC_c_odd*1e-5 :.5f} bar')
    
    print(f'Pressure_drop_single_network new hyd. cond. net_a_odd: {Pressure_drop_network_new_HC_a_odd*1e-5 :.5f} bar')
    print(f'Pressure_drop_total_network new hyd. cond. net_a_odd: {Pressure_drop_network_new_HC_a_odd*1e-5*17 :.5f} bar')
    print(f'maximum_pressure_network new hyd. cond. net_a_odd: {maximum_pressure_network_new_HC_a_odd*1e-5 :.5f} bar')
    print(f'mean_pressure_inlet new hyd. cond. net_a_odd: {mean_pressure_inlet_new_HC_a_odd*1e-5 :.5f} bar')
    print(f'mean_pressure_outlet new hyd. cond. net_a_odd: {mean_pressure_outlet_new_HC_a_odd*1e-5 :.5f} bar')
        
    '''------- Plotting and exporting fluid field for odd networks ---------'''
    # Plotting histogram of:
    #     1. Hydraulic conductance 
    #     2. Flow rate
    #     3. Throat velocity
    #     4. Pore velocity
    #  to compare OpenPNM and Larachi fluid field
        
    # Plotting hydraulic conductance histrogram
    plt.figure(0)
    bins = np.linspace(-.25e-11, .75e-11, 200)
    plt.hist(x=Hydraulic_conductance_OpenPNM_cat[throats_internal_c], bins=bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. odd.')
    plt.hist(x=catholyte['throat.hydraulic_conductance'][throats_internal_c], bins=bins, alpha=0.6, rwidth=0.85, label='Larachi et al. physics cat. odd.')
    plt.hist(x=Hydraulic_conductance_OpenPNM_ano[throats_internal_a], bins=bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics ano. odd.')
    plt.hist(x=anolyte['throat.hydraulic_conductance'][throats_internal_a], bins=bins, alpha=0.6, rwidth=0.85, label='Larachi et al. physics ano. odd.')    
    plt.title(label = 'Throat hydraulic conductance Histogram (Internal)')
    plt.xlabel('Throat hydraulic conductance')
    plt.ylabel('# Throats')    
    plt.legend()
    
    # PLotting flow rate histogram
    plt.figure(1)
    # flow_rate_bins = np.linspace(0, 1e-10, 150)   # 10 cms
    flow_rate_bins = np.linspace(0, .3e-10, 150)   # 1.5 cms
    plt.hist(x=Q_throats_OpenPNM_physics_c_odd[throats_internal_c], bins=flow_rate_bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. odd.')
    plt.hist(x=Q_throats_new_HC_c_odd[throats_internal_c], bins=flow_rate_bins, alpha=0.6, rwidth=0.85, label='Larachi et al. physics cat. odd.')
    plt.hist(x=Q_throats_OpenPNM_physics_a_odd[throats_internal_a], bins=flow_rate_bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics ano. odd.')
    plt.hist(x=Q_throats_new_HC_a_odd[throats_internal_a], bins=flow_rate_bins, alpha=0.6, rwidth=0.85, label='Larachi et al. physics ano. odd.')
    plt.title(label = 'Throat absolute flow rate Histogram (Internal)')
    plt.xlabel('Absolute flow rate [m-3/s]')
    plt.ylabel('# Throats')    
    plt.legend()
    
    # Plotting Throat velocity histogram
    plt.figure(2)
    # bins_throat_velocity = np.linspace(0, 1.0 ,200)
    bins_throat_velocity = np.linspace(0, .25 ,200)
    plt.hist(x=net_c['throat.absolute_velocity_OpenPNM_odd_network'][throats_internal_c], bins=bins_throat_velocity, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. odd.')
    plt.hist(x=net_c['throat.absolute_velocity_odd_network'][throats_internal_c], bins=bins_throat_velocity, alpha=0.6, rwidth=0.85, label='Larachi et al. physics cat. odd.')
    plt.hist(x=net_a['throat.absolute_velocity_OpenPNM_odd_network'][throats_internal_a], bins=bins_throat_velocity, alpha=0.6, rwidth=0.85, label='OpenPNM physics ano. odd.')
    plt.hist(x=net_a['throat.absolute_velocity_odd_network'][throats_internal_a], bins=bins_throat_velocity, alpha=0.6, rwidth=0.85, label='Larachi et al. physics ano. odd.')
    plt.title(label = 'Throat absolute throat velocity Histogram (Internal)')
    plt.xlabel('Throat absolute throat velocity [m/s]')
    plt.ylabel('# Throats')    
    plt.legend()
    
    # Plotting pore velocity histogram    
    plt.figure(3)
    bins_pore_velocity = np.linspace(0, 0.40 ,50)
    plt.hist(x=net_c['pore.velocity_magnitude_OpenPNM_odd_network'][net_c.pores('internal')], bins=bins_pore_velocity, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. odd.')
    plt.hist(x=net_c['pore.velocity_magnitude_odd_network'][net_c.pores('internal')], bins=bins_pore_velocity, alpha=0.6, rwidth=0.85, label='Larachi et al. physics cat. odd.')
    plt.hist(x=net_a['pore.velocity_magnitude_OpenPNM_odd_network'][net_c.pores('internal')], bins=bins_pore_velocity, alpha=0.6, rwidth=0.85, label='OpenPNM physics ano. odd.')
    plt.hist(x=net_a['pore.velocity_magnitude_odd_network'][net_c.pores('internal')], bins=bins_pore_velocity, alpha=0.6, rwidth=0.85, label='Larachi et al. physics ano. odd.')
    plt.title(label = 'Pore absolute velocity magnitude Histogram (Internal)')
    plt.xlabel('Pore absolute throat velocity [m/s]')
    plt.ylabel('# Pores')    
    plt.legend()
        
    # Outputting the stokesflow results of OpenPNM and Larachi for the odd networks
    # This is done via a network and its associated phase object
    Folder = '\\Input_Fluid_field\\'  # Folder where the phase is saved containing the fluid field
    project_sf_catholyte_odd_export = op.network.Network()
    project_sf_anolyte_odd_export = op.network.Network()
    sf_catholyte_odd = net_c.project['stokes_01']
    sf_anolyte_odd = net_a.project['stokes_01']
    cf_pres_fit.Export_stokesflow_algo(sf_catholyte_odd, project_sf_catholyte_odd_export,1)
    cf_pres_fit.Export_stokesflow_algo(sf_anolyte_odd, project_sf_anolyte_odd_export,1)
      
    file_name_sf_c_odd = 'SF_' + str(v_in_c*100).replace(".", "_") + 'cms_catholyte_odd_network'
    file_name_sf_a_odd = 'SF_' + str(v_in_c*100).replace(".", "_") + 'cms_anolyte_odd_network'
    path_sf_c_odd = cwd + Folder + file_name_sf_c_odd
    path_sf_a_odd = cwd + Folder + file_name_sf_a_odd
    
    op.Workspace().save_project(project = project_sf_catholyte_odd_export.project,          
                                filename = path_sf_c_odd)
    op.Workspace().save_project(project = project_sf_anolyte_odd_export.project,          
                                filename = path_sf_a_odd)
    
    # output to Paraview VTK file for visualization of results in Paraview
    pn_cat_vtk = op.network.Network()
    pn_cat_vtk.update(net_c)
    catholyte_vtk = op.phase.Phase(network=pn_cat_vtk, name='catholyte_vtk')
    catholyte_vtk.update(catholyte)
    op.io.project_to_vtk(project = pn_cat_vtk.project, filename = path_sf_c_odd)

    pn_an_vtk = op.network.Network()
    pn_an_vtk.update(net_a)
    anolyte_vtk = op.phase.Phase(network=pn_an_vtk, name='anolyte_vtk')
    anolyte_vtk.update(anolyte)
    op.io.project_to_vtk(project = pn_an_vtk.project, filename = path_sf_a_odd)
    
    '''------------------------ Running even networks ----------------------'''
    # Regenerate OpenPNM hydraulic conductance 
    catholyte.regenerate_models()
    anolyte.regenerate_models()
    
    ################# Compute stokesflow for OpenPNM physics ##################
    # Assign value of n, m, Gamma and Init to the network 
    net_c['pore.Fitting_parameter_n'] = n_factor_fit
    net_c['pore.Fitting_parameter_m'] = m_factor_fit
    net_c['pore.Fitting_parameter_Gamma'] = Gamma
    net_c['pore.parameter_Init'] = 0      # Indicate we work with hydraulic conductance of OpenPNM    
    
    net_a['pore.Fitting_parameter_n'] = n_factor_fit
    net_a['pore.Fitting_parameter_m'] = m_factor_fit
    net_a['pore.Fitting_parameter_Gamma'] = Gamma
    net_a['pore.parameter_Init'] = 0      # Indicate we work with hydraulic conductance of OpenPNM    
    
    # Run electrolyte transport algorithms for even networks
    sf_c_even.run()
    sf_a_even.run()
       
    # Find the throat flowrates and mean flowrate of odd flow
    catholyte.update(sf_c_even.soln)
    anolyte.update(sf_a_even.soln)
    
    # Find the throat flowrates and velocity, and the pore velocity of even flow networks and save them to phase
    Q_throats_OpenPNM_physics_c_even, Abs_pres_diff_OpenPNM_physics_c_even = cf_pres_fit.Throat_flowrate_total_hydraulic_conductance(net_c, 
                                                                                                                                     throat_area='throat.area',
                                                                                                                                     pore_diameter='pore.diameter',
                                                                                                                                     throat_diameter='throat.diameter',
                                                                                                                                     Hydraulic_conductance='throat.hydraulic_conductance')
    Q_throats_OpenPNM_physics_c_even_mean = Q_throats_OpenPNM_physics_c_even.mean()
    cf_pres_fit.Throat_pore_velocity_extraction(net_c, Q_throats_OpenPNM_physics_c_even, 0, 0, sf_c_even)
    Q_throats_OpenPNM_physics_a_even, Abs_pres_diff_OpenPNM_physics_a_even = cf_pres_fit.Throat_flowrate_total_hydraulic_conductance(net_a, 
                                                                                                                                     throat_area='throat.area',
                                                                                                                                     pore_diameter='pore.diameter',
                                                                                                                                     throat_diameter='throat.diameter',
                                                                                                                                     Hydraulic_conductance='throat.hydraulic_conductance')
    Q_throats_OpenPNM_physics_a_even_mean = Q_throats_OpenPNM_physics_a_even.mean()
    cf_pres_fit.Throat_pore_velocity_extraction(net_a, Q_throats_OpenPNM_physics_a_even, 0, 0, sf_a_even)
    
    sf_c_even['throat.Hydraulic_conductance_OpenPNM'] = Hydraulic_conductance_OpenPNM_cat.copy() 
    sf_a_even['throat.Hydraulic_conductance_OpenPNM'] = Hydraulic_conductance_OpenPNM_ano.copy()  
    
    ######## Compute pressure drop over network with OpenPNM physics ##########        
    # Catholyte even network
    Inlet_surface_pores_c_even = net_c.find_connected_pores(throats = inlet_throats_c_even)[:,0]
    if Flow_field == 0:
        outlet_throats_c_even = net_c.find_neighbor_throats(pores = net_c.pores('flow_inlet'))
    if Flow_field == 1:
        outlet_throats_c_even = net_c.find_neighbor_throats(pores = net_c.pores('flow_outlet'))
    Outlet_surface_pores_c_even = net_c.find_connected_pores(throats = outlet_throats_c_even)[:,0]
    # Mean pressure inlet and outlet
    mean_pressure_inlet_c_even = sf_c_even['pore.pressure'][Inlet_surface_pores_c_even].mean()
    mean_pressure_outlet_c_even = sf_c_even['pore.pressure'][Outlet_surface_pores_c_even].mean()
    Pressure_drop_network_c_even = mean_pressure_inlet_c_even - mean_pressure_outlet_c_even
    maximum_pressure_network_c_even = sf_c_even['pore.pressure'].max()
    
    # Anolyte odd network
    Inlet_surface_pores_a_even = net_a.find_connected_pores(throats = inlet_throats_a_even)[:,0]  
    if Flow_field == 0:
        outlet_throats_a_even = net_a.find_neighbor_throats(pores = net_a.pores('flow_inlet'))
    if Flow_field == 1:
        outlet_throats_a_even = net_a.find_neighbor_throats(pores = net_a.pores('flow_outlet'))
    Outlet_surface_pores_a_even = net_a.find_connected_pores(throats = outlet_throats_a_even)[:,0]
    # Mean pressure inlet and outlet
    mean_pressure_inlet_a_even = sf_a_even['pore.pressure'][Inlet_surface_pores_a_even].mean()
    mean_pressure_outlet_a_even = sf_a_even['pore.pressure'][Outlet_surface_pores_a_even].mean()
    Pressure_drop_network_a_even = mean_pressure_inlet_a_even - mean_pressure_outlet_a_even
    maximum_pressure_network_a_even = sf_a_even['pore.pressure'].max()
    
    print('\n ---------------------- Even NETWORK OpenPNM physics---------------------------\n')
    ################### Outputting Stokes flow algorithm data #####################
    print(f'Pressure_drop_single_network_net_c_even: {Pressure_drop_network_c_even*1e-5 :.5f} bar')
    print(f'Pressure_drop_total_network_net_c_even: {Pressure_drop_network_c_even*1e-5*17 :.5f} bar')
    print(f'maximum_pressure_network_net_c_even: {maximum_pressure_network_c_even*1e-5 :.5f} bar')
    print(f'mean_pressure_inlet_net_c_even: {mean_pressure_inlet_c_even*1e-5 :.5f} bar')
    print(f'mean_pressure_outlet_net_c_even: {mean_pressure_outlet_c_even*1e-5 :.5f} bar')
    Init_guess_c_even = sf_c_even['pore.pressure'].copy() + 1e-3
    Init_guess_c_even[np.where(Init_guess_c_even == 1e-3)[0]] = 0    
    
    print(f'Pressure_drop_single_network_net_a_even: {mean_pressure_inlet_a_even*1e-5 :.5f} bar')
    print(f'Pressure_drop_total_network_net_a_even: {mean_pressure_inlet_a_even*1e-5*17 :.5f} bar')
    print(f'maximum_pressure_network_net_a_even: {maximum_pressure_network_a_even*1e-5 :.5f} bar')
    print(f'mean_pressure_inlet_net_a_even: {mean_pressure_inlet_a_even*1e-5 :.5f} bar')
    print(f'mean_pressure_outlet_net_a_even: {mean_pressure_outlet_a_even*1e-5 :.5f} bar')
    Init_guess_a_even = sf_a_even['pore.pressure'].copy() + 1e-3
    Init_guess_a_even[np.where(Init_guess_a_even == 1e-3)[0]] = 0   
    
    """--- Setting and solving Stokes-flow algorithm with new Hydraulic model---"""
    # Run stokesflow algorithm with the new hydraulic model
    net_c['pore.parameter_Init'] = 1            # Indicate we work with hydraulc conductance of Larachi et al.
    # sf_c_even.run(x0=Init_guess_c_even)       # Compute with an intitial pressure guess based on the converged results from OpenPNM physics stokesflow
    sf_c_even.run()                             # Compute with an initial pressure guess of zero
    
    net_a['pore.parameter_Init'] = 1            # Indicate we work with hydraulc conductance of Larachi et al.
    # sf_a_even.run(x0=Init_guess_a_even)       # Compute with an intitial pressure guess based on the converged results from OpenPNM physics stokesflow
    sf_a_even.run()                             # Compute with an initial pressure guess of zero

    # Find the throat flowrates and velocity, and the pore velocity of odd flow networks and save them to phase
    catholyte.update(sf_c_even.soln)
    anolyte.update(sf_a_even.soln)
    # Find the throat flowrates and mean flowrate
    Q_throats_new_HC_c_even, Abs_pres_diff_new_HC_c_even = cf_pres_fit.Throat_flowrate_total_hydraulic_conductance(net_c, 
                                                                                       throat_area='throat.area',
                                                                                       pore_diameter='pore.diameter',
                                                                                       throat_diameter='throat.diameter',
                                                                                       Hydraulic_conductance='throat.hydraulic_conductance')
    
    Q_throats_new_HC_c_even_mean = Q_throats_new_HC_c_even.mean()
    cf_pres_fit.Throat_pore_velocity_extraction(net_c, Q_throats_new_HC_c_even, 1, 0, sf_c_even)
    Q_throats_new_HC_a_even, Abs_pres_diff_new_HC_a_even = cf_pres_fit.Throat_flowrate_total_hydraulic_conductance(net_a, 
                                                                                       throat_area='throat.area',
                                                                                       pore_diameter='pore.diameter',
                                                                                       throat_diameter='throat.diameter',
                                                                                       Hydraulic_conductance='throat.hydraulic_conductance')
    
    Q_throats_new_HC_a_even_mean = Q_throats_new_HC_a_even.mean()
    cf_pres_fit.Throat_pore_velocity_extraction(net_a, Q_throats_new_HC_a_even, 1, 0, sf_a_even)    
    
    sf_c_even['throat.Hydraulic_conductance_Larachi'] = catholyte['throat.hydraulic_conductance'].copy()
    sf_a_even['throat.Hydraulic_conductance_Larachi'] = anolyte['throat.hydraulic_conductance'].copy()  
    
    ############## Pressure drop over network Larachi Physics #################
    # Catholyte odd network
    Inlet_surface_pores_c_even = net_c.find_connected_pores(throats = inlet_throats_c_even)[:,0]   
    if Flow_field == 0:
        outlet_throats_c_even = net_c.find_neighbor_throats(pores = net_c.pores('flow_inlet'))
    if Flow_field == 1:
        outlet_throats_c_even = net_c.find_neighbor_throats(pores = net_c.pores('flow_outlet'))
    Outlet_surface_pores_c_even = net_c.find_connected_pores(throats = outlet_throats_c_even)[:,0]
    # Mean pressure inlet and outlet
    mean_pressure_inlet_new_HC_c_even = sf_c_even['pore.pressure'][Inlet_surface_pores_c_even].mean()
    mean_pressure_outlet_new_HC_c_even = sf_c_even['pore.pressure'][Outlet_surface_pores_c_even].mean()
    Pressure_drop_network_new_HC_c_even = mean_pressure_inlet_new_HC_c_even - mean_pressure_outlet_new_HC_c_even
    maximum_pressure_network_new_HC_c_even = sf_c_even['pore.pressure'].max()
    
    # Anolyte odd network
    Inlet_surface_pores_a_even = net_a.find_connected_pores(throats = inlet_throats_a_even)[:,0]    
    if Flow_field == 0:
        outlet_throats_a_even = net_a.find_neighbor_throats(pores = net_a.pores('flow_inlet'))
    if Flow_field == 1:
        outlet_throats_a_even = net_a.find_neighbor_throats(pores = net_a.pores('flow_outlet'))
    Outlet_surface_pores_a_even = net_a.find_connected_pores(throats = outlet_throats_a_even)[:,0]    
    # Mean pressure inlet and outlet
    mean_pressure_inlet_new_HC_a_even = sf_a_even['pore.pressure'][Inlet_surface_pores_a_even].mean()
    mean_pressure_outlet_new_HC_a_even = sf_a_even['pore.pressure'][Outlet_surface_pores_a_even].mean()
    Pressure_drop_network_new_HC_a_even = mean_pressure_inlet_new_HC_a_even - mean_pressure_outlet_new_HC_a_even
    maximum_pressure_network_new_HC_a_even = sf_a_even['pore.pressure'].max()
    
    print('\n ---------------------- Even NETWORK Larachi physics---------------------------\n')
    # ################### Outputting Stokes flow algorithm data #####################
    print(f'Pressure_drop_single_network new hyd. cond. net_c_even: {Pressure_drop_network_new_HC_c_even*1e-5 :.5f} bar')
    print(f'Pressure_drop_total_network new hyd. cond. net_c_even: {Pressure_drop_network_new_HC_c_even*1e-5*17 :.5f} bar')
    print(f'maximum_pressure_network new hyd. cond. net_c_even: {maximum_pressure_network_new_HC_c_even*1e-5 :.5f} bar')
    print(f'mean_pressure_inlet new hyd. cond. net_c_even: {mean_pressure_inlet_new_HC_c_even*1e-5 :.5f} bar')
    print(f'mean_pressure_outlet new hyd. cond. net_c_even: {mean_pressure_outlet_new_HC_c_even*1e-5 :.5f} bar')
    
    print(f'Pressure_drop_single_network new hyd. cond. net_a_even: {Pressure_drop_network_new_HC_a_even*1e-5 :.5f} bar')
    print(f'Pressure_drop_total_network new hyd. cond. net_a_even: {Pressure_drop_network_new_HC_a_even*1e-5*17 :.5f} bar')
    print(f'maximum_pressure_network new hyd. cond. net_a_even: {maximum_pressure_network_new_HC_a_even*1e-5 :.5f} bar')
    print(f'mean_pressure_inlet new hyd. cond. net_a_even: {mean_pressure_inlet_new_HC_a_even*1e-5 :.5f} bar')
    print(f'mean_pressure_outlet new hyd. cond. net_a_even: {mean_pressure_outlet_new_HC_a_even*1e-5 :.5f} bar')
    
    '''------- Plotting and exporting fluid field for even networks --------'''
    # Plotting histogram of:
    #     1. Hydraulic conductance 
    #     2. Flow rate
    #     3. Throat velocity
    #     4. Pore velocity
    #  to compare OpenPNM and Larachi fluid field

    # Plotting hydraulic conductance histrogram
    plt.figure(4)
    bins = np.linspace(-.25e-11, .75e-11, 200)
    plt.hist(x=Hydraulic_conductance_OpenPNM_cat[throats_internal_c], bins=bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. even.')
    plt.hist(x=catholyte['throat.hydraulic_conductance'][throats_internal_c], bins=bins, alpha=0.6, rwidth=0.85, label='Larachi et al. physics cat. even.')
    plt.hist(x=Hydraulic_conductance_OpenPNM_ano[throats_internal_a], bins=bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics ano. even.')
    plt.hist(x=anolyte['throat.hydraulic_conductance'][throats_internal_a], bins=bins, alpha=0.6, rwidth=0.85, label='Larachi et al. physics ano. even.')    
    plt.title(label = 'Throat hydraulic conductance Histogram (Internal)')
    plt.xlabel('Throat hydraulic conductance')
    plt.ylabel('# Throats')    
    plt.legend()
    
    # PLotting flow rate histogram
    plt.figure(5)
    # flow_rate_bins = np.linspace(0, 1e-10, 150)   # 10 cms
    flow_rate_bins = np.linspace(0, .3e-10, 150)   # 1.5 cms
    plt.hist(x=Q_throats_OpenPNM_physics_c_even[throats_internal_c], bins=flow_rate_bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. even.')
    plt.hist(x=Q_throats_new_HC_c_even[throats_internal_c], bins=flow_rate_bins, alpha=0.6, rwidth=0.85, label='Larachi et al. physics cat. even.')
    plt.hist(x=Q_throats_OpenPNM_physics_a_even[throats_internal_a], bins=flow_rate_bins, alpha=0.6, rwidth=0.85, label='OpenPNM physics ano. even.')
    plt.hist(x=Q_throats_new_HC_a_even[throats_internal_a], bins=flow_rate_bins, alpha=0.6, rwidth=0.85, label='Larachi et al. physics ano. even.')
    plt.title(label = 'Throat absolute flow rate Histogram (Internal)')
    plt.xlabel('Absolute flow rate [m-3/s]')
    plt.ylabel('# Throats')    
    plt.legend()
    
    # Plotting Throat velocity histogram
    plt.figure(6)
    # bins_throat_velocity = np.linspace(0, 1.0 ,200)
    bins_throat_velocity = np.linspace(0, .25 ,200)
    plt.hist(x=net_c['throat.absolute_velocity_OpenPNM_even_network'][throats_internal_c], bins=bins_throat_velocity, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. even.')
    plt.hist(x=net_c['throat.absolute_velocity_even_network'][throats_internal_c], bins=bins_throat_velocity, alpha=0.6, rwidth=0.85, label='Larachi et al. physics cat. even.')
    plt.hist(x=net_a['throat.absolute_velocity_OpenPNM_even_network'][throats_internal_a], bins=bins_throat_velocity, alpha=0.6, rwidth=0.85, label='OpenPNM physics ano. even.')
    plt.hist(x=net_a['throat.absolute_velocity_even_network'][throats_internal_a], bins=bins_throat_velocity, alpha=0.6, rwidth=0.85, label='Larachi et al. physics ano. even.')
    plt.title(label = 'Throat absolute throat velocity Histogram (Internal)')
    plt.xlabel('Throat absolute throat velocity [m/s]')
    plt.ylabel('# Throats')    
    plt.legend()
    
    # Plotting pore velocity histogram    
    plt.figure(7)
    bins_pore_velocity = np.linspace(0, 0.40 ,50)
    plt.hist(x=net_c['pore.velocity_magnitude_OpenPNM_even_network'][net_c.pores('internal')], bins=bins_pore_velocity, alpha=0.6, rwidth=0.85, label='OpenPNM physics cat. even.')
    plt.hist(x=net_c['pore.velocity_magnitude_even_network'][net_c.pores('internal')], bins=bins_pore_velocity, alpha=0.6, rwidth=0.85, label='Larachi et al. physics cat. even.')
    plt.hist(x=net_a['pore.velocity_magnitude_OpenPNM_even_network'][net_c.pores('internal')], bins=bins_pore_velocity, alpha=0.6, rwidth=0.85, label='OpenPNM physics ano. even.')
    plt.hist(x=net_a['pore.velocity_magnitude_even_network'][net_c.pores('internal')], bins=bins_pore_velocity, alpha=0.6, rwidth=0.85, label='Larachi et al. physics ano. even.')
    plt.title(label = 'Pore absolute velocity magnitude Histogram (Internal)')
    plt.xlabel('Pore absolute throat velocity [m/s]')
    plt.ylabel('# Pores')    
    plt.legend()
    
    # Outputting the stokesflow results of OpenPNM and Larachi for the even networks
    # This is done via a network and its associated phase object
    project_sf_catholyte_even_export = op.network.Network()
    project_sf_anolyte_even_export = op.network.Network()
    sf_catholyte_even = net_c.project['stokes_02']
    sf_anolyte_even = net_a.project['stokes_02']
    cf_pres_fit.Export_stokesflow_algo(sf_catholyte_even, project_sf_catholyte_even_export,0)
    cf_pres_fit.Export_stokesflow_algo(sf_anolyte_even, project_sf_anolyte_even_export,0)
        
    Folder = '\\Input_Fluid_field\\'
    file_name_sf_c_even = 'SF_' + str(v_in_c*100).replace(".", "_") + 'cms_catholyte_even_network'
    file_name_sf_a_even = 'SF_' + str(v_in_c*100).replace(".", "_") + 'cms_anolyte_even_network'
    path_sf_c_even = cwd + Folder + file_name_sf_c_even
    path_sf_a_even = cwd + Folder + file_name_sf_a_even
    
    op.Workspace().save_project(project = project_sf_catholyte_even_export.project,          
                                filename = path_sf_c_even)
    op.Workspace().save_project(project = project_sf_anolyte_even_export.project,          
                                filename = path_sf_a_even)
    
    # output to Paraview VTK file for visualization of results in Paraview
    pn_cat_vtk_even = op.network.Network()
    pn_cat_vtk_even.update(net_c)
    catholyte_vtk_even = op.phase.Phase(network=pn_cat_vtk_even, name='catholyte_vtk_even')
    catholyte_vtk_even.update(catholyte)
    op.io.project_to_vtk(project = pn_cat_vtk_even.project, filename = path_sf_c_even)

    pn_an_vtk_even = op.network.Network()
    pn_an_vtk_even.update(net_a)
    anolyte_vtk_even = op.phase.Phase(network=pn_an_vtk_even, name='anolyte_vtk_even')
    anolyte_vtk_even.update(anolyte)
    op.io.project_to_vtk(project = pn_an_vtk_even.project, filename = path_sf_a_even)
    
    '''----------- Exporting the altered odd and even network --------------'''
    Folder = '\\Input_networks\\'
    file_name_net_c = network_name + '_catholyte_network_throat_diameter_scaled'
    file_name_net_a = network_name + '_anolyte_network_throat_diameter_scaled'
    path_net_c = cwd + Folder + file_name_net_c
    path_net_a = cwd + Folder + file_name_net_a
    
    Export_network_cat = op.network.Network()
    Export_network_cat.update(net_c)
    
    Export_network_ano = op.network.Network()
    Export_network_ano.update(net_a)
    
    op.Workspace().save_project(project = Export_network_cat.project,          
                                filename = path_net_c)
    op.Workspace().save_project(project = Export_network_ano.project,          
                                filename = path_net_a)