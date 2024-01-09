import numpy as np
import openpnm as op
import matplotlib.pyplot as plt
import numpy as np

def validate_face_labels(net):
    r"""
    This function validates if the network has the correct face labels.    
    Usage example: validate_face_labels(net=net_c)
    """
    # Update face labels if the network is extracted using the SNOW algorithm:
    if 'pore.left' and 'pore.front' and 'pore.top' in net:
        print('Network contains the correct labels, please make sure that the'
              ' labels are consistent with:\n'
                  '\'left\' and \'right\' - x-dimension\n'
                  '\'front\' and \'back\' - y-dimension\n'
                  '\'top\' and \'bottom\' - z-dimension\n')
    else: 
        raise KeyError('Network does not contain the correct boundary pore labels.\n'
                       'Please assign the following network labels:\n'
                       '\'left\' and \'right\' - x-dimension\n'
                       '\'front\' and \'back\' - y-dimension\n'
                       '\'top\' and \'bottom\' - z-dimension\n')

def assign_boundary_pores(net, W_dim, L_dim):
    r"""
    This function assigns the right labels to the boundary pores for the given network.
    
    The SNOW algorithm adds the labels 'left' and 'right' to the x-dimension
                                       'front' and 'back' to the y-dimension
                                       'top' and 'bottom' to the z-dimension
    Because we invert the network in the thickness, we can assign the same boundary
    pores for the membrane in both the anode and the cathode in this step.
    """
    
    if W_dim == 0:
        net['pore.membrane'] = net['pore.left']
        net['pore.current_collector'] = net['pore.right']
    elif W_dim == 1:
        net['pore.membrane'] = net['pore.front']
        net['pore.current_collector'] = net['pore.back']
    elif W_dim == 2:
        net['pore.membrane'] = net['pore.top']
        net['pore.current_collector'] =  net['pore.bottom']
    
    if L_dim == 0:
        net['pore.flow_inlet'] = net['pore.left']
        net['pore.flow_outlet'] = net['pore.right']
    elif L_dim == 1:
        net['pore.flow_inlet'] = net['pore.front']
        net['pore.flow_outlet'] = net['pore.back']
    elif L_dim == 2:
        net['pore.flow_inlet'] = net['pore.bottom']
        net['pore.flow_outlet'] = net['pore.top']


def assign_boundary_pores_IDFF(net, W_dim, L_dim, H_dim, H, H_min, Inlet_channel):
    r"""
    This function assigns the right labels to the boundary pores for the given 
    network with an interdigitated flow field
    
    The SNOW algorithm adds the labels 'left' and 'right' to the x-dimension
                                       'front' and 'back' to the y-dimension
                                       'top' and 'bottom' to the z-dimension
    Because we invert the network in the thickness, we can assign the same boundary
    pores for the membrane in both the anode and the cathode in this step. 
    
    The inlet and outlet is assigned with boundary conditions or via flow field (FF) pores:
        
    1. Inlet via boundary conditions:
        The network will be sliced in 3 parts via the height along the length and 
        thickness of the electrode:     Inlet channel - rib - Outlet channel
        This will be done by assigning masks based on the height (pore y-coordinate)
        domain and the by selecting all pores on the left boundaries.
            Inlet channel:  H_min < Pore y-coordinate <= Quarter_heigth + H_min
            rib:            H_min + Quarter_heigth < Pore y-coordinate <= H_min + 3 * Quarter_heigth
            Outlet channel: H_min + 3 * Quarter_heigth < Pore y-coordinate <= H_min + 4 * Quarter_heigth
        The boundary pores of the inlet channel will be assigned the label: 'pore.flow_inlet'.
        The boundary pores of the outlet channel will be assigned the label: 'pore.flow_outlet'
        
    2. Inlet via FF pores:
        Bottom boundary pores of the inlet flow field are labeled as inlet
        Top boundary pores of the outlet flow field are labeled as outlet
        
    The boundary pores of the rib will be assigned the label: 'pore.current_collector'
    """
    
    # In case of using Boundary conditions
    if Inlet_channel == 0:              
        total_heigth = H - H_min
        Quarter_heigth = total_heigth/4
        
        Heigth_0 = H_min
        Heigth_1 = H_min + Quarter_heigth
        Heigth_2 = H_min + 3 * Quarter_heigth
        Heigth_3 = H_min + 4* Quarter_heigth
    
        # Masking: note "xxx_mask_boolean" and "xxx_mask" contain ALL pores at this height -> we need to
        # select only those present at the left surface of the electrode (this excludes the other boundary pores like bottom and top. 
        # The pores at the left AND at the specefic heights are then selected in 'xxx_mask_left':
        Inlet_mask_boolean = np.logical_and(net['pore.coords'][net.pores(), H_dim] >= Heigth_0,
                                    net['pore.coords'][net.pores(), H_dim] <= Heigth_1)
        Inlet_mask = np.where(Inlet_mask_boolean)[0]
        Inlet_heigth_pores = net.pores('all')[Inlet_mask]
        net.set_label(label='height_1', pores = Inlet_heigth_pores)
        Inlet_mask_left = net.pores(['height_1', 'left'], mode='and')
        Inlet_pores = net.pores('all')[Inlet_mask_left]
        net.set_label(label = 'flow_inlet', pores=Inlet_pores)
    
        Rib_mask_boolean = np.logical_and(net['pore.coords'][net.pores(), H_dim] > Heigth_1,
                                    net['pore.coords'][net.pores(), H_dim] <= Heigth_2)
        Rib_mask = np.where(Rib_mask_boolean)[0]
        Rib_heigth_pores = net.pores('all')[Rib_mask]
        net.set_label(label='height_2', pores = Rib_heigth_pores)
        Rib_mask_left = net.pores(['height_2', 'left'], mode='and')
        Rib_pores = net.pores('all')[Rib_mask_left]
        net.set_label(label = 'pore.current_collector', pores=Rib_pores)
    
        Outlet_mask_boolean = np.logical_and(net['pore.coords'][net.pores(), H_dim] > Heigth_2,
                                    net['pore.coords'][net.pores(), H_dim] <= Heigth_3)
        Outlet_mask = np.where(Outlet_mask_boolean)[0]
        Outlet_heigth_pores = net.pores('all')[Outlet_mask]
        net.set_label(label='height_3', pores = Outlet_heigth_pores)
        Outlet_mask_left = net.pores(['height_3', 'left'], mode='and')
        Outlet_pores = net.pores('all')[Outlet_mask_left]
        net.set_label(label = 'pore.flow_outlet', pores=Outlet_pores)
        
    # In case of using Flow field pores
    elif Inlet_channel == 1:           
        net.set_label(label = 'flow_inlet', pores=net.pores('Boundary_FF_inlet_bottom'))
        net.set_label(label = 'pore.flow_outlet', pores=net.pores('Boundary_FF_outlet_top'))
        
    net['pore.membrane'] = net['pore.right']
    
    
def add_throat_surface_area_to_pores(net):
    r"""
    This function updatse the internal surface area of every pore with half the area of 
    the connecting throats to account for surface area of the throats.
    
    Usage example: add_throat_surface_area_to_pores(net=net_c)
    """
    
    net['throat.surface_area'] = op.models.geometry.throat_surface_area.cylinder(net, throat_diameter='throat.diameter', throat_length='throat.length')

    for pore in net.Ps:
        connected_throats = net.find_neighbor_throats(pores=pore)
        net['pore.surface_area'][pore] += 1 / 2 * np.sum(net['throat.surface_area'][connected_throats])

        
# Formulation of the Butler-Volmer equation
def bv_rate_constant_oc_c(c, eta, Ai_c, rp_c, j0, km_factor, km_pore, cnst):
    r"""
    Calculates A2 in linear kinetics for OhmicConduction algorithm in the cathode.
    """
    
    c1 = j0 * Ai_c * c / cnst['conc_ref_c']
    c2 = j0 / (cnst['F'] * cnst['conc_ref_c'] * km_factor * km_pore)
    c3 = j0 / (cnst['F'] * cnst['conc_ref_c'] * km_factor * km_pore)
    # c2 = j0 / (cnst['F'] * cnst['conc_ref_c'] * km_factor * cnst['D_c'] / rp_c)
    # c3 = j0 / (cnst['F'] * cnst['conc_ref_c'] * km_factor * cnst['D_a'] / rp_c)
    arg1 = cnst['alpha_a_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
    arg2 = -cnst['alpha_c_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
    
    return (c1 * (np.exp(arg1) - np.exp(arg2)))/(1 + c2 * np.exp(arg1) + c3 * np.exp(arg2))


def bv_rate_constant_oc_a(c, eta, Ai_a, rp_a, j0, km_factor, km_pore, cnst):
    r"""
    Calculates A2 in linear kinetics for OhmicConduction algorithm in the anode.
    """
    c1 = j0 * Ai_a * c / cnst['conc_ref_a']
    c2 = j0 / (cnst['F'] * cnst['conc_ref_a'] * km_factor * km_pore)
    c3 = j0 / (cnst['F'] * cnst['conc_ref_a'] * km_factor * km_pore)
    # c2 = j0 / (cnst['F'] * cnst['conc_ref_a'] * km_factor * cnst['D_c'] / rp_a)
    # c3 = j0 / (cnst['F'] * cnst['conc_ref_a'] * km_factor * cnst['D_a'] / rp_a)
    arg1 = -cnst['alpha_a_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
    arg2 = cnst['alpha_c_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
    
    return (c1 * (np.exp(arg1) - np.exp(arg2)))/(1 + c2 * np.exp(arg1) + c3 * np.exp(arg2))


def bv_rate_derivative_oc_c(conc_c, eta_c, Ai_c, rp_c, j0, km_factor, km_pore, cnst ):
    c1 = j0 * Ai_c * conc_c / (cnst['conc_ref_c'])
    c2 = j0 / (cnst['F'] * cnst['conc_ref_c'] * km_factor * km_pore)
    c3 = j0 / (cnst['F'] * cnst['conc_ref_c'] * km_factor * km_pore)
    # c2 = j0 / (cnst['F'] * cnst['conc_ref_c'] * km_factor * cnst['D_c'] / rp_c)
    # c3 = j0 / (cnst['F'] * cnst['conc_ref_c'] * km_factor * cnst['D_a'] / rp_c)
    m1 = cnst['alpha_a_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T'])
    m2 = cnst['alpha_c_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T'])
    
    nom1 = -c3*m2*np.exp(m1*eta_c-m2*eta_c)
    nom2 = -c2*m2*np.exp(m1*eta_c-m2*eta_c)
    nom3 = -c2*m1*np.exp(m1*eta_c-m2*eta_c)
    nom4 = -c3*m1*np.exp(m1*eta_c-m2*eta_c)
    nom5 = -m2*np.exp(-m2*eta_c)
    nom6 = -m1*np.exp(m1*eta_c)
    
    den1 = c2*np.exp(m1*eta_c)
    den2 = c3*np.exp(-m2*eta_c)
    
    nom = c1*(nom1+nom2+nom3+nom4+nom5+nom6)
    den = (1+den1+den2)**2
    
    return nom/den
    

def bv_rate_derivative_oc_a(conc_a, eta_a, Ai_a, rp_a, j0, km_factor, km_pore, cnst):
    c1 = j0 * Ai_a * conc_a / cnst['conc_ref_a']
    c2 = j0 / (cnst['F'] * cnst['conc_ref_a'] * km_factor * km_pore)
    c3 = j0 / (cnst['F'] * cnst['conc_ref_a'] * km_factor * km_pore)
    # c2 = j0 / (cnst['F'] * cnst['conc_ref_a'] * km_factor * cnst['D_c'] / rp_a)
    # c3 = j0 / (cnst['F'] * cnst['conc_ref_a'] * km_factor * cnst['D_a'] / rp_a)
    m1 = cnst['alpha_a_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T'])
    m2 = cnst['alpha_c_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T'])
    
    nom1 = -c2*m1*np.exp(m2*eta_a-m1*eta_a)
    nom2 = -c3*m1*np.exp(m2*eta_a-m1*eta_a)
    nom3 = -c3*m2*np.exp(m2*eta_a-m1*eta_a)
    nom4 = -c2*m2*np.exp(m2*eta_a-m1*eta_a)
    nom5 = -m1*np.exp(-m1*eta_a)
    nom6 = -m2*np.exp(m2*eta_a)
    
    den1 = c2*np.exp(-m1*eta_a)
    den2 = c3*np.exp(m2*eta_a)
    
    nom = c1*(nom1+nom2+nom3+nom4+nom5+nom6)
    den = (1+den1+den2)**2

    return nom/den 


def bv_rate_constant_ad_c(eta, Ai_c, rp_c, j0, km_factor, km_pore, cnst):
    r"""
    Calculate A1 in linear kinetics for AdvectionDiffusion algorithm in the cathode.
    """
    # c = 1.0 is a workaround so that A1 is rate "constant" not the actual rate
    return bv_rate_constant_oc_c(c=1, eta=eta, Ai_c=Ai_c, rp_c=rp_c, j0=j0, km_factor=km_factor, km_pore = km_pore, cnst=cnst) / (cnst['F'] * cnst['val_c'])
    

def bv_rate_constant_ad_a(eta, Ai_a, rp_a, j0, km_factor, km_pore, cnst):
    r"""
    Calculates A1 in linear kinetics for AdvectionDiffusion algorithm in the anode.
    """
    # c = 1.0 is a workaround so that A1 is rate "constant" not the actual rate
    return bv_rate_constant_oc_a(c=1, eta=eta, Ai_a=Ai_a, rp_a=rp_a, j0=j0, km_factor=km_factor, km_pore = km_pore, cnst=cnst) / (cnst['F'] * cnst['val_a'])


def rel_error(current_new, current):
    r"""
    Calculates the relative error of the current estimation with the previous estimation."""
    # Check if it's safe to calculate relative error (division by 0)
    if current != 0.0:
        rel_err = abs((current_new - current) / current)
        
    # Solution has converged, but relative error --> infinity (division by 0)
    elif current_new == current == 0.0:
        rel_err = 0.0
        
    # Solution has not converged, set relative error to an arbitrary high value
    else:
        rel_err = 1e3
        
    return rel_err


def find_eta_act(eta, cell, i_actual, Ai, j0, cnst):
    r'''
    find_eta_act is passed onto the least squares optimization to compute the activation 
    overpotential required for a given current. It is used to separate the contribution
    of concentration and activation overpotential on the cell performance. 
    
    Parameters
    ----------
    eta : Guess for the activation overpotential [V].
    cell : cell type (anode or cathode).
    i_actual : current that the pore should generate.
    Ai : Internal surface area of the considered pore.
    p : Internal pore radius of the considered pore.
            
    Returns
    -------
    err : relative error between the calculated current and i_actual
    for a guessed overpotential 
    '''
    
    if cell == 'anode':
        # Compute the current that can be generated without concentration effects
        c1 = j0 * Ai * cnst['conc_in_a'] / cnst['conc_ref_a']
        arg1 = -cnst['alpha_a_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
        arg2 = cnst['alpha_c_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta        
        
        ideal_current = c1 * (np.exp(arg1) - np.exp(arg2))
        # Calculate the error between the ideal current and the actual current
        err = (ideal_current - i_actual) / i_actual          

    elif cell == 'cathode':
        # Compute the current that can be generated without concentration effects
        c1 = j0 * Ai * cnst['conc_in_c'] / cnst['conc_ref_c']
        arg1 = cnst['alpha_a_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
        arg2 = -cnst['alpha_c_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
        
        ideal_current = c1 * (np.exp(arg1) - np.exp(arg2))
        err = (ideal_current - i_actual) / i_actual
                             
    return err


def inlet_pressure(P_in, alg, Q_desired, inlet_throats, inlet_pores, net, phase, phys):
    r''' inlet_pressure is used as a target function to find the required inlet pressure
    in every pore to match the desired inlet flow rate Q_desired'''
    
    alg.set_value_BC(values=P_in, pores=inlet_pores, mode='overwrite') # Dirichlet boundary inlet condition
    alg.run()
    
    phase.update(alg.results())
    
    pore_1 = net['throat.conns'][inlet_throats][:, 0]
    pore_2 = net['throat.conns'][inlet_throats][:, 1]
    
    delta_P = abs(phase['pore.pressure'][pore_1] - phase['pore.pressure'][pore_2])
    Q_tot = (phys['throat.hydraulic_conductance'][inlet_throats] * delta_P).sum()
    
    return (Q_tot - Q_desired)/Q_tot


def inlet_pressure_V3(P_in, alg, Q_desired, inlet_throats, inlet_pores, net, phase, outlet_pores):
    r''' inlet_pressure is used as a target function to find the required inlet pressure
    in every pore to match the desired inlet flow rate Q_desired'''
    
    alg.set_value_BC(pores = net.pores('internal', mode = 'nor'), mode='remove')
    alg.set_value_BC(values = 0.0, pores=outlet_pores, mode='overwrite') # Dirichlet boundary inlet condition
    alg.set_value_BC(values=P_in, pores=inlet_pores, mode='overwrite') # Dirichlet boundary inlet condition
    alg.run()
    
    phase.update(alg.soln)
    
    pore_1 = net['throat.conns'][inlet_throats][:, 0]
    pore_2 = net['throat.conns'][inlet_throats][:, 1]
    
    delta_P = (phase['pore.pressure'][pore_2] - phase['pore.pressure'][pore_1])
    Q_tot = (phase['throat.hydraulic_conductance'][inlet_throats] * delta_P).sum()
    
    return (Q_tot - Q_desired)/Q_tot


def Flowrate(outlet_throats, net, phase):
    r'''
    Function computes the outlet flow rate of the network (per outlet throat), 
    and is used for coupling the networks. 

    Parameters
    ----------
    inlet_throats : throat
        Throats connecting the inlet flow boundary pores to the electrode/flow channel
    outlet_throats : TYPE
        Throats connecting the outlet flow boundary pores to the electrode/flow channel
    net : OpenPNM network
    phase : OpenPNM phase

    Returns
    -------
    Returns the outlet flow rate of the network.
    '''
    
    pore_network = net['throat.conns'][outlet_throats][:, 0]    
    pore_boundary = net['throat.conns'][outlet_throats][:, 1]
    delta_P = abs(phase['pore.pressure'][pore_network] - phase['pore.pressure'][pore_boundary])
    Q_outlet = phase['throat.hydraulic_conductance'][outlet_throats] * delta_P
    
    # If you would like to plot the outlet pores, for verification.
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    # op.visualization.plot_coordinates(net, pore_network, c='green', 
    #                   markersize=10, ax=ax) 
    # op.visualization.plot_coordinates(net, pore_boundary, c='red', 
    #                   markersize=10, ax=ax)
    
    return Q_outlet


def Membrane_conductivity(Electrolyte_conductivity, param):
    r'''
    Computes membrane resisitivity based on the conductivity 
    of the electrolyte with Bruggemans correlation.

    Parameters
    ----------
    Electrolyte_conductivity : The conductivity of the electrolyte [S/m]
        
    param : Dictionairy     

    Returns
    -------
    membrane_resistivity : TYPE
        The resisitivty of the membrane [Ohm m-2]
    '''
    
    Membrane_thickness = param["membrane_thickness"] 
    Membrane_porosity = param["membrane_porosity"] 
    
    Effective_conductivity_membrane = Electrolyte_conductivity * (Membrane_porosity)**(1.5)
    membrane_resistivity = Membrane_thickness / Effective_conductivity_membrane 
     
    return membrane_resistivity


def costum_diffusivity(target , prop):
    Diff = prop
    
    return Diff    