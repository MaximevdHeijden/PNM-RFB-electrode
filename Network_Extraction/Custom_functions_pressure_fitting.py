"""
Script contains functions for:
    1. Computing throat flowrate and pressure drop
    2. Computing Hydraulic conductance 
    3. Extracting pore and throat velocity
    4. Exporting/Importing the stokesflow algorithms
    5. Updating a phase with the stokesflow algorithm
    6. Computing mass transfer coefficients
"""

import numpy as np
import openpnm as op
import sys

def Throat_flowrate_total_hydraulic_conductance(target, 
                                                throat_area='throat.area',
                                                pore_diameter='pore.diameter',
                                                throat_diameter='throat.diameter',
                                                Hydraulic_conductance='throat.hydraulic_conductance'):
    r'''
    Function 
        1. Assigns the contraction and expansion pore per conduit throat and and the direction of the flow
            If there is no flow between the two conduit pores, a value of *-1* is assigned.
        2. Computes absolute pressure difference between the two connected pores per conduit
        3. Computes absolute throat flowrate beween the two connected pores per conduit

    Parameters
    ----------
    target : OpenPNM network
    
    throat_area : string
        Dictionary key of throat cross-sectional area.
        
    pore_diameter : string
        Dictionary key of the pore diameter values.

    throat_diameter : string
        Dictionary key of the throat diameter values.
        
    Hydraulic_conductance : ndarray, float
        The throat hydraulic conductancte 

    Returns
    -------
    Q_throats : ndarray, float
        Absolute flow rate in throat of each conduit [m3/s].
        
    Abs_pressure_difference : ndarray, float
        Absolute pressure drop in throat of each conduit [Pa].
    '''
    
    # Obtain network and corresponding phase
    np.warnings.filterwarnings('ignore', category=RuntimeWarning)
    network = target.project.network
    phase = network.project.phases.copy()[0]
    
    # Conduit throat hydraulic conductance
    throats = network.throats('all')
    cn = network['throat.conns'][throats]
    Hyd_con_throats = phase[Hydraulic_conductance][throats]
    
    # Pressure in each pore of conduit and pressure drop over conduit
    Pressure_pore_1 = phase['pore.pressure'][cn[:, 0]]
    Pressure_pore_2 = phase['pore.pressure'][cn[:, 1]]
    Pres_diff = Pressure_pore_1 - Pressure_pore_2
    
    # Absolute pressure drop and flow rate in conduit
    Abs_pressure_difference = np.abs(Pres_diff)
    Q_throats = Abs_pressure_difference * Hyd_con_throats
    abs_Q_throats = np.abs(Q_throats)
    
    ''' Assigning per conduit throat the contraction and expansion pore, 
        and the direction of the flow. The following methodology is maintained:
            
        Per conduit number the pressure of its connecting pores (termed here 
        Pore 1 and 2) are checked. Their difference is set as: 
            Pres_diff = Pressure_pore_1 - Pressure_pore_2
        
        If conduit Pressure_pore_1 < Pressure_pore_2 -> flow from pore 2 to pore 1
            We assign network['throat.Flow_direction'] = 12 (flow from pore 1 -> pore 2).
            network['throat.Pore_Contraction']  = Pore number 2
            network['throat.Pore_Expansion'][throats_idx] = Pore number 1
            
         If conduit Pressure_pore_2 < Pressure_pore_1 -> flow from pore 1 to pore 2
             We assign network['throat.Flow_direction'] = 21 (flow from pore 2 -> pore 1).
             network['throat.Pore_Contraction']  = Pore number 1
             network['throat.Pore_Expansion']= Pore_2
             
        If there is no pressure difference between them, there will also be no
        flow between pore 1 and 2. Then the network['throat.Pore_Contraction'], 
        ['throat.Pore_Contraction'] and ['throat.Pore_Expansion'] are assigned 
        a value of *-1*. 
        '''

    # Initializing the contraction, expansion and flow direction arrays for each pore-throat-pore element
    network['throat.Pore_Contraction'] = np.nan
    network['throat.Pore_Expansion'] = np.nan
    network['throat.Flow_direction'] = np.nan
    
    # Asign contraction (inflow) and expansion (outflow) pores for each pore-throat-pore element
    # Negative pressure drop
    # If conduit Pressure_pore_1 < Pressure_pore_2 -> flow from pore 2 to pore 1
    idx_neg_pres = np.where(Pres_diff<0)[0]
    Pore_1_idx_neg_pres = cn[idx_neg_pres, 0]
    Pore_2_idx_neg_pres = cn[idx_neg_pres, 1]
    network['throat.Pore_Contraction'][idx_neg_pres] = Pore_2_idx_neg_pres
    network['throat.Pore_Expansion'][idx_neg_pres] = Pore_1_idx_neg_pres
    network['throat.Flow_direction'][idx_neg_pres] = 21        

    # Positive pressure drop
    # If conduit Pressure_pore_1 > Pressure_pore_2 -> flow from pore 1 to pore 2
    idx_pos_pres = np.where(Pres_diff>0)[0]
    Pore_1_idx_pos_pres = cn[idx_pos_pres, 0]
    Pore_2_idx_pos_pres = cn[idx_pos_pres, 1]
    network['throat.Pore_Contraction'][idx_pos_pres] = Pore_1_idx_pos_pres
    network['throat.Pore_Expansion'][idx_pos_pres] = Pore_2_idx_pos_pres
    network['throat.Flow_direction'][idx_pos_pres] = 12 
    
    # No pressure drop
    # If conduit Pressure_pore_1 == Pressure_pore_2 -> No flow
    idx_no_pres = np.where(Pres_diff == 0)[0]
    network['throat.Pore_Contraction'][idx_no_pres] = -1
    network['throat.Pore_Expansion'][idx_no_pres] = -1
    network['throat.Flow_direction'][idx_no_pres] = -1

    # Converting labels from float to integers.
    network['throat.Pore_Contraction'] = network['throat.Pore_Contraction'].astype(int)
    network['throat.Pore_Expansion'] = network['throat.Pore_Expansion'].astype(int)
    network['throat.Flow_direction'] = network['throat.Flow_direction'].astype(int)
    
    return abs_Q_throats, Abs_pressure_difference

    
def Reynolds_throat_total_hydraulic_conductance(target, 
                                                throat_diameter='throat.diameter'):
    r'''
    Function computes the throat reynolds of each conduit following the method
    of  F. Larachi et al. https://doi.org/10.1016/j.cej.2013.11.077 (eqn 13) 
    "X-ray micro-tomography and pore network modeling of single-phase fixed-bed reactors"
    
    Parameters
    ----------
    target : OpenPNM network
        
    parameter_script : dict
        Dictionairy with phase properties

    throat_diameter : ndarray, float
        Extracted throat diameter from X-ray CT image. The default is 'throat.diameter'.

    Returns
    -------
    Re_throat : ndarray, float
        Throat Reynolds.
    '''
    
    np.warnings.filterwarnings('ignore', category=RuntimeWarning)
    network = target.project.network
    phase = network.project.phases.copy()[0]
    
    throats = network.throats('all')
    throat_diam = network[throat_diameter][throats]
    Q_throats, Abs_pressure_difference = Throat_flowrate_total_hydraulic_conductance(network) 
       
    viscosity = phase['pore.viscosity'][0]      # [Pa s]
    density = phase['pore.density'][0]          # [kg m-3]
    
    Re_throat = (4 * density * Q_throats) / (np.pi * viscosity * throat_diam)
    
    return Re_throat

    
def Total_hydraulic_conductance_inv(target,
                                    n,
                                    m,
                                    Gamma
                                    ):   
    r'''
    Computes the constriction dissipation (A), Constriction dissipation (C) and
    Expansion dissipation (E) loss factor according to:
     
    Liu et al - AIChE J., 66 (9) (2020), p. e16258. "A pore network model for 
     calculating pressure drop in packed beds of arbitrary-shaped particles."
    
    Than combines these into a single hydraulic conductance per conduit as given 
    by equation 3.4 in Thesis. Note that the velocity head is neglected in the 
    computation.
    
    Parameters
    ----------
    target : OpenPNM network
        
    n : float
        contraction curvature parameter (fitting factor)
    
    m : float
        expansion curvature parameter (fitting factor)
    
    Gamma : float
        flow pattern constant
    
    Returns
    -------
    Hydraulic_conductance_network : Ndarray, float
        Inverse of constriction dissipation loss factors for each throat.
    '''
    
    # Obtain network and corresponding phase, and set labels 
    np.warnings.filterwarnings('ignore', category=RuntimeWarning)
    network = target.project.network
    phase = network.project.phases.copy()[0]
    throat_diameter = 'throat.diameter'
    throat_cond_length = 'throat.conduit_lengths.throat'    
    
    # Set the laminar constants and obtain the physical properties of the fluid
    C0 = 26                                     # Laminar constant for contraction
    E0 = 27                                     # Laminar constant for expansion
    viscosity = phase['pore.viscosity'][0]      # [Pa s]
    density = phase['pore.density'][0]          # [kg m-3]    
    
    # Obtain the flowrate and the throat Reynolds
    Q_throats, Abs_pressure_difference = Throat_flowrate_total_hydraulic_conductance(network)
    Throat_rey = Reynolds_throat_total_hydraulic_conductance(target, 
                                                             throat_diameter='throat.diameter')    
    
    # Extract throat geometrical and flow directional properties
    throats = network.throats('all')    
    throat_rad = network[throat_diameter][throats]/2
    rad_con = network['pore.diameter'][network['throat.Pore_Contraction']]/2  # Pore radius of the conduit constriction pores of a pore-throat-pore element
    rad_exp = network['pore.diameter'][network['throat.Pore_Expansion']]/2    # Pore radius of the conduit expansion pores of a pore-throat-pore element
    Throat_conduit_length = network[throat_cond_length]                       # Conduit length of the throat-pore-throat element    

    ############ Inverse of constriction dissipation loss factor ##############
    A_inv = (np.pi * throat_rad**4)/(8 * viscosity * Throat_conduit_length)
   
    ####### Inverse of contraction & Expansion dissipation loss factor ########
    # Initialize the inverse of the contraction and expansion dissipation loss factor
    C_inv = np.zeros(shape = [len(network.throats('all'))])
    E_inv = np.zeros(shape = [len(network.throats('all'))])
    
    # Select throats having a non-zero flowrate (avoid infs in computation of C_inv and E_inv)
    Internal_throats = network.throats('internal')
    Non_zero_flow_direction = np.where(network['throat.Flow_direction'] > 0)[0]
    int_throats_nonzero_flow = np.intersect1d(Internal_throats, Non_zero_flow_direction)
    
    # Compute inverse of Contraction dissipation loss factor
    Non_linear_term_Cont_1 = (C0/Throat_rey[int_throats_nonzero_flow])**n + 1/(2**n) * (1 - (throat_rad[int_throats_nonzero_flow]**2)/(rad_con[int_throats_nonzero_flow]**2))**n
    C_inv[int_throats_nonzero_flow] = ((2 * np.pi**2)/density ) * (throat_rad[int_throats_nonzero_flow]**4)/(Q_throats[int_throats_nonzero_flow]) * 1/Non_linear_term_Cont_1
    
    # Compute inverse of Expansion dissipation loss factor    
    Non_linear_term_Expa_1 = (E0/Throat_rey[int_throats_nonzero_flow])**m + (1 - (throat_rad[int_throats_nonzero_flow]**2)/(rad_exp[int_throats_nonzero_flow]**2))**(2*m)
    E_inv[int_throats_nonzero_flow] = ((2 * np.pi**2)/density ) * (throat_rad[int_throats_nonzero_flow]**4)/(Q_throats[int_throats_nonzero_flow]) * 1/Non_linear_term_Expa_1
    
    # Assemble total hydraulic conductance
    Hydraulic_conductance_network = A_inv + C_inv + E_inv 
    
    return Hydraulic_conductance_network


def Throat_pore_velocity_extraction(target, Q_profile, New_HC, Odd, Stokes_flow_algorithm):
    r'''
    Function computes:
        1. Throat flow rate and absolute velocity
        2. Pore absolute velocity and directional velocity
    and assigns these to the loaded in network. This is done following the approach
    presented by F. Larachi et al.   - http://dx.doi.org/10.1016/j.cej.2013.11.077 

    Parameters
    ----------
    target : 
        Network
    Q_profile : ndarray, float
        Computed throat flow rate through network.
    New_HC : 
        Indicator for physics package:
            if New_HC == 0 -> You used OpenPNM physics
            if New_HC == 1 -> You used the physics of Larachi et al.
    Odd : 
        Indicator for physics package:
            if Odd == 0 -> Even network
            if Odd == 1 -> Odd network
            
    Stokes_flow_algorithm: OpenPNM object
        The stokesflow algorithm containing the fluid field data
        
    Returns
    -------
    None.
    '''
    
    # Obtain network and corresponding phase
    np.warnings.filterwarnings('ignore', category=RuntimeWarning)
    net = target.project.network
    phase = net.project.phases.copy()[0]
    
    # Initializing pore/throat velocity vector [x, y, z]
    Pore_vel_matrix = np.zeros(shape = [len(net.pores('all')), 3])
    Pore_vel_magnitude_vector = np.zeros(shape = [len(net.pores('all')), 1])
    Pore_vel_magnitude_vector_check = np.zeros(shape = [len(net.pores('all')), 1])
    Throat_vel_vector = np.zeros(shape = [len(net.throats('all')), 1])
    Throat_vel_vector = Q_profile/(np.pi * (net['throat.diameter']/2)**2)

    # Throat velocity OpenPNM physics
    if New_HC == 0:
        if Odd == 0:
             net['throat.absolute_velocity_OpenPNM_even_network'] = Throat_vel_vector.copy()
             net['throat.absolute_flow_rate_OpenPNM_even_network'] = Q_profile.copy()
             Stokes_flow_algorithm['throat.absolute_velocity_OpenPNM_even_network'] = Throat_vel_vector.copy()
             Stokes_flow_algorithm['throat.absolute_flow_rate_OpenPNM_even_network'] = Q_profile.copy()
        else:
             net['throat.absolute_velocity_OpenPNM_odd_network'] = Throat_vel_vector.copy()
             net['throat.absolute_flow_rate_OpenPNM_odd_network'] = Q_profile.copy()
             Stokes_flow_algorithm['throat.absolute_velocity_OpenPNM_odd_network'] = Throat_vel_vector.copy()
             Stokes_flow_algorithm['throat.absolute_flow_rate_OpenPNM_odd_network'] = Q_profile.copy()
             
    # Throat velocity Larachi physics         
    elif New_HC == 1:
        if Odd == 0:
            net['throat.absolute_velocity_even_network'] = Throat_vel_vector.copy()
            net['throat.absolute_flow_rate_even_network'] = Q_profile.copy()
            Stokes_flow_algorithm['throat.absolute_velocity_even_network'] = Throat_vel_vector.copy()
            Stokes_flow_algorithm['throat.absolute_flow_rate_even_network'] = Q_profile.copy()
        else:
            net['throat.absolute_velocity_odd_network'] = Throat_vel_vector.copy()
            net['throat.absolute_flow_rate_odd_network'] = Q_profile.copy()
            Stokes_flow_algorithm['throat.absolute_velocity_odd_network'] = Throat_vel_vector.copy()
            Stokes_flow_algorithm['throat.absolute_flow_rate_odd_network'] = Q_profile
            
    # Initializing pore directional vector components [x, y, z]
    U_x_i_alt = np.zeros(shape = [len(net.throats('all')), 1])
    U_y_i_alt = np.zeros(shape = [len(net.throats('all')), 1])
    U_z_i_alt = np.zeros(shape = [len(net.throats('all')), 1])

    # Compute Center to Center length
    throats_all = net.throats('all')
    throat_conns = net['throat.conns'][throats_all]     # Throat connections
    C1 = net['pore.coords'][throat_conns[:, 0]]         # Coords of pore 1
    C2 = net['pore.coords'][throat_conns[:, 1]]         # Coords of pore 2
    L_ctc = np.sqrt(((C1 - C2)**2).sum(axis=1))
    
    # Loop over the internal pores (Only internal pores contain non-zero reaction source term)
    for idx in net.pores('internal'):
        
        # Finding discharging throats (channels)
        Outflow_throats_i = np.where(net['throat.Pore_Contraction'] == idx)[0]
        
        # Velocity can only be computed for those pores which contain discharging throats
        if len(Outflow_throats_i) > 0: 
            # Total discharged volumetric flow from pore i
            Q_outflow_throats_i_sum = Q_profile[Outflow_throats_i].sum()
            
            # Extracting coordinates of pore i
            x_i = net['pore.coords'][idx][0]; y_i = net['pore.coords'][idx][1]; z_i = net['pore.coords'][idx][2]
                  
            # Initializing the [x, y, z]-component of the unit vector dictating the 
            # flow direction within pore i. 
            U_i_x = 0; U_i_y = 0; U_i_z = 0; delta_x = 0; delta_y = 0; delta_z = 0; U_x_i =0; U_y_i = 0; U_z_i =0;
            
            for it in Outflow_throats_i:
                Outflow_pore_j = net['throat.Pore_Expansion'][it]
                # Extracting unit vector flow direction for each cartesion component
                # per discharging flow rate. Note that due to numerical accuracy 
                # (accuracy of 10e-20 required!) we compute U_i_x (used for pore directional vector, || U_{i} ||)
                # and U_x_i (used for pore directional vector *components*) seperately. 
                
                delta_x = net['pore.coords'][Outflow_pore_j][0] - x_i
                delta_y = net['pore.coords'][Outflow_pore_j][1] - y_i
                delta_z = net['pore.coords'][Outflow_pore_j][2] - z_i
                            
                # Computing components of pore directional vector || U_{i} || (eqn 29)
                Q_outflow_throats_j = Q_profile[it]
                U_i_x = U_i_x + (Q_outflow_throats_j**2)/L_ctc[it] * delta_x 
                U_i_y = U_i_y + (Q_outflow_throats_j**2)/L_ctc[it] * delta_y
                U_i_z = U_i_z + (Q_outflow_throats_j**2)/L_ctc[it] * delta_z
                
            # Computing pore directional vector squared|| U_{i} || (eqn 29)
            U_i_pow2 = (U_i_x**2 + U_i_y**2 + U_i_z**2) / (Q_outflow_throats_i_sum**4)        
            U_i_sqrt = np.sqrt(U_i_pow2)
            
            # Computing [x, y, z]-component of the unit vector dictating the 
            # flow direction within pore i. 
            U_x_i_alt[idx] = U_i_x / (Q_outflow_throats_i_sum**2 * U_i_sqrt)  
            U_y_i_alt[idx] = U_i_y / (Q_outflow_throats_i_sum**2 * U_i_sqrt)  
            U_z_i_alt[idx] = U_i_z / (Q_outflow_throats_i_sum**2 * U_i_sqrt) 
            
            # Computing unit vector length (check if pore directional vector is computed well)
            Lenght_pore_directional_vector = np.sqrt((U_x_i_alt[idx])**2 + (U_y_i_alt[idx])**2 + (U_z_i_alt[idx])**2)
            Deviation_from_unity = 1.0 - Lenght_pore_directional_vector
            if Deviation_from_unity[0] > 1.0e-10:
                print(f'Pore {idx} directional vector is {Lenght_pore_directional_vector[0]:.5f} -> NOT EQUAL TO UNITY!')
            
            # Computing pore velocity components [x, y, z]
            Pore_vel_matrix[idx,0] = Q_outflow_throats_i_sum / (np.pi * (net['pore.diameter'][idx]/2)**2) * U_x_i_alt[idx]
            Pore_vel_matrix[idx,1] = Q_outflow_throats_i_sum / (np.pi * (net['pore.diameter'][idx]/2)**2) * U_y_i_alt[idx]
            Pore_vel_matrix[idx,2] = Q_outflow_throats_i_sum / (np.pi * (net['pore.diameter'][idx]/2)**2) * U_z_i_alt[idx]

            Pore_vel_magnitude_vector[idx] = np.sqrt( (Pore_vel_matrix[idx,0])**2 + (Pore_vel_matrix[idx,1])**2 + (Pore_vel_matrix[idx,2])**2)
            Pore_vel_magnitude_vector_check[idx] = Q_outflow_throats_i_sum / (np.pi * (net['pore.diameter'][idx]/2)**2)
   
    # Assign pore velocity for OpenPNM physics    
    if New_HC ==0:
        if Odd == 0:
            # Writing cartesian velocity components and the velocity magnitude to the network
            net['pore.velocity_components_OpenPNM_even_network'] = net['pore.centroid'].copy()
            net['pore.velocity_components_OpenPNM_even_network'][:,0] = Pore_vel_matrix[:,0].copy()                         # X-velocity component
            net['pore.velocity_components_OpenPNM_even_network'][:,1] = Pore_vel_matrix[:,1].copy()                         # Y-velocity component
            net['pore.velocity_components_OpenPNM_even_network'][:,2] = Pore_vel_matrix[:,2].copy()                         # Z-velocity component
            net['pore.velocity_magnitude_OpenPNM_even_network'] = Pore_vel_magnitude_vector.copy()                          # Velocity magnitude
            Stokes_flow_algorithm['pore.velocity_components_OpenPNM_even_network'] = net['pore.centroid'].copy()
            Stokes_flow_algorithm['pore.velocity_components_OpenPNM_even_network'][:,0] = Pore_vel_matrix[:,0].copy()       # X-velocity component
            Stokes_flow_algorithm['pore.velocity_components_OpenPNM_even_network'][:,1] = Pore_vel_matrix[:,1].copy()       # Y-velocity component
            Stokes_flow_algorithm['pore.velocity_components_OpenPNM_even_network'][:,2] = Pore_vel_matrix[:,2].copy()       # Z-velocity component
            Stokes_flow_algorithm['pore.velocity_magnitude_OpenPNM_even_network'] = Pore_vel_magnitude_vector.copy()        # Velocity magnitude
        else:
            # Writing cartesian velocity components and the velocity magnitude to the network
            net['pore.velocity_components_OpenPNM_odd_network'] = net['pore.centroid'].copy()
            net['pore.velocity_components_OpenPNM_odd_network'][:,0] = Pore_vel_matrix[:,0].copy()                          # X-velocity component
            net['pore.velocity_components_OpenPNM_odd_network'][:,1] = Pore_vel_matrix[:,1].copy()                          # Y-velocity component
            net['pore.velocity_components_OpenPNM_odd_network'][:,2] = Pore_vel_matrix[:,2].copy()                          # Z-velocity component
            net['pore.velocity_magnitude_OpenPNM_odd_network'] = Pore_vel_magnitude_vector.copy()                           # Velocity magnitude
            Stokes_flow_algorithm['pore.velocity_components_OpenPNM_odd_network'] = net['pore.centroid'].copy()
            Stokes_flow_algorithm['pore.velocity_components_OpenPNM_odd_network'][:,0] = Pore_vel_matrix[:,0].copy()       # X-velocity component
            Stokes_flow_algorithm['pore.velocity_components_OpenPNM_odd_network'][:,1] = Pore_vel_matrix[:,1].copy()       # Y-velocity component
            Stokes_flow_algorithm['pore.velocity_components_OpenPNM_odd_network'][:,2] = Pore_vel_matrix[:,2].copy()       # Z-velocity component
            Stokes_flow_algorithm['pore.velocity_magnitude_OpenPNM_odd_network'] = Pore_vel_magnitude_vector.copy()        # Velocity magnitude
    
    # Assign pore velocity for Larachi physics    
    elif New_HC == 1:
        if Odd == 0:
            # Writing cartesian velocity components and the velocity magnitude to the network
            net['pore.velocity_components_even_network'] = net['pore.centroid'].copy()
            net['pore.velocity_components_even_network'][:,0] = Pore_vel_matrix[:,0].copy()                                 # X-velocity component
            net['pore.velocity_components_even_network'][:,1] = Pore_vel_matrix[:,1].copy()                                 # Y-velocity component
            net['pore.velocity_components_even_network'][:,2] = Pore_vel_matrix[:,2].copy()                                 # Z-velocity component
            net['pore.velocity_magnitude_even_network'] = Pore_vel_magnitude_vector.copy()                                  # Velocity magnitude
            Stokes_flow_algorithm['pore.velocity_components_even_network'] = net['pore.centroid'].copy()
            Stokes_flow_algorithm['pore.velocity_components_even_network'][:,0] = Pore_vel_matrix[:,0].copy()               # X-velocity component
            Stokes_flow_algorithm['pore.velocity_components_even_network'][:,1] = Pore_vel_matrix[:,1].copy()               # Y-velocity component
            Stokes_flow_algorithm['pore.velocity_components_even_network'][:,2] = Pore_vel_matrix[:,2].copy()               # Z-velocity component
            Stokes_flow_algorithm['pore.velocity_magnitude_even_network'] = Pore_vel_magnitude_vector.copy()                # Velocity magnitude
        else:
            # Writing cartesian velocity components and the velocity magnitude to the network
            net['pore.velocity_components_odd_network'] = net['pore.centroid'].copy()
            net['pore.velocity_components_odd_network'][:,0] = Pore_vel_matrix[:,0].copy()                                  # X-velocity component
            net['pore.velocity_components_odd_network'][:,1] = Pore_vel_matrix[:,1].copy()                                  # Y-velocity component
            net['pore.velocity_components_odd_network'][:,2] = Pore_vel_matrix[:,2].copy()                                  # Z-velocity component
            net['pore.velocity_magnitude_odd_network'] = Pore_vel_magnitude_vector.copy()                                   # Velocity magnitude  
            Stokes_flow_algorithm['pore.velocity_components_odd_network'] = net['pore.centroid'].copy()
            Stokes_flow_algorithm['pore.velocity_components_odd_network'][:,0] = Pore_vel_matrix[:,0].copy()                # X-velocity component
            Stokes_flow_algorithm['pore.velocity_components_odd_network'][:,1] = Pore_vel_matrix[:,1].copy()                # Y-velocity compone
            Stokes_flow_algorithm['pore.velocity_components_odd_network'][:,2] = Pore_vel_matrix[:,2].copy()                # Z-velocity component
            Stokes_flow_algorithm['pore.velocity_magnitude_odd_network'] = Pore_vel_magnitude_vector.copy()                 # Velocity magnitude  

def Export_stokesflow_algo(SF_algorithm, export_project, Odd):
    r'''
    Function writes all properties of the SF_algorithm to the to be exported 
    export_project.

    Parameters
    ----------
    SF_algorithm : 
        The stokes flow algorithm.
    export_project : TYPE
        The to be exported file.
    Odd : TYPE
        Indicator for odd and even network
            Odd = 0 -> Even network.
            Odd = 1 -> Odd network.
    -------    '''    

    export_project['pore.bc.rate'] = SF_algorithm['pore.bc.rate']
    export_project['pore.bc.value'] = SF_algorithm['pore.bc.value']
    export_project['pore.pressure'] = SF_algorithm['pore.pressure']
    export_project['pore.initial_guess'] = SF_algorithm['pore.initial_guess']
    
    if Odd == 0:
        export_project['throat.absolute_velocity_OpenPNM_even_network'] = SF_algorithm['throat.absolute_velocity_OpenPNM_even_network']
        export_project['throat.absolute_flow_rate_OpenPNM_even_network'] = SF_algorithm['throat.absolute_flow_rate_OpenPNM_even_network']
        export_project['pore.velocity_components_OpenPNM_even_network'] = SF_algorithm['pore.velocity_components_OpenPNM_even_network']
        export_project['pore.velocity_magnitude_OpenPNM_even_network'] = SF_algorithm['pore.velocity_magnitude_OpenPNM_even_network']
    else:
        export_project['throat.absolute_velocity_OpenPNM_odd_network'] = SF_algorithm['throat.absolute_velocity_OpenPNM_odd_network']
        export_project['throat.absolute_flow_rate_OpenPNM_odd_network'] = SF_algorithm['throat.absolute_flow_rate_OpenPNM_odd_network']
        export_project['pore.velocity_components_OpenPNM_odd_network'] = SF_algorithm['pore.velocity_components_OpenPNM_odd_network']
        export_project['pore.velocity_magnitude_OpenPNM_odd_network'] = SF_algorithm['pore.velocity_magnitude_OpenPNM_odd_network']
    
    export_project['throat.Hydraulic_conductance_OpenPNM'] = SF_algorithm['throat.Hydraulic_conductance_OpenPNM']
    
    if Odd == 0:
        export_project['throat.absolute_velocity_even_network'] = SF_algorithm['throat.absolute_velocity_even_network']
        export_project['throat.absolute_flow_rate_even_network'] = SF_algorithm['throat.absolute_flow_rate_even_network']
        export_project['pore.velocity_components_even_network'] = SF_algorithm['pore.velocity_components_even_network']
        export_project['pore.velocity_magnitude_even_network'] = SF_algorithm['pore.velocity_magnitude_even_network']
    else:
        export_project['throat.absolute_velocity_odd_network'] = SF_algorithm['throat.absolute_velocity_odd_network']
        export_project['throat.absolute_flow_rate_odd_network'] = SF_algorithm['throat.absolute_flow_rate_odd_network']
        export_project['pore.velocity_components_odd_network'] = SF_algorithm['pore.velocity_components_odd_network']
        export_project['pore.velocity_magnitude_odd_network'] = SF_algorithm['pore.velocity_magnitude_odd_network']
       
    export_project['throat.Hydraulic_conductance_Larachi'] = SF_algorithm['throat.Hydraulic_conductance_Larachi']
    pores_all = SF_algorithm.pores('all')
    export_project.set_label(label = 'all', pores = pores_all)
    throats_all = SF_algorithm.throats('all')
    export_project.set_label(label = 'all', throats = throats_all) 


def Import_stokesflow_algo(SF_algorithm, Odd, path_sf):
    r'''
    Function writes all properties of the imported_project to the to newly 
    created stokes flow algorithm
    
    Parameters
    ----------
    SF_algorithm : 
        The stokes flow algorithm on which you want to write the data
    Odd : TYPE
        Indicator for odd and even network
            Odd = 0 -> Even network.
            Odd = 1 -> Odd network.
    path_sf:
        Path of to be imported stokesflow
    -------    '''    
    
    sf_project_import_proj = op.Workspace().load_project(filename = path_sf)
    Import_SF = sf_project_import_proj['net']
    
    SF_algorithm['pore.bc.rate'] = Import_SF['pore.bc.rate']
    SF_algorithm['pore.bc.value'] = Import_SF['pore.bc.value']
    SF_algorithm['pore.pressure'] = Import_SF['pore.pressure']
    SF_algorithm['pore.initial_guess'] = Import_SF['pore.initial_guess']
    
    if Odd == 0:
        SF_algorithm['throat.absolute_velocity_OpenPNM_even_network'] = Import_SF['throat.absolute_velocity_OpenPNM_even_network']
        SF_algorithm['throat.absolute_flow_rate_OpenPNM_even_network'] = Import_SF['throat.absolute_flow_rate_OpenPNM_even_network']
        SF_algorithm['pore.velocity_components_OpenPNM_even_network'] = Import_SF['pore.velocity_components_OpenPNM_even_network']
        SF_algorithm['pore.velocity_magnitude_OpenPNM_even_network'] = Import_SF['pore.velocity_magnitude_OpenPNM_even_network']
    else:
        SF_algorithm['throat.absolute_velocity_OpenPNM_odd_network'] = Import_SF['throat.absolute_velocity_OpenPNM_odd_network']
        SF_algorithm['throat.absolute_flow_rate_OpenPNM_odd_network'] = Import_SF['throat.absolute_flow_rate_OpenPNM_odd_network']
        SF_algorithm['pore.velocity_components_OpenPNM_odd_network'] = Import_SF['pore.velocity_components_OpenPNM_odd_network']
        SF_algorithm['pore.velocity_magnitude_OpenPNM_odd_network'] = Import_SF['pore.velocity_magnitude_OpenPNM_odd_network']
    
    SF_algorithm['throat.Hydraulic_conductance_OpenPNM'] = Import_SF['throat.Hydraulic_conductance_OpenPNM']
    
    if Odd == 0:
        SF_algorithm['throat.absolute_velocity_even_network'] = Import_SF['throat.absolute_velocity_even_network']
        SF_algorithm['throat.absolute_flow_rate_even_network'] = Import_SF['throat.absolute_flow_rate_even_network']
        SF_algorithm['pore.velocity_components_even_network'] = Import_SF['pore.velocity_components_even_network']
        SF_algorithm['pore.velocity_magnitude_even_network'] = Import_SF['pore.velocity_magnitude_even_network']
    else:
        SF_algorithm['throat.absolute_velocity_odd_network'] = Import_SF['throat.absolute_velocity_odd_network']
        SF_algorithm['throat.absolute_flow_rate_odd_network'] = Import_SF['throat.absolute_flow_rate_odd_network']
        SF_algorithm['pore.velocity_components_odd_network'] = Import_SF['pore.velocity_components_odd_network']
        SF_algorithm['pore.velocity_magnitude_odd_network'] = Import_SF['pore.velocity_magnitude_odd_network']
       
    SF_algorithm['throat.Hydraulic_conductance_Larachi'] = Import_SF['throat.Hydraulic_conductance_Larachi']
    pores_all = Import_SF.pores('all')
    SF_algorithm.set_label(label = 'all', pores = pores_all)
    throats_all = Import_SF.throats('all')
    SF_algorithm.set_label(label = 'all', throats = throats_all) 


def Update_phase_with_SF(SF_algorithm, phase, Odd, Larachi):
    r'''
    Function writes all properties of the SF_algorithm to the to phase.

    Parameters
    ----------
    SF_algorithm : 
        The stokes flow algorithm.
    phase : 
        The phase which you want the stokes flow algorithm data to be written to
    Odd : TYPE
        Indicator for odd and even network
            Odd = 0 -> Even network.
            Odd = 1 -> Odd network.
    Larachi : 
        Indicator for which physics you want to use for the hydraulic conductance
            Larachi = 0 -> OpenPNM physics
            Larachi = 1 -> Larachi physics
    Returns
    -------
    None.
    '''
    
    phase['pore.pressure'] = SF_algorithm['pore.pressure']
    if Larachi == 0:
        phase['throat.hydraulic_conductance'] = SF_algorithm['throat.Hydraulic_conductance_OpenPNM']
    elif Larachi == 1:
        phase['throat.hydraulic_conductance'] = SF_algorithm['throat.Hydraulic_conductance_Larachi']      
        
    if Odd == 0:
        phase['throat.absolute_velocity_OpenPNM_even_network'] = SF_algorithm['throat.absolute_velocity_OpenPNM_even_network']
        phase['throat.absolute_flow_rate_OpenPNM_even_network'] = SF_algorithm['throat.absolute_flow_rate_OpenPNM_even_network']
        phase['pore.velocity_components_OpenPNM_even_network'] = SF_algorithm['pore.velocity_components_OpenPNM_even_network']
        phase['pore.velocity_magnitude_OpenPNM_even_network'] = SF_algorithm['pore.velocity_magnitude_OpenPNM_even_network']
    else:
        phase['throat.absolute_velocity_OpenPNM_odd_network'] = SF_algorithm['throat.absolute_velocity_OpenPNM_odd_network']
        phase['throat.absolute_flow_rate_OpenPNM_odd_network'] = SF_algorithm['throat.absolute_flow_rate_OpenPNM_odd_network']
        phase['pore.velocity_components_OpenPNM_odd_network'] = SF_algorithm['pore.velocity_components_OpenPNM_odd_network']
        phase['pore.velocity_magnitude_OpenPNM_odd_network'] = SF_algorithm['pore.velocity_magnitude_OpenPNM_odd_network']
    
    if Odd == 0:
        phase['throat.absolute_velocity_even_network'] = SF_algorithm['throat.absolute_velocity_even_network']
        phase['throat.absolute_flow_rate_even_network'] = SF_algorithm['throat.absolute_flow_rate_even_network']
        phase['pore.velocity_components_even_network'] = SF_algorithm['pore.velocity_components_even_network']
        phase['pore.velocity_magnitude_even_network'] = SF_algorithm['pore.velocity_magnitude_even_network']
    else:
        phase['throat.absolute_velocity_odd_network'] = SF_algorithm['throat.absolute_velocity_odd_network']
        phase['throat.absolute_flow_rate_odd_network'] = SF_algorithm['throat.absolute_flow_rate_odd_network']
        phase['pore.velocity_components_odd_network'] = SF_algorithm['pore.velocity_components_odd_network']
        phase['pore.velocity_magnitude_odd_network'] = SF_algorithm['pore.velocity_magnitude_odd_network']
    
    
def Mass_transfer_coefficient(network, phase, Odd, Vel_dependent, Larachi, param):    
    r'''
    Function computes pore local mass transfer coefficient for:
        1. Velocity-independent correlation : km_{local} = R_{pore, i} / D_{active species}
        2. Velocity-dependent correlation : km_{local} = C1 * velocity_{pore}^{C2}.
        
    Parameters
    ----------
    network : OpenPNM network

    phase : OpenPNM phase object associated with network

    Odd : int
        Indicate odd or even network.
            Odd = 0 -> Even network
            Odd = 1 -> Odd network
            
    Vel_dependent: int
        Indicator for use of velocity dependent or independent local mass transfer coefficient
        Vel_dependent == 0: Velocity independent mass transfer 
            correlation km_{local} = R_{pore, i} / D_{active species}. Correlation 
            is adapted from Van der Heijden et al. DOI: 10.1149/1945-7111/ac5e46 
            
        Vel_dependent == 1: Velocity dependent mass transfer 
            correlation km_{local} = C1 * velocity_{pore}^{C2}. Here C1 and C2
            are fitted, and velocity_{pore} the velocity in the *pores*, which 
            is derived following the methodology presented by Larachi et al. 
            (http://dx.doi.org/10.1016/j.cej.2013.11.077 ).
            
    Larachi_HC: int
        Indicate odd or even network.
            Larachi == 1: OpenPNM physics are used (i.e. shape factor)
            Larachi == 1: Physics based on mechanical energy balance are used (Larachi et al.)  
            
    parameter_script : dict
        Dictionairy with phase properties
    '''
    
    # Velocity independent:
    if Vel_dependent == 0:
        rad_pore = network['pore.diameter'] / 2             # Pore radii in the cathode [m]
        phase['pore.km_exp'] = param['D_c'] / rad_pore
    
    # Velocity dependent:
    if Vel_dependent == 1:
        C1 = param['MT_coeff_C1']
        C2 = param['MT_coeff_C2']
        
        if Larachi == 0:
            if Odd == 0: 
                pore_vel = phase['pore.velocity_magnitude_OpenPNM_even_network'][:,0]
            elif Odd == 1:
                pore_vel = phase['pore.velocity_magnitude_OpenPNM_odd_network'][:,0]
        else:
            if Odd == 0: 
                pore_vel = phase['pore.velocity_magnitude_even_network'][:,0]
            elif Odd == 1:
                pore_vel = phase['pore.velocity_magnitude_odd_network'][:,0]
        
        phase['pore.km_exp'] = C1 * pore_vel**(C2)
        
        # Account for internal pores with zero pore velocity. We assign these pores the 
        # Velocity in-dependent km (according to the thin film theory)
        Zero_vel_pores = network.pores('internal')[np.where(pore_vel[network.pores('internal')] == 0)]
        rad_pore = network['pore.diameter'] / 2
        phase['pore.km_exp'][Zero_vel_pores] = param['D_c'] / rad_pore[Zero_vel_pores]
        
        # Check if there are still pores with a zero MT-coefficient
        No_internal_pores_zero_MT = len(np.where(phase['pore.km_exp'][network.pores('internal')] == 0)[0])
        
        if No_internal_pores_zero_MT == 0:
                pass
        else:
            print("Some internal pores are assigned a zero mass transfer coefficient!")
            sys.exit()
            
            
def Local_mass_transfer_coefficient_fitting(network, phase, Odd, alpha, beta):  
    r'''
    Function for fitting LOCAL mass transfer coefficient:
        km_local = alpha * (pore-velocity) ^ {Beta}
        
    * Note that currently the km of a pore without velocity is computed as zero.
        
    Parameters
    ----------
    network : OpenPNM network

    phase : OpenPNM phase object associated with network

    Odd : 
        Indicate odd or even network.
            Odd = 0 -> Even network
            Odd = 1 -> Odd network
    alpha : 
        Pre-exponential fitting constant.
    beta : 
        Exponential fitting constant.
    '''
    
    # Make sure that the values are properly overwritten
    phase['pore.km_loc'] = np.nan
    
    if Odd == 0: 
        pore_vel = phase['pore.velocity_magnitude_even_network'][:,0]
    elif Odd == 1:
        pore_vel = phase['pore.velocity_magnitude_odd_network'][:,0]

    phase['pore.km_loc'] = alpha * pore_vel**(beta)    
    
    
'''---------- OpenPNM physics mass transfer computation functions ----------'''    
# Note that these functions are copied from the functions above. However, certain
# parts are removed as they are not required for OpenPNM physics. 
 
def Throat_flowrate_OpenPNM(target, 
                            throat_area='throat.area',
                            pore_diameter='pore.diameter',
                            throat_diameter='throat.diameter',
                            Hydraulic_conductance='throat.hydraulic_conductance'):
    r'''
    Function 
        1. Assigns the contraction and expansion pore per conduit throat and and the direction of the flow
            If there is no flow between the two conduit pores, a value of *-1* is assigned.
        2. Computes absolute pressure difference between the two connected pores per conduit
        3. Computes absolute throat flowrate beween the two connected pores per conduit

    Parameters
    ----------
    target : OpenPNM network
    
    throat_area : string
        Dictionary key of throat cross-sectional area.
        
    pore_diameter : string
        Dictionary key of the pore diameter values.

    throat_diameter : string
        Dictionary key of the throat diameter values.
        
    Hydraulic_conductance : ndarray, float
        The throat hydraulic conductancte 

    Returns
    -------
    Q_throats : ndarray, float
        Absolute flow rate in throat of each conduit [m3/s].
        
    Abs_pressure_difference : ndarray, float
        Absolute pressure drop in throat of each conduit [Pa].
    '''
    
    # Obtain network and corresponding phase
    np.warnings.filterwarnings('ignore', category=RuntimeWarning)
    network = target.project.network
    phase = network.project.phases.copy()[0]
    
    # Conduit throat hydraulic conductance
    throats = network.throats('all')
    cn = network['throat.conns'][throats]
    Hyd_con_throats = phase[Hydraulic_conductance][throats]
    
    # Pressure in each pore of conduit and pressure drop over conduit
    Pressure_pore_1 = phase['pore.pressure'][cn[:, 0]]
    Pressure_pore_2 = phase['pore.pressure'][cn[:, 1]]
    Pres_diff = Pressure_pore_1 - Pressure_pore_2
    
    # Absolute pressure drop and flow rate in conduit
    Abs_pressure_difference = np.abs(Pres_diff)
    Q_throats = Abs_pressure_difference * Hyd_con_throats
    abs_Q_throats = np.abs(Q_throats)
    
    ''' Assigning per conduit throat the contraction and expansion pore, 
        and the direction of the flow. The following methodology is maintained:
            
        Per conduit number the pressure of its connecting pores (termed here 
        Pore 1 and 2) are checked. Their difference is set as: 
            Pres_diff = Pressure_pore_1 - Pressure_pore_2
        
        If conduit Pressure_pore_1 < Pressure_pore_2 -> flow from pore 2 to pore 1
            We assign network['throat.Flow_direction'] = 12 (flow from pore 1 -> pore 2).
            network['throat.Pore_Contraction']  = Pore number 2
            network['throat.Pore_Expansion'][throats_idx] = Pore number 1
            
         If conduit Pressure_pore_2 < Pressure_pore_1 -> flow from pore 1 to pore 2
             We assign network['throat.Flow_direction'] = 21 (flow from pore 2 -> pore 1).
             network['throat.Pore_Contraction']  = Pore number 1
             network['throat.Pore_Expansion']= Pore_2
             
        If there is no pressure difference between them, there will also be no
        flow between pore 1 and 2. Then the network['throat.Pore_Contraction'], 
        ['throat.Pore_Contraction'] and ['throat.Pore_Expansion'] are assigned 
        a value of *-1*. 
        '''
    
    # Initializing the contraction, expansion and flow direction arrays for each pore-throat-pore element
    network['throat.Pore_Contraction'] = np.nan
    network['throat.Pore_Expansion'] = np.nan
    network['throat.Flow_direction'] = np.nan
    
    # Asign contraction (inflow) and expansion (outflow) pores for each pore-throat-pore element
    # Negative pressure drop
    # If conduit Pressure_pore_1 < Pressure_pore_2 -> flow from pore 2 to pore 1
    idx_neg_pres = np.where(Pres_diff<0)[0]
    Pore_1_idx_neg_pres = cn[idx_neg_pres, 0]
    Pore_2_idx_neg_pres = cn[idx_neg_pres, 1]
    network['throat.Pore_Contraction'][idx_neg_pres] = Pore_2_idx_neg_pres
    network['throat.Pore_Expansion'][idx_neg_pres] = Pore_1_idx_neg_pres
    network['throat.Flow_direction'][idx_neg_pres] = 21        

    # Positive pressure drop
    # If conduit Pressure_pore_1 > Pressure_pore_2 -> flow from pore 1 to pore 2
    idx_pos_pres = np.where(Pres_diff>0)[0]
    Pore_1_idx_pos_pres = cn[idx_pos_pres, 0]
    Pore_2_idx_pos_pres = cn[idx_pos_pres, 1]
    network['throat.Pore_Contraction'][idx_pos_pres] = Pore_1_idx_pos_pres
    network['throat.Pore_Expansion'][idx_pos_pres] = Pore_2_idx_pos_pres
    network['throat.Flow_direction'][idx_pos_pres] = 12 
    
    # No pressure drop
    # If conduit Pressure_pore_1 == Pressure_pore_2 -> No flow
    idx_no_pres = np.where(Pres_diff == 0)[0]
    network['throat.Pore_Contraction'][idx_no_pres] = -1
    network['throat.Pore_Expansion'][idx_no_pres] = -1
    network['throat.Flow_direction'][idx_no_pres] = -1

    # Converting labels from float to integers.
    network['throat.Pore_Contraction'] = network['throat.Pore_Contraction'].astype(int)
    network['throat.Pore_Expansion'] = network['throat.Pore_Expansion'].astype(int)
    network['throat.Flow_direction'] = network['throat.Flow_direction'].astype(int)
    
    return abs_Q_throats, Abs_pressure_difference


def Throat_pore_velocity_extraction_OpenPNM(target, Odd, Stokes_flow_algorithm):
    r'''
    Function computes:
        1. Throat flow rate and absolute velocity
        2. Pore absolute velocity 
    and assigns these to the loaded in network. The pore velocity is computed
    following the approach presented by F. Larachi et al.   
    - http://dx.doi.org/10.1016/j.cej.2013.11.077 

    Parameters
    ----------
    target : 
        Network
    Odd : 
        Indicator for network:
            if Odd == 0 -> Even network
            if Odd == 1 -> Odd network
            
    Stokes_flow_algorithm: OpenPNM object
        The stokesflow algorithm containing the fluid field data
        
    Returns
    -------
    None.
    '''
    
    # Obtain network and corresponding phase
    np.warnings.filterwarnings('ignore', category=RuntimeWarning)
    net = target.project.network
    phase = net.project.phases.copy()[0]
    
    # Compute the throat flowrate and extract the direction of the flow 
    Q_profile, Abs_pres_diff = Throat_flowrate_OpenPNM(net, 
                                                       throat_area='throat.area',
                                                       pore_diameter='pore.diameter',
                                                       throat_diameter='throat.diameter',
                                                       Hydraulic_conductance='throat.hydraulic_conductance')

    # Initializing pore/throat velocity vector [x, y, z]
    Pore_vel_magnitude_vector = np.zeros(shape = [len(net.pores('all')), 1])
    Throat_vel_vector = np.zeros(shape = [len(net.throats('all')), 1])
    Throat_vel_vector = Q_profile/(np.pi * (net['throat.diameter']/2)**2)

    # Throat velocity:
    if Odd == 0:
         net['throat.absolute_velocity_even_network'] = Throat_vel_vector.copy()
         net['throat.absolute_flow_rate_even_network'] = Q_profile.copy()
         Stokes_flow_algorithm['throat.absolute_velocity_even_network'] = Throat_vel_vector.copy()
         Stokes_flow_algorithm['throat.absolute_flow_rate_even_network'] = Q_profile.copy()
    else:
         net['throat.absolute_velocity_odd_network'] = Throat_vel_vector.copy()
         net['throat.absolute_flow_rate_odd_network'] = Q_profile.copy()
         Stokes_flow_algorithm['throat.absolute_velocity_odd_network'] = Throat_vel_vector.copy()
         Stokes_flow_algorithm['throat.absolute_flow_rate_odd_network'] = Q_profile.copy()
               
    # Compute and assign the magnitude of the pore velocity. The pore velocity 
    # is computed based on sum of flow leaving each pore divided by the cross-sectional
    # area of each pore:
    # vel_{pore, i} = sum(Q_{pore i, outflow}) / Cross-sectional_area_{pore i}
    
    # Loop over the internal pores:
    for idx in net.pores('internal'):
        
        # Identifying discharging throats (channels) that are connected to a 
        # contraction pore i. Fluid flows from a contraction pore i to an 
        # expansion pore j.  
        Outflow_throats_i = np.where(net['throat.Pore_Contraction'] == idx)[0]
        
        # Velocity can only be computed for those pores which contain discharging throats
        if len(Outflow_throats_i) > 0: 
            
            # Total discharged volumetric flow from pore i
            Q_outflow_throats_i_sum = Q_profile[Outflow_throats_i].sum()
            
            Pore_vel_magnitude_vector[idx] = Q_outflow_throats_i_sum / (np.pi * (net['pore.diameter'][idx]/2)**2)
   
    # Assign pore velocity for OpenPNM physics    
    if Odd == 0:
        # Writing cartesian velocity components and the velocity magnitude to the network
        net['pore.velocity_magnitude_even_network'] = Pore_vel_magnitude_vector.copy()                          # Velocity magnitude
        Stokes_flow_algorithm['pore.velocity_magnitude_even_network'] = Pore_vel_magnitude_vector.copy()        # Velocity magnitude
    else:
        # Writing cartesian velocity components and the velocity magnitude to the network
        net['pore.velocity_magnitude_odd_network'] = Pore_vel_magnitude_vector.copy()                           # Velocity magnitude
        Stokes_flow_algorithm['pore.velocity_magnitude_odd_network'] = Pore_vel_magnitude_vector.copy()         # Velocity magnitude
    
    
def Update_phase_with_SF_OpenPNM(SF_algorithm, phase, Odd):
    r'''
    Function writes all properties of the stokes flow algorithm (SF_algorithm)
    to the to inputted phase.

    Parameters
    ----------
    SF_algorithm : 
        The stokes flow algorithm.
    phase : 
        The phase which you want the stokes flow algorithm data to be written to
    Odd : TYPE
        Indicator for odd and even network
            Odd = 0 -> Even network.
            Odd = 1 -> Odd network.
    Returns
    -------
    None.
    '''
    
    phase['pore.pressure'] = SF_algorithm['pore.pressure']
        
    if Odd == 0:
        phase['throat.absolute_velocity_even_network'] = SF_algorithm['throat.absolute_velocity_even_network']
        phase['throat.absolute_flow_rate_even_network'] = SF_algorithm['throat.absolute_flow_rate_even_network']
        phase['pore.velocity_magnitude_even_network'] = SF_algorithm['pore.velocity_magnitude_even_network']
    else:
        phase['throat.absolute_velocity_odd_network'] = SF_algorithm['throat.absolute_velocity_odd_network']
        phase['throat.absolute_flow_rate_odd_network'] = SF_algorithm['throat.absolute_flow_rate_odd_network']
        phase['pore.velocity_magnitude_odd_network'] = SF_algorithm['pore.velocity_magnitude_odd_network']
    
    
def Mass_transfer_coefficient_OpenPNM(network, phase, Odd, Vel_dependent, param):    
    r'''
    Function computes pore local mass transfer coefficient for:
        1. Velocity-independent correlation : km_{local} = R_{pore, i} / D_{active species}
        2. Velocity-dependent correlation : km_{local} = C1 * velocity_{pore}^{C2}.
        
    Parameters
    ----------
    network : OpenPNM network

    phase : OpenPNM phase object associated with network

    Odd : int
        Indicate odd or even network.
            Odd = 0 -> Even network
            Odd = 1 -> Odd network
            
    Vel_dependent: int
        Indicator for use of velocity dependent or independent local mass transfer coefficient
        Vel_dependent == 0: Velocity independent mass transfer 
            correlation km_{local} = R_{pore, i} / D_{active species}. Correlation 
            is adapted from Van der Heijden et al. DOI: 10.1149/1945-7111/ac5e46 
            
        Vel_dependent == 1: Velocity dependent mass transfer 
            correlation km_{local} = C1 * velocity_{pore}^{C2}. Here C1 and C2
            are fitted, and velocity_{pore} the velocity in the *pores*, which 
            is derived following the methodology presented by Larachi et al. 
            (http://dx.doi.org/10.1016/j.cej.2013.11.077 ).
            
    parameter_script : dict
        Dictionairy with phase properties
    '''
    
    # Velocity independent:
    if Vel_dependent == 0:
        rad_pore = network['pore.diameter'] / 2             # Pore radii in the cathode [m]
        phase['pore.km_exp'] = param['D_c'] / rad_pore
    
    # Velocity dependent:
    if Vel_dependent == 1:
        C1 = param['MT_coeff_C1']
        C2 = param['MT_coeff_C2']
        
        if Odd == 0: 
            pore_vel = phase['pore.velocity_magnitude_even_network'][:,0]
        elif Odd == 1:
            pore_vel = phase['pore.velocity_magnitude_odd_network'][:,0]
    
        phase['pore.km_exp'] = C1 * pore_vel**(C2)
        
        # Account for internal pores with zero pore velocity. We assign these pores the 
        # Velocity in-dependent km (according to the thin film theory)
        Zero_vel_pores = network.pores('internal')[np.where(pore_vel[network.pores('internal')] == 0)]
        rad_pore = network['pore.diameter'] / 2
        phase['pore.km_exp'][Zero_vel_pores] = param['D_c'] / rad_pore[Zero_vel_pores]
        
        # Check if there are still pores with a zero MT-coefficient
        No_internal_pores_zero_MT = len(np.where(phase['pore.km_exp'][network.pores('internal')] == 0)[0])
        
        if No_internal_pores_zero_MT == 0:
                pass
        else:
            print("Some internal pores are assigned a zero mass transfer coefficient!")
            sys.exit()


def Euler_number(target, Gamma):   
    r'''
    Function computes Euler number (ratio of pressure to kinetic contribution 
    in fluid total head) following the approach presented by F. Larachi et al. 
    http://dx.doi.org/10.1016/j.cej.2013.11.077 

    Parameters
    ----------
    target : OpenPNM network
        
    Gamma : float
        flow pattern constant.

    Returns
    -------
    Eu_inv : ndarray, float
        Inverse of the throat Euler number
    '''
    
    # Obtain network and corresponding phase
    np.warnings.filterwarnings('ignore', category=RuntimeWarning)
    network = target.project.network
    phase = network.project.phases.copy()[0]   
    diam_con = network['pore.diameter'][network['throat.Pore_Contraction']]
    diam_exp = network['pore.diameter'][network['throat.Pore_Expansion']]     
    density = phase['pore.density'][0]      # [kg m-3]

    # Extract the throat flowrate and the index of throats with a non-zero flowrate
    Q_throats, Abs_pressure_difference = Throat_flowrate_total_hydraulic_conductance(network)
    Internal_throats = network.throats('internal')
    Non_zero_flow_direction = np.where(network['throat.Flow_direction'] > 0)[0]
    int_throats_nonzero_flow = np.intersect1d(Internal_throats, Non_zero_flow_direction)
    throat_all = network.throats('all')
    arr3 = np.where(diam_con[throat_all] != diam_exp[throat_all])[0]
    int_throats_nonzero_flow_unequal_radius = np.intersect1d(int_throats_nonzero_flow, arr3)
    
    # Initialization of arrays to save data
    Denom = np.zeros(shape = [len(network.throats('all'))])
    Eu = np.zeros(shape = [len(network.throats('all'))])
    Eu_inv = np.zeros(shape = [len(network.throats('all'))])
    
    # Compute inverse Euler number
    Denom[int_throats_nonzero_flow_unequal_radius] = Gamma/(2 * np.pi**2) * density * (Q_throats[int_throats_nonzero_flow_unequal_radius])**2  * ( (1/(diam_con[int_throats_nonzero_flow_unequal_radius]))**4 - (1 / (diam_exp[int_throats_nonzero_flow_unequal_radius]))**4)
    Eu[int_throats_nonzero_flow_unequal_radius] = np.abs( Abs_pressure_difference[int_throats_nonzero_flow_unequal_radius] / Denom[int_throats_nonzero_flow_unequal_radius])
    Eu_inv[int_throats_nonzero_flow_unequal_radius] = 1/Eu[int_throats_nonzero_flow_unequal_radius]    
    
    return Eu_inv