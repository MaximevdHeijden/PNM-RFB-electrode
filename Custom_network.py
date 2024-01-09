"""
Creating a new network to correctly re-assign the properties and the labels 
of the network.
"""

import openpnm as op
import numpy as np
from scipy import sqrt as _sqrt
import scipy as _sp


def Re_assign_network(network):
    r'''
    Function creates a new network with correctly assigned properties and 
    labels from the imported .CSV network. 

    Parameters
    ----------
    network : OpenPNM.network object
        The imported .CSV network.

    Returns
    -------
    network_new : OpenPNM.network object
        New network object with correctly assigned properties and labels.
    '''
    
    # Create the network
    network_new = op.network.Network()
    
    # Assign the properties
    network_new['pore.area'] = network['pore.area']
    network_new['pore.coords'] = network['pore.coords']
    network_new['throat.conns'] = network['throat.conns']
    network_new['pore.centroid'] = network['pore.centroid']
    network_new['pore.diameter'] = network['pore.diameter']
    network_new['pore.equivalent_diameter'] = network['pore.equivalent_diameter']
    network_new['pore.extended_diameter'] = network['pore.extended_diameter']
    network_new['pore.inscribed_diameter'] = network['pore.inscribed_diameter']
    network_new['pore.surface_area'] = network['pore.surface_area']
    network_new['pore.volume'] = network['pore.volume']
    network_new['throat.area'] = network['throat.area']
    network_new['throat.centroid'] = network['throat.centroid']
    network_new['throat.conduit_lengths.pore1'] = network['throat.conduit_lengths.pore1']
    network_new['throat.conduit_lengths.pore2'] = network['throat.conduit_lengths.pore2']
    network_new['throat.conduit_lengths.throat'] = network['throat.conduit_lengths.throat']
    network_new['throat.diameter'] = network['throat.diameter']
    network_new['throat.direct_length'] = network['throat.direct_length']
    network_new['throat.endpoints.head'] = network['throat.endpoints.head']
    network_new['throat.endpoints.tail'] = network['throat.endpoints.tail']
    network_new['throat.equivalent_diameter'] = network['throat.equivalent_diameter']
    network_new['throat.inscribed_diameter'] = network['throat.inscribed_diameter']
    network_new['throat.length'] = network['throat.length']
    network_new['throat.perimeter'] = network['throat.perimeter']
    network_new['throat.total_length'] = network['throat.total_length']
    network_new['throat.volume'] = network['throat.volume']

    # Assigning the labels
    All = network.pores('all')
    network_new.set_label(label = 'pore.all', pores = All)
    Back = network.pores('back')
    network_new.set_label(label = 'pore.back', pores = Back)
    Bottom = network.pores('bottom')
    network_new.set_label(label = 'pore.bottom', pores = Bottom)
    Boundary = network.pores('pore.boundary')
    network_new.set_label(label = 'boundary', pores = Boundary)
    internal = network.pores('pore.internal')
    network_new.set_label(label = 'internal', pores = internal)
    front = network.pores('front')
    network_new.set_label(label = 'pore.front', pores = front)
    label = network.pores('label')
    network_new.set_label(label = 'pore.label', pores = label)
    left = network.pores('left')
    network_new.set_label(label = 'pore.left', pores = left)
    right = network.pores('right')
    network_new.set_label(label = 'pore.right', pores = right)
    top = network.pores('top')
    network_new.set_label(label = 'pore.top', pores = top)
    tall = network.throats('all')
    network_new.set_label(label = 'throat.all', throats = tall)
    tboundary = network.throats('boundary')
    network_new.set_label(label = 'throat.boundary', throats = tboundary)   
    tinternal = network.throats('internal')
    network_new.set_label(label = 'throat.internal', throats = tinternal)   
    
    return network_new


def Re_assign_network_V2(network):
    r'''
    Function creates a new network with correctly assigned properties and 
    labels from the imported .CSV network. 

    Parameters
    ----------
    network : OpenPNM.network object
        The imported .CSV network.

    Returns
    -------
    network_new : OpenPNM.network object
        New network object with correctly assigned properties and labels.
    '''
    
    # Create the network
    network_new = op.network.Network(conns = network['throat.conns'], 
                                     coords = network['pore.coords'] )
    
    # Assign the properties
    network_new['pore.area'] = network['pore.area']
    network_new['pore.centroid'] = network['pore.centroid']
    network_new['pore.diameter'] = network['pore.diameter']
    network_new['pore.equivalent_diameter'] = network['pore.equivalent_diameter']
    network_new['pore.extended_diameter'] = network['pore.extended_diameter']
    network_new['pore.inscribed_diameter'] = network['pore.inscribed_diameter']
    network_new['pore.surface_area'] = network['pore.surface_area']
    network_new['pore.volume'] = network['pore.volume']
    network_new['throat.area'] = network['throat.area']
    network_new['throat.centroid'] = network['throat.centroid']
    network_new['throat.conduit_lengths.pore1'] = network['throat.conduit_lengths.pore1']
    network_new['throat.conduit_lengths.pore2'] = network['throat.conduit_lengths.pore2']
    network_new['throat.conduit_lengths.throat'] = network['throat.conduit_lengths.throat']
    network_new['throat.diameter'] = network['throat.diameter']
    network_new['throat.direct_length'] = network['throat.direct_length']
    network_new['throat.endpoints.head'] = network['throat.endpoints.head']
    network_new['throat.endpoints.tail'] = network['throat.endpoints.tail']
    network_new['throat.equivalent_diameter'] = network['throat.equivalent_diameter']
    network_new['throat.inscribed_diameter'] = network['throat.inscribed_diameter']
    network_new['throat.length'] = network['throat.length']
    network_new['throat.perimeter'] = network['throat.perimeter']
    network_new['throat.total_length'] = network['throat.total_length']
    network_new['throat.volume'] = network['throat.volume']

    # Assigning the labels
    All = network.pores('all')
    network_new.set_label(label = 'pore.all', pores = All)
    Back = network.pores('back')
    network_new.set_label(label = 'pore.back', pores = Back)
    Bottom = network.pores('bottom')
    network_new.set_label(label = 'pore.bottom', pores = Bottom)
    Boundary = network.pores('pore.boundary')
    network_new.set_label(label = 'boundary', pores = Boundary)
    internal = network.pores('pore.internal')
    network_new.set_label(label = 'internal', pores = internal)
    front = network.pores('front')
    network_new.set_label(label = 'pore.front', pores = front)
    label = network.pores('label')
    network_new.set_label(label = 'pore.label', pores = label)
    left = network.pores('left')
    network_new.set_label(label = 'pore.left', pores = left)
    right = network.pores('right')
    network_new.set_label(label = 'pore.right', pores = right)
    top = network.pores('top')
    network_new.set_label(label = 'pore.top', pores = top)
    tall = network.throats('all')
    network_new.set_label(label = 'throat.all', throats = tall)
    tboundary = network.throats('boundary')
    network_new.set_label(label = 'throat.boundary', throats = tboundary)   
    tinternal = network.throats('internal')
    network_new.set_label(label = 'throat.internal', throats = tinternal)   
    
    return network_new


def Center_to_center_length(target):
    r"""
    Calculate throat length assuming point-like pores, i.e. center-to-center
    distance between pores. Also, this model assumes that pores and throat
    centroids are colinear.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    Returns
    -------
    value : NumPy ndarray
        Array containing throat length values.
    """
    
    network = target.project.network
    # throats = network.map_throats(throats=target.Ts, origin=target)
    throats = network.throats('all')
    cn = network['throat.conns'][throats]
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    value = _sqrt(((C1 - C2)**2).sum(axis=1))
    
    return value


def cylinder_throat_cross_sectional_area(network, throat_diameter='throat.diameter'):
    r"""
    Calculate throat cross-sectional area for a cylindrical throat
    """
    
    diams = network[throat_diameter]
    value = np.pi / 4 * diams**2
    
    return value


def spherical_pores_conduit_lengths(target, throat_endpoints='throat.endpoints',
                    throat_length='throat.length'):
    r"""
    Calculate conduit lengths. A conduit is defined as half pore + throat
    + half pore.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_endpoints : string
        Dictionary key of the throat endpoint values.

    throat_diameter : string
        Dictionary key of the throat length values.

    throat_length : string (optional)
        Dictionary key of the throat length values.  If not given then the
        direct distance bewteen the two throat end points is used.

    Returns
    -------
    Dictionary containing conduit lengths, which can be accessed via the dict
    keys 'pore1', 'pore2', and 'throat'.
    """
    
    network = target.project.network
    # throats = network.map_throats(throats=target.Ts, origin=target)
    throats = network.throats('all')
    cn = network['throat.conns'][throats]
    # Get pore coordinates
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    # Get throat endpoints and length
    EP1 = network[throat_endpoints + '.head'][throats]
    EP2 = network[throat_endpoints + '.tail'][throats]
    try:
        # Look up throat length if given
        Lt = network[throat_length][throats]
    except KeyError:
        # Calculate throat length otherwise
        Lt = _sqrt(((EP1 - EP2)**2).sum(axis=1))
    # Calculate conduit lengths
    L1 = _sqrt(((C1 - EP1)**2).sum(axis=1))
    L2 = _sqrt(((C2 - EP2)**2).sum(axis=1))
    return {'pore1': L1, 'throat': Lt, 'pore2': L2}


def end_points_spherical_pores(target, pore_diameter='pore.diameter',
                    throat_diameter='throat.diameter',
                    throat_centroid='throat.centroid'):
    r"""
    Calculate the coordinates of throat endpoints, assuming spherical pores.
    This model accounts for the overlapping lens between pores and throats.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : string
        Dictionary key of the pore diameter values.

    throat_diameter : string
        Dictionary key of the throat diameter values.

    throat_centroid : string, optional
        Dictionary key of the throat centroid values. See the notes.

    Returns
    -------
    EP : dictionary
        Coordinates of throat endpoints stored in Dict form. Can be accessed
        via the dict keys 'head' and 'tail'.

    Notes
    -----
    (1) This model should not be applied to true 2D networks. Use
    `circular_pores` model instead.

    (2) By default, this model assumes that throat centroid and pore
    coordinates are colinear. If that's not the case, such as in extracted
    networks, `throat_centroid` could be passed as an optional argument, and
    the model takes care of the rest.
    """
    
    network = target.project.network
    # throats = network.map_throats(throats=target.Ts, origin=target)
    throats = network.throats('all')
    xyz = network['pore.coords']
    cn = network['throat.conns'][throats]
    L = Center_to_center_length(target=target) + 1e-15
    Dt = network[throat_diameter][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    L1 = _sp.zeros_like(L)
    L2 = _sp.zeros_like(L)
    # Handle the case where Dt > Dp
    mask = Dt > D1
    L1[mask] = 0.5 * D1[mask]
    L1[~mask] = _sp.sqrt(D1[~mask]**2 - Dt[~mask]**2) / 2
    mask = Dt > D2
    L2[mask] = 0.5 * D2[mask]
    L2[~mask] = _sp.sqrt(D2[~mask]**2 - Dt[~mask]**2) / 2
    # Handle non-colinear pores and throat centroids
    try:
        TC = network[throat_centroid][throats]
        LP1T = _sp.linalg.norm(TC - xyz[cn[:, 0]], axis=1) + 1e-15
        LP2T = _sp.linalg.norm(TC - xyz[cn[:, 1]], axis=1) + 1e-15
        unit_vec_P1T = (TC - xyz[cn[:, 0]]) / LP1T[:, None]
        unit_vec_P2T = (TC - xyz[cn[:, 1]]) / LP2T[:, None]
    except KeyError:
        unit_vec_P1T = (xyz[cn[:, 1]] - xyz[cn[:, 0]]) / L[:, None]
        unit_vec_P2T = -1 * unit_vec_P1T
    # Find throat endpoints
    EP1 = xyz[cn[:, 0]] + L1[:, None] * unit_vec_P1T
    EP2 = xyz[cn[:, 1]] + L2[:, None] * unit_vec_P2T
    # Handle throats w/ overlapping pores
    L1 = (4*L**2 + D1**2 - D2**2) / (8*L)
    # L2 = (4*L**2 + D2**2 - D1**2) / (8*L)
    h = (2*_sp.sqrt(D1**2/4 - L1**2)).real
    overlap = L - 0.5 * (D1+D2) < 0
    mask = overlap & (Dt < h)
    EP1[mask] = EP2[mask] = (xyz[cn[:, 0]] + L1[:, None] * unit_vec_P1T)[mask]
    return {'head': EP1, 'tail': EP2}


def Shrink_throat_radius_pressure_fitting(target, 
                                          pore_diameter='pore.diameter',
                                          throat_diameter='throat.diameter'):
    r'''
    Function assigns throats with a diameter greater than the diameter of their connected
    pores (i.e. exceeding the upper physical limit of 1 to be considered a constriction)
    the minimum diameter of its two connected pores. Throats that do not adhere 
    to this are unadjusted.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
        
    pore_diameter : string
        Dictionary key of the pore diameter values.

    throat_diameter : string
        Dictionary key of the throat diameter values.

    Returns
    -------
    Dt : ndarray, float
        Array containing throat diameter, where throat with a diameter greater 
        than the diameter of their connected pores are assigned the minumum diameter
        of its two connecting pores.
    '''
        
    np.warnings.filterwarnings('ignore', category=RuntimeWarning)
    network = target.project.network
    throats = network.throats('internal')

    cn = network['throat.conns'][throats]
    Dt = network[throat_diameter][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    
    for idx in range(0,len(Dt)):
        
        if (Dt[idx] > D1[idx]) or (Dt[idx] > D2[idx]):
            Dt[idx] = np.min([D1[idx],D2[idx]])
    
    return Dt


def Shrink_throat_radius(target, 
                         Adjust_SA, 
                         Extracted_throat_diameter,
                         pore_diameter='pore.diameter',
                         throat_diameter='throat.diameter'):
        
    r'''
    Function assigns throats with a diameter greater than the diameter of their connected
    pores (i.e. exceeding the upper physical limit of 1 to be considered a constriction)
    the minimum diameter of its two connected pores. Throats that do not adhere 
    to this are unadjusted. Aditionally, it alters the surface area of the pores 
    to reflect these changes in the throat diameter. 

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    
    Adjust_SA: int 
        Indicator if the pore area should be adjusted to reflect changes in the throat diameter.
            If Adjust_SA == 0 -> do not alter the Surface area
            If Adjust_SA == 1 -> do alter the Surface area
        
    Extracted_throat_diameter: ndarray, float
        Values of the throat diameter attained by extraction of porous medium 
        with the SNOW algorithm.
    
    pore_diameter : string
        Dictionary key of the pore diameter values.

    throat_diameter : string
        Dictionary key of the throat diameter values.

    Returns
    -------
    Dt : ndarray, float
        Array containing throat diameter, where throat with a diameter greater 
        than the diameter of their connected pores are assigned the minumum diameter
        of its two connecting pores.
    '''        
    
    np.warnings.filterwarnings('ignore', category=RuntimeWarning)
    network = target.project.network
    throats = network.throats('internal')

    # Obtain throat and pore diameters of all pore-throat-pore elements
    cn = network['throat.conns'][throats]
    Dt = network[throat_diameter][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    
    # Geometrical throat area
    Throat_SA_original_copy = Extracted_throat_diameter * np.pi * network['throat.length'][throats].copy()

    for idx in range(0,len(Dt)):
        
        # Change throat diameter for the condition stated in function description
        if (Dt[idx] > D1[idx]) or (Dt[idx] > D2[idx]):
            Dt[idx] = np.min([D1[idx],D2[idx]])
        
        # Adjust surface area for altered throat diameter
        if Adjust_SA == 1:
            Throat_SA_original = Throat_SA_original_copy[idx]
            Throat_SA_adjusted = Dt[idx] * np.pi * network['throat.length'][idx]
            
            assign_0 = network['pore.surface_area'][cn[:, 0]][idx] - Throat_SA_original/2 + Throat_SA_adjusted/2
            assign_1 = network['pore.surface_area'][cn[:, 1]][idx] - Throat_SA_original/2 + Throat_SA_adjusted/2
            idx_0 = cn[:, 0][idx]
            idx_1 = cn[:, 1][idx]
            network['pore.surface_area'][idx_0] = assign_0
            network['pore.surface_area'][idx_1] = assign_1
            
    return Dt