"""
Creating a new network to correctly re-assign the properties and the labels 
of the network.
"""

import openpnm as op
import numpy as np

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
    # network_new['pore.equivalent_diameter'] = network['pore.equivalent_diameter']
    # network_new['pore.extended_diameter'] = network['pore.extended_diameter']
    # network_new['pore.inscribed_diameter'] = network['pore.inscribed_diameter']
    network_new['pore.surface_area'] = network['pore.surface_area']
    network_new['pore.volume'] = network['pore.volume']
    network_new['throat.area'] = network['throat.area']
    network_new['throat.centroid'] = network['throat.centroid']
    network_new['throat.conduit_lengths.pore1'] = network['throat.conduit_lengths.pore1']
    network_new['throat.conduit_lengths.pore2'] = network['throat.conduit_lengths.pore2']
    network_new['throat.conduit_lengths.throat'] = network['throat.conduit_lengths.throat']
    network_new['throat.diameter'] = network['throat.diameter']
    # network_new['throat.direct_length'] = network['throat.direct_length']
    network_new['throat.endpoints.head'] = network['throat.endpoints.head']
    network_new['throat.endpoints.tail'] = network['throat.endpoints.tail']
    # network_new['throat.equivalent_diameter'] = network['throat.equivalent_diameter']
    # network_new['throat.inscribed_diameter'] = network['throat.inscribed_diameter']
    # network_new['throat.length'] = network['throat.length']
    # network_new['throat.perimeter'] = network['throat.perimeter']
    # network_new['throat.total_length'] = network['throat.total_length']
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
    network_new['pore.surface_area'] = network['pore.surface_area']
    network_new['pore.volume'] = network['pore.volume']
    network_new['throat.area'] = network['throat.area']
    network_new['throat.centroid'] = network['throat.centroid']
    network_new['throat.conduit_lengths.pore1'] = network['throat.conduit_lengths.pore1']
    network_new['throat.conduit_lengths.pore2'] = network['throat.conduit_lengths.pore2']
    network_new['throat.conduit_lengths.throat'] = network['throat.conduit_lengths.throat']
    network_new['throat.diameter'] = network['throat.diameter']
    network_new['throat.endpoints.head'] = network['throat.endpoints.head']
    network_new['throat.endpoints.tail'] = network['throat.endpoints.tail']
    network_new['throat.length'] = network['throat.length']
    #network_new['throat.perimeter'] = network['throat.perimeter']
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


def Re_assign_network_V2_IDFF(network):
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
    network_new['pore.surface_area'] = network['pore.surface_area']
    network_new['pore.volume'] = network['pore.volume']
    network_new['throat.area'] = network['throat.area']
    network_new['throat.centroid'] = network['throat.centroid']
    network_new['throat.conduit_lengths.pore1'] = network['throat.conduit_lengths.pore1']
    network_new['throat.conduit_lengths.pore2'] = network['throat.conduit_lengths.pore2']
    network_new['throat.conduit_lengths.throat'] = network['throat.conduit_lengths.throat']
    network_new['throat.diameter'] = network['throat.diameter']
    network_new['throat.endpoints.head'] = network['throat.endpoints.head']
    network_new['throat.endpoints.tail'] = network['throat.endpoints.tail']
    network_new['throat.length'] = network['throat.length']
    #network_new['throat.perimeter'] = network['throat.perimeter']
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
    
    Flow_inlet_pores = network.pores('pore.flow_inlet')
    network_new.set_label(label = 'flow_inlet', pores = Flow_inlet_pores)
    Flow_outlet_pores = network.pores('pore.flow_outlet')
    network_new.set_label(label = 'flow_outlet', pores = Flow_outlet_pores)
    Membrane_pores = network.pores('pore.membrane')
    network_new.set_label(label = 'membrane', pores = Membrane_pores)
        
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


def Re_assign_phase_V2(network, net_new, phase_name):
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
    phase_new = op.phase.Water(network=net_new, name = phase_name)
    
    # Assign the properties
    phase_new['pore.activation_overpotential'] = network['pore.activation_overpotential']
    phase_new['pore.concentration'] = network['pore.concentration']
    phase_new['pore.concentration_overpotential'] = network['pore.concentration_overpotential']
    phase_new['pore.current'] = network['pore.current']
    phase_new['pore.density'] = network['pore.density']
    phase_new['pore.diffusivity'] = network['pore.diffusivity']
    phase_new['pore.electrical_conductivity'] = network['pore.electrical_conductivity']
    phase_new['pore.molar_density'] = network['pore.molar_density']
    phase_new['pore.molecular_weight'] = network['pore.molecular_weight']
    phase_new['pore.pressure'] = network['pore.pressure']
    phase_new['pore.surface_tension'] = network['pore.surface_tension']
    phase_new['pore.temperature'] = network['pore.temperature']
    phase_new['pore.thermal_conductivity'] = network['pore.thermal_conductivity']
    phase_new['pore.vapor_pressure'] = network['pore.vapor_pressure']
    phase_new['pore.viscosity'] = network['pore.viscosity']
    phase_new['pore.voltage'] = network['pore.voltage']
    
    return phase_new