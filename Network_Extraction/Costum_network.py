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


def Rotating_network(network, angle, shift, shift_value):
    r'''
    Function rotates the network specified angle around x-axes and re-assigns 
    the correct labels from a positive z-direction reference. 

    Parameters
    ----------
    network : pn.Network object
    angle : Radial angle of rotation
    shift: shifting the network with in y-direction with specified value such that its boundary plane coincides with the origin (e.g. y_min = 0)
    Returns
    -------
    network : TYPE
        DESCRIPTION.
    '''
    
    # Copying "old" (original coordinates)
    coords_x = network['pore.coords'][:,0]
    coords_y = network['pore.coords'][:,1]
    coords_z = network['pore.coords'][:,2]
    
    # Shifting the coordinates 90 degrees around the x-axis
    # angle = np.pi/2
    x_new = coords_x
    y_new = coords_y*np.cos(angle) - coords_z*np.sin(angle)
    z_new = coords_y*np.sin(angle) + coords_z*np.cos(angle)
    
    # Assigning new coordinates
    network['pore.coords'][:,0] = x_new
    network['pore.coords'][:,1] = y_new
    network['pore.coords'][:,2] = z_new
    
    if shift == 1:
        network['pore.coords'][:,1] = network['pore.coords'][:,1] + shift_value
    
    # Re-assigning the labels
    back_pores = network.pores('back')      # Currently top
    front_pores = network.pores('front')    # Currently bottom
    top_pores =  network.pores('top')       # Currently front
    bottom_pores = network.pores('bottom')  # Currently back
    
    network.set_label('pore.back', pores=bottom_pores,  mode='overwrite')
    network.set_label('pore.front', pores=top_pores,  mode='overwrite')
    network.set_label('pore.top', pores=back_pores,  mode='overwrite')
    network.set_label('pore.bottom', pores=front_pores,  mode='overwrite')
    
    return network