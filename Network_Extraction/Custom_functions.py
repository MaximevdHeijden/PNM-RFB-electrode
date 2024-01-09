"""
Custom function module containing the following functions:
    1. Flow_shape_factors_ball_and_stick: Computing the shape (size) factors
       for hydraulic conductance.
    2. Hydraulic_conductance_Hagen_Poiseuille: Computing the Hydraulic 
       conductance used in Stokesflow.
    3. Poisson_shape_factors_ball_and_stick: conduit shape factors for throat 
       conductance associated with diffusion-like physics.
    4. Diffusive_conductance_mixed_diffusion: computing diffusive conductance 
       used in Stokesflow.
    5. Advection_diffusion: Calculate the advective-diffusive conductance of
       conduits in network. 

These functions are retrieved from older OpenPNM versions (such as OpenPNM version 2.6)
"""

from scipy import pi as _pi
from scipy import arctanh as _atanh
import scipy as _sp
import numpy as _np
import scipy.constants as const

def Flow_shape_factors_ball_and_stick(target, pore_area='pore.area',
                   throat_area='throat.area',
                   pore_diameter='pore.diameter',
                   throat_diameter='throat.diameter',
                   conduit_lengths='throat.conduit_lengths'):
    r"""
    Calculate conduit shape factors for hydraulic conductance, assuming
    pores and throats are spheres (balls) and constant cross-section
    cylinders (sticks).

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_diameter : string
        Dictionary key of the pore diameter values

    throat_diameter : string
        Dictionary key of the throat diameter values

    conduit_lengths : string
        Dictionary key of the conduit lengths' values

    Returns
    -------
    SF : dictionary
        Dictionary containing conduit shape factors to be used in hagen-
        poiseuille hydraulic conductance model. Shape factors are accessible
        via the keys: 'pore1', 'pore2' and 'throat'.

    Notes
    -----
    (1) This model accounts for the variable cross-section area in spheres.

    (2) WARNING: This model could break if `conduit_lengths` does not
    correspond to an actual ball and stick! Example: pore length is greater
    than pore radius 

    References
    ----------
    Akbari, M., Sinton, D., & Bahrami, M. (2011). Viscous flow in variable
    cross-section microchannels of arbitrary shapes. International Journal of
    Heat and Mass Transfer, 54(17-18), 3970-3978.
    """
    
    _np.warnings.filterwarnings('ignore', category=RuntimeWarning)
    network = target.project.network
    #throats = network.map_throats(throats=target.Ts, origin=target)
    throats = network.throats('all')
    cn = network['throat.conns'][throats]
    # Get pore diameter
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    Dt = network[throat_diameter][throats]
    # Get conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    # Get pore/throat baseline areas (the one used in generic conductance)
    A1 = network[pore_area][cn[:, 0]]
    A2 = network[pore_area][cn[:, 1]]
    At = network[throat_area][throats]
    # Preallocating F, SF
    # F is INTEGRAL(1/A^2) dx , x : 0 --> L
    F1, F2, Ft = _sp.zeros((3, len(Lt)))
    SF1, SF2, SFt = _sp.ones((3, len(Lt)))
    # Setting SF to 1 when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    SF1[~m1] = SF2[~m2] = SFt[~mt] = 1
    # Handle the case where Dt >= Dp
    M1, M2 = [(Di <= Dt) & mi for Di, mi in zip([D1, D2], [m1, m2])]
    F1[M1] = 16/3 * (L1*(D1**2 + D1*Dt + Dt**2) / (D1**3 * Dt**3 * _pi**2))[M1]
    F2[M2] = 16/3 * (L2*(D2**2 + D2*Dt + Dt**2) / (D2**3 * Dt**3 * _pi**2))[M2]
    # Handle the rest (true balls and sticks)
    N1, N2 = [(Di > Dt) & mi for Di, mi in zip([D1, D2], [m1, m2])]
    F1[N1] = (4/(D1**3*_pi**2) * ((2*D1*L1) / (D1**2-4*L1**2) + _atanh(2*L1/D1)))[N1]
    F2[N2] = (4/(D2**3*_pi**2) * ((2*D2*L2) / (D2**2-4*L2**2) + _atanh(2*L2/D2)))[N2]
    Ft[mt] = (Lt / At**2)[mt]
    # Calculate conduit shape factors
    SF1[m1] = (L1 / (A1**2 * F1))[m1]
    SF2[m2] = (L2 / (A2**2 * F2))[m2]
    SFt[mt] = (Lt / (At**2 * Ft))[mt]
    _np.warnings.filterwarnings('default', category=RuntimeWarning)
    
    return {'pore1': SF1, 'throat': SFt, 'pore2': SF2}


def Hydraulic_conductance_Hagen_Poiseuille(target,
                     pore_area='pore.area',
                     throat_area='throat.area',
                     pore_viscosity='pore.viscosity',
                     throat_viscosity='throat.viscosity',
                     conduit_lengths='throat.conduit_lengths',
                     conduit_shape_factors='throat.flow_shape_factors'):
    r"""
    Calculate the hydraulic conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_viscosity : string
        Dictionary key of the pore viscosity values

    throat_viscosity : string
        Dictionary key of the throat viscosity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    Returns
    -------
    g : ndarray
        Array containing hydraulic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper conduit_shape_factors argument.

    (4) shape_factor depends on the physics of the problem, i.e. diffusion-like
    processes and fluid flow need different shape factors.
    """
    
    network = target.project.network
    #throats = network.map_throats(throats=target.Ts, origin=target)
    throats = network.throats('all')
    #phase = target.project.find_phase(target)
    phase = network.project.phases.copy()[0]
    cn = network['throat.conns'][throats]
    # Getting equivalent areas
    A1 = network[pore_area][cn[:, 0]]
    At = network[throat_area][throats]
    A2 = network[pore_area][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    # Preallocating g
    g1, g2, gt = _sp.zeros((3, len(Lt)))
    # Setting g to inf when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    g1[~m1] = g2[~m2] = gt[~mt] = _sp.inf
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors+'.pore1'][throats]
        SFt = phase[conduit_shape_factors+'.throat'][throats]
        SF2 = phase[conduit_shape_factors+'.pore2'][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Interpolate pore phase property values to throats
    try:
        Dt = phase[throat_viscosity][throats]
    except KeyError:
        Dt = phase.interpolate_data(propname=pore_viscosity)[throats]
    try:
        D1 = phase[pore_viscosity][cn[:, 0]]
        D2 = phase[pore_viscosity][cn[:, 1]]
    except KeyError:
        D1 = phase.interpolate_data(propname=throat_viscosity)[cn[:, 0]]
        D2 = phase.interpolate_data(propname=throat_viscosity)[cn[:, 1]]
    # Find g for half of pore 1, throat, and half of pore 2
    pi = _sp.pi
    g1[m1] = A1[m1]**2 / (8*pi*D1*L1)[m1]
    g2[m2] = A2[m2]**2 / (8*pi*D2*L2)[m2]
    gt[mt] = At[mt]**2 / (8*pi*Dt*Lt)[mt]
    # Apply shape factors and calculate the final conductance
    
    return (1/gt/SFt + 1/g1/SF1 + 1/g2/SF2)**(-1)


def Poisson_shape_factors_ball_and_stick(target, pore_area='pore.area',
                   throat_area='throat.area',
                   pore_diameter='pore.diameter',
                   throat_diameter='throat.diameter',
                   conduit_lengths='throat.conduit_lengths'):
    r"""
    Calculate conduit shape factors for throat conductance associated with
    diffusion-like physics (ex. thermal/diffusive/electrical conductance),
    assuming pores and throats are spheres (balls) and constant cross-section
    cylinders (sticks).

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_diameter : string
        Dictionary key of the pore diameter values

    throat_diameter : string
        Dictionary key of the throat diameter values

    conduit_lengths : string
        Dictionary key of the conduit lengths' values

    Returns
    -------
    SF : dictionary
        Dictionary containing conduit shape factors to be used in conductance
        models associated with diffusion-like physics. Shape factors are
        accessible via the keys: 'pore1', 'pore2' and 'throat'.

    Notes
    -----
    (1) This model accounts for the variable cross-section area in spheres.

    (2) WARNING: This model could break if `conduit_lengths` does not
    correspond to an actual ball and stick! Example: pore length is greater
    than pore radius 
    """
    
    _np.warnings.filterwarnings('ignore', category=RuntimeWarning)
    network = target.project.network
    #throats = network.map_throats(throats=target.Ts, origin=target)
    throats = network.throats('all')
    cn = network['throat.conns'][throats]
    # Get pore diameter
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    Dt = network[throat_diameter][throats]
    # Get conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    # Get pore/throat baseline areas (the one used in generic conductance)
    A1 = network[pore_area][cn[:, 0]]
    A2 = network[pore_area][cn[:, 1]]
    At = network[throat_area][throats]
    # Preallocating F, SF
    # F is INTEGRAL(1/A) dx , x : 0 --> L
    F1, F2, Ft = _sp.zeros((3, len(Lt)))
    SF1, SF2, SFt = _sp.ones((3, len(Lt)))
    # Setting SF to 1 when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    SF1[~m1] = SF2[~m2] = SFt[~mt] = 1
    # Handle the case where Dt >= Dp
    M1, M2 = [(Di <= Dt) & mi for Di, mi in zip([D1, D2], [m1, m2])]
    F1[M1] = (4*L1/(D1*Dt*_pi))[M1]
    F2[M2] = (4*L2/(D2*Dt*_pi))[M2]
    # Handle the rest (true balls and sticks)
    N1, N2 = [(Di > Dt) & mi for Di, mi in zip([D1, D2], [m1, m2])]
    F1[N1] = (2/(D1*_pi) * _atanh(2*L1/D1))[N1]
    F2[N2] = (2/(D2*_pi) * _atanh(2*L2/D2))[N2]
    Ft[mt] = (Lt/At)[mt]
    # Calculate conduit shape factors
    SF1[m1] = (L1 / (A1*F1))[m1]
    SF2[m2] = (L2 / (A2*F2))[m2]
    SFt[mt] = (Lt / (At*Ft))[mt]
    _np.warnings.filterwarnings('default', category=RuntimeWarning)
    
    return {'pore1': SF1, 'throat': SFt, 'pore2': SF2}


def Diffusive_conductance_mixed_diffusion(
        target,
        pore_area='pore.area',
        throat_area='throat.area',
        pore_diameter='pore.diameter',
        throat_diameter='throat.diameter',
        pore_diffusivity='pore.diffusivity',
        pore_temperature='pore.temperature',
        molecular_weight='pore.molecular_weight',
        throat_diffusivity='throat.diffusivity',
        conduit_lengths='throat.conduit_lengths',
        conduit_shape_factors='throat.poisson_shape_factors'):
    r"""
    Calculate the diffusive conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ), assuming Knudsen
    diffusivity. See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_diffusivity : string
        Dictionary key of the pore diffusivity values

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    Returns
    -------
    g : ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (0) This function is ONLY suitable for dilute mixtures and NOT those with
    concentrated species.

    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper conduit_shape_factors argument.

    (4) shape_factor depends on the physics of the problem, i.e. diffusion-like
    processes and fluid flow need different shape factors.
    """
    
    network = target.project.network
    #throats = network.map_throats(throats=target.Ts, origin=target)
    throats = network.throats('all')
    #phase = target.project.find_phase(target)
    phase = network.project.phases.copy()[0]
    cn = network['throat.conns'][throats]
    # Getting equivalent areas
    A1, A2 = network[pore_area][cn].T
    At = network[throat_area][throats]
    # Getting conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors+'.pore1'][throats]
        SFt = phase[conduit_shape_factors+'.throat'][throats]
        SF2 = phase[conduit_shape_factors+'.pore2'][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Interpolate pore phase property values to throats
    D1, D2 = phase[pore_diffusivity][cn].T
    #Dt = phase.interpolate_data(propname=pore_diffusivity)[throats]
    Dt = phase['throat.diffusivity']
    # Calculating Knudsen diffusivity
    d1, d2 = network[pore_diameter][cn].T
    dt = network[throat_diameter][throats]
    MW1, MW2 = phase[molecular_weight][cn].T
    #MWt = phase.interpolate_data(propname=molecular_weight)[throats]
    MWt = phase['throat.molecular_weight']      # Interpolating pore -> throat
    T1, T2 = phase[pore_temperature][cn].T
    #Tt = phase.interpolate_data(pore_temperature)
    Tt = phase['throat.temperature']            # Interpolating pore -> throat 
    DK1 = d1/3 * (8*const.R*T1/const.pi/MW1)**0.5
    DK2 = d2/3 * (8*const.R*T2/const.pi/MW2)**0.5
    DKt = dt/3 * (8*const.R*Tt/const.pi/MWt)**0.5
    # Calculate mixed diffusivity
    D1e = (1/DK1 + 1/D1)**(-1)
    D2e = (1/DK2 + 1/D2)**(-1)
    Dte = (1/DKt + 1/Dt)**(-1)
    # Find g for half of pore 1, throat, and half of pore 2
    g1 = (D1e * A1) / L1
    g2 = (D2e * A2) / L2
    gt = (Dte * At) / Lt
    # Apply shape factors and calculate the final conductance
    
    return (1/gt/SFt + 1/g1/SF1 + 1/g2/SF2)**(-1)


def Advection_diffusion(target,
           conduit_lengths='throat.conduit_lengths',
           pore_pressure='pore.pressure',
           throat_hydraulic_conductance='throat.hydraulic_conductance',
           throat_diffusive_conductance='throat.diffusive_conductance',
           s_scheme='powerlaw'):
    r"""
    Calculate the advective-diffusive conductance of conduits in network, where
    a conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    conduit_lengths : string
        Dictionary key of the conduit length values

    pore_pressure : string
        Dictionary key of the pore pressure values

   throat_hydraulic_conductance : string
       Dictionary key of the throat hydraulic conductance values

   throat_diffusive_conductance : string
       Dictionary key of the throat diffusive conductance values

   s_scheme : string
       Name of the space discretization scheme to use

    Returns
    -------
    g : ndarray
        Array containing advective-diffusive conductance values for conduits in
        the geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper conduit_shape_factors argument when
    computig the diffusive and hydraulic conductances.

    (4) shape_factor depends on the physics of the problem, i.e. diffusion-like
    processes and fluid flow need different shape factors.
    """
    
    network = target.project.network
    #throats = network.map_throats(throats=target.Ts, origin=target)
    throats = network.throats('all')
    #phase = target.project.find_phase(target)
    phase = network.project.phases.copy()[0]
    cn = network['throat.conns'][throats]
    # Getting conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    # Preallocating g
    g1, g2, gt = _sp.zeros((3, len(Lt)))
    # Setting g to inf when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    g1[~m1] = g2[~m2] = gt[~mt] = _sp.inf
    # Find g for half of pore 1, throat, and half of pore 2
    P = phase[pore_pressure]
    gh = phase[throat_hydraulic_conductance][throats]
    gd = phase[throat_diffusive_conductance][throats]
    gd = _sp.tile(gd, 2)

    Qij = -gh*_sp.diff(P[cn], axis=1).squeeze()
    Qij = _sp.append(Qij, -Qij)

    Peij = Qij / gd
    Peij[(Peij < 1e-10) & (Peij >= 0)] = 1e-10
    Peij[(Peij > -1e-10) & (Peij <= 0)] = -1e-10

    # Export Peclet values (half only since Peij = -Peji)
    phase['throat.peclet.ad'] = _sp.nan
    phase['throat.peclet.ad'][throats] = _sp.absolute(Peij[0:len(Lt)])

    # Correct the flow rate
    Qij = Peij * gd

    if s_scheme == 'upwind':
        w = gd + _sp.maximum(0, -Qij)
    elif s_scheme == 'hybrid':
        w = _sp.maximum(0, _sp.maximum(-Qij, gd-Qij/2))
    elif s_scheme == 'powerlaw':
        w = gd * _sp.maximum(0, (1 - 0.1*_sp.absolute(Peij))**5) + \
            _sp.maximum(0, -Qij)
    elif s_scheme == 'exponential':
        w = -Qij / (1 - _sp.exp(Peij))
    else:
        raise Exception('Unrecognized discretization scheme: ' + s_scheme)
    w = _sp.reshape(w, (throats.size, 2), order='F')
    
    return w


def Electrical_conductance_series_resistors(target,
                     pore_area='pore.area',
                     throat_area='throat.area',
                     pore_conductivity='pore.electrical_conductivity',
                     throat_conductivity='throat.electrical_conductivity',
                     conduit_lengths='throat.conduit_lengths',
                     conduit_shape_factors='throat.poisson_shape_factors'):
    r"""
    Calculate the electrical conductance of conduits (throat.electrical_conductance) 
    in network, where a conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_conductivity : string
        Dictionary key of the pore thermal conductivity values

    throat_conductivity : string
        Dictionary key of the throat thermal conductivity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    Returns
    -------
    g : ndarray
        Array containing electrical conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper conduit_shape_factors argument.

    (4) shape_factor depends on the physics of the problem, i.e. diffusion-like
    processes and fluid flow need different shape factors.
    """
    
    network = target.project.network
    #throats = network.map_throats(throats=target.Ts, origin=target)
    throats = network.throats('all')
    #phase = target.project.find_phase(target)
    phase = network.project.phases.copy()[0]
    cn = network['throat.conns'][throats]
    # Getting equivalent areas
    A1 = network[pore_area][cn[:, 0]]
    At = network[throat_area][throats]
    A2 = network[pore_area][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    # Preallocating g
    g1, g2, gt = _sp.zeros((3, len(Lt)))
    # Setting g to inf when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    g1[~m1] = g2[~m2] = gt[~mt] = _sp.inf
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors+'.pore1'][throats]
        SFt = phase[conduit_shape_factors+'.throat'][throats]
        SF2 = phase[conduit_shape_factors+'.pore2'][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Interpolate pore phase property values to throats
    try:
        Dt = phase[throat_conductivity][throats]
    except KeyError:
        Dt = phase.interpolate_data(propname=pore_conductivity)[throats]
    try:
        D1 = phase[pore_conductivity][cn[:, 0]]
        D2 = phase[pore_conductivity][cn[:, 1]]
    except KeyError:
        D1 = phase.interpolate_data(propname=throat_conductivity)[cn[:, 0]]
        D2 = phase.interpolate_data(propname=throat_conductivity)[cn[:, 1]]
    # Find g for half of pore 1, throat, and half of pore 2
    g1[m1] = (D1*A1)[m1] / L1[m1]
    g2[m2] = (D2*A2)[m2] / L2[m2]
    gt[mt] = (Dt*At)[mt] / Lt[mt]
    # Apply shape factors and calculate the final conductance
    return (1/gt/SFt + 1/g1/SF1 + 1/g2/SF2)**(-1)


def custom_diffusivity(target , prop):
    Diff = prop
    return Diff    