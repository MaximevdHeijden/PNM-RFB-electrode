# -*- coding: utf-8 -*-
"""
Custom functions phase properties

OpenPNM automatically assigns water models to update anolyte and catholyte 
properties. This model overwrites these models and asssigns truly constants
"""

def custom_viscosity(target, parameter_script , prop):
    mu = parameter_script[prop]
    return mu    


def custom_density(target, parameter_script , prop):
    rho = parameter_script[prop]
    return rho    


def custom_diffusivity(target, parameter_script , prop):
    Diff = parameter_script[prop]
    return Diff    


def custom_electrical_conductivity(target, parameter_script , prop):
    EC = parameter_script[prop]
    return EC    


def custom_electrical_conductivity_fit(target, conductivity_fitted):
    EC = conductivity_fitted
    return EC    