# Import functions
import openpnm as op
import numpy as np
from matplotlib import rc
from matplotlib import pyplot as plt
from matplotlib import cm as cm
from matplotlib import colors
import openpyxl
import pandas as pd
import pdb                                                  
import Custom_network_plotting as cn

'''
Plot current density over specified axis:
    Height_width = 0 -> plot over heigth 
    Height_width = 1 -> plot over width
'''
Height_width = 0

# Type of flow field under investigation: 
Flow_field = 0

# Number of networks in thickness
n = 1

# Load in network and phase via csv file. Note that only the network can be loaded
# in with .csv files. We thus copy the phase properties to the network, and all
# will be saved under project_c.

# Input data
output_name = 'File_name_plot'
path = './input_plots/'
name_c = 'cathode_-1_0V_network0'
name_a = 'cathode_-1_0V_network0'

project_c = op.io.network_from_csv(filename = path + name_c + '.csv')
project_a = op.io.network_from_csv(filename = path + name_a + '.csv')

# Creating a new phase and network:
if Flow_field == 0:
    net_c = cn.Re_assign_network(project_c)
    net_a = cn.Re_assign_network(project_a)
if Flow_field == 1:
    net_c = cn.Re_assign_network_V2_IDFF(project_c)
    net_a = cn.Re_assign_network_V2_IDFF(project_a)
catholyte = cn.Re_assign_phase_V2(project_c, net_c, 'catholyte')
anolyte = cn.Re_assign_phase_V2(project_a, net_a, 'anolyte')

# Pressure
p = catholyte['pore.pressure'][net_c.pores('internal')]
p_max = p.max()
p_min = p.min()
print(p_max)

# Concentration
conc = catholyte['pore.concentration'][net_c.pores('internal')]
c_min = conc.min()
c_avg = conc.mean()
print(c_avg)

# General plot settings:
font_dict = {'sans-serif': 'Arial',
            'size': '12', 'weight': 'normal'}
rc('font', **font_dict) # Controls default text sizes
rc('axes', labelsize=18, linewidth=2, labelweight='normal')
rc('xtick', labelsize=18)
rc('xtick.major', width=2)
rc('ytick', labelsize=18)
rc('ytick.major', width=2)
rc('legend', fontsize=18)
rc('savefig', dpi=600)

# For the IDFF configuration we use the following axes. The current density will
# plotted over the width (the dimension (x) corresponding to the sides of the electrode)
def find_coords(net): 
    x = net['pore.coords'][net.pores('internal'), 1]
    y = net['pore.coords'][net.pores('internal'), 2]
    
    # correct network coordinates
    y = y+abs(y.min())
    z = net['pore.coords'][net.pores('internal'), 0]
    return x, y, z

# Get coordinates of anode and cathode
x_c, y_c, z_c = find_coords(net=net_c)
x_a, y_a, z_a = find_coords(net=net_a)

# correct network coordinates.
z_c = z_c - z_c.min(); 
z_a = z_a - z_a.min(); 

# Find width, length and thickness of anode and cathode
w_c = x_c.max() - x_c.min()
l_c = y_c.max() - y_c.min() 
t_c = z_c.max() - z_c.min()
w_a = x_a.max() - x_a.min()
l_a = y_a.max() - y_a.min()
t_a = z_a.max() - z_a.min() 
max_length = min(l_c, l_a)

# Define the performance indicators
# Current field [A]
i_c = abs(catholyte['pore.current'][net_c.pores('internal')])
i_a = abs(anolyte['pore.current'][net_a.pores('internal')])

# Convert units
gridsize = 16

# Current distribition in 1D:
def current_distribution(bin_width, z, i, i_cm3, pores):     
    z_2 = z 
    norm_i = i 
    norm_icm3 = i_cm3
    total_i_cm3 = i_cm3.sum()
    
    # Set up iterating variables
    max_z = int(max(np.ceil(z_2/bin_width)*bin_width))
    z_range = range(0, max_z, bin_width) 
    step = max_z / bin_width
    norm_z_range = np.linspace(0, 1, int(step))
    norm_ii = np.zeros(len(z_range)) 
    norm_iicm3 = np.zeros(len(z_range))
    bins = np.zeros(len(z_range))
    j = 0
    
    # Divide pores over bins
    for size_bin in z_range:
        for pore in pores:
            if size_bin < z_2[pore] <= size_bin + bin_width:
                bins[j] += 1
                norm_ii[j] += norm_i[pore]
                norm_iicm3[j] += norm_icm3[pore]
        j += 1
        
    cd = {}
    cd['z_range'] = z_range
    cd['norm_z_range'] = norm_z_range 
    cd['bins'] = bins
    cd['norm_ii'] = norm_ii
    cd['norm_iicm3'] = norm_iicm3 / total_i_cm3
    return cd

# Input parameters
i_cm3=(i_c)/(w_c*l_c*t_c)

# Current distribution
z_c_final = z_c*1e6
cd = current_distribution(bin_width = 5, z=z_c_final, i=i_c, i_cm3=i_cm3, pores=net_c.pores('internal'))
norm_z = z_c_final / z_c_final.max()

if Height_width == 0:
    # Plot length-thickness hexagonal binning plots for the current field
    fig, ax = plt.subplots(figsize=[6*n,8])
    hb = ax.hexbin(z_c_final, y_c*1e6, i_cm3, gridsize=(gridsize*n+n-1,9), cmap='viridis', edgecolors=None, bins='log', reduce_C_function=np.sum)
    cb = fig.colorbar(hb, ax=ax, orientation="horizontal",pad=0.18)
    ax.set_xlabel('Thickness [μm]')
    ax.set_ylabel('Length [μm]')
    ax.xaxis.set_ticks([0, z_c_final.max()/2, z_c_final.max()])
    cb.set_label('Current [A cm$^-$$^3$]')
    cb.mappable.set_clim(0.9e3, 1.1e4)
    # cb.set_ticklabels([1e4, 1e5])

    ax2 = ax.twinx()
    ax2.plot(cd['z_range'], cd['norm_iicm3'], 'k', linewidth=4)
    ax2.set_xlabel('Thickness [-]')
    ax2.set_ylabel('Normalized current [-]', rotation=270)
    ax2.xaxis.set_ticks([0, z_c_final.max()/2, z_c_final.max()])
    ax2.yaxis.set_ticks([])
    #plt.ticklabel_format(axis="y", style="sci", scilimits = (0,0))
    ax2.yaxis.set_label_coords(1.04, .5)
    plt.savefig(fname='.\\plots\\current_field' + output_name + '.tif', bbox_inches='tight')
    plt.show()
    
elif Height_width == 1:    
    # Plot length-thickness hexagonal binning plots for the current field
    fig, ax = plt.subplots(figsize=[6*n,8])
    hb = ax.hexbin(z_c_final, x_c*1e6, i_cm3, gridsize=(gridsize*n+n-1,9), cmap='viridis', edgecolors=None, bins='log', reduce_C_function=np.sum)
    cb = fig.colorbar(hb, ax=ax, orientation="horizontal",pad=0.18)
    ax.set_xlabel('Thickness [μm]')
    ax.set_ylabel('Width [μm]')
    ax.xaxis.set_ticks([0, z_c_final.max()/2, z_c_final.max()])
    cb.set_label('Current [A cm$^-$$^3$]')
    cb.mappable.set_clim(0.9e3, 1.1e4)
    # cb.set_ticklabels([1e4, 1e5])
    
    ax2 = ax.twinx()
    ax2.plot(cd['z_range'], cd['norm_iicm3'], 'k', linewidth=4)
    ax2.set_xlabel('Thickness [-]')
    ax2.set_ylabel('Normalized current [-]', rotation=270)
    ax2.xaxis.set_ticks([0, z_c_final.max()/2, z_c_final.max()])
    ax2.yaxis.set_ticks([])
    #plt.ticklabel_format(axis="y", style="sci", scilimits = (0,0))
    ax2.yaxis.set_label_coords(1.05, .5)
    plt.savefig(fname='.\\plots\\current_field' + output_name + '.tif', bbox_inches='tight')
    plt.show()

# Current profile:
i_tot_array = cd['norm_iicm3']*i_cm3.sum()
factor = int(len(i_tot_array) / n)

if n == 1:
    I_network_1 = i_tot_array[0:].sum()
    print(I_network_1)
    
if n == 2:
    I_network_1 = i_tot_array[factor:].sum()
    I_network_2 = i_tot_array[0:factor].sum()
    I_n_1_percentage = I_network_1 / i_cm3.sum() * 100
    I_n_2_percentage = I_network_2 / i_cm3.sum() * 100
    print(i_cm3.sum())
    print(I_n_1_percentage)
    print(I_n_2_percentage)

if n == 3:
    I_network_1 = i_tot_array[factor*2:].sum()
    I_network_2 = i_tot_array[factor:factor*2].sum()
    I_network_3 = i_tot_array[0:factor].sum()
    I_n_1_percentage = I_network_1 / i_cm3.sum() * 100
    I_n_2_percentage = I_network_2 / i_cm3.sum() * 100
    I_n_3_percentage = I_network_3 / i_cm3.sum() * 100
    print(i_cm3.sum())
    print(I_n_1_percentage)
    print(I_n_2_percentage)
    print(I_n_3_percentage)