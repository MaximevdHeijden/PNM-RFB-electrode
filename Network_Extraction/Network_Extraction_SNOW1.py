# Import functions:
import porespy as ps
import scipy as sp
import matplotlib.pyplot as plt
import openpnm as op
from skimage import io
from skimage import util
import time
import numpy as np
import imageio
import os

# Measure the time the programme takes to completes the snow network extraction
start_time = time.time()
op.Workspace().clear()

# Load in image
cwd = os.getcwd()                               # Curent Working directory
network_path = '\\input\\'                      # Folder with stored file
file_name ='Freudenberg'                        # Name of the file  
directory = cwd + network_path + file_name + '.tif'
im = io.imread(directory)

# SNOW algorithm is relatively fast up to 300x300x300 voxels
print('(x,y,z) shape:', im.shape)

# Specify the scanning voxel size of the binarized image
voxel_size = 1.65e-6                            # [m]

# SNOW algorithm accepts Binary image in the Boolean form with True’s as void phase and False’s as solid phase.
im_void = util.invert(im)                       # Invert void to True and solid to False
 
im_void1 = im_void > 0.5                        # Try the other way around too. 
 
# Perform snow extraction
# Accepts  Binary image in the Boolean form with True’s as void phase and False’s as solid phase.
snow_output = ps.networks.snow(im_void1,
                               voxel_size=voxel_size,
                               marching_cubes_area = False)
del im_void1
 
# Convert to PNM network
pn = op.network.GenericNetwork()
pn.update(snow_output)
proj = pn.project

# Depending on the image, the pore diameter type can be changed to get a more accurate extraction
proj.network['pore.diameter'] = proj.network['pore.inscribed_diameter']     # equivalent/inscribed
proj.network['throat.diameter'] = proj.network['throat.inscribed_diameter']

# Trim isolated pores
net = proj.network
print(net.check_network_health())
health = net.check_network_health()
trim_these_pores = health['trim_pores']
op.topotools.trim(network=net, pores=trim_these_pores)

# Generate an image that can overlay with the network in Paraview
#im = ps.tools.align_image_with_openpnm(im)     # Align the image with the coords used in openPNM

# Save data
network_path = '\\output\\'
directory = cwd + network_path + file_name
op.io.CSV.save(network = pn, filename = directory)
proj.save_project(directory)
proj.export_data(filename=directory, filetype='vtk')
imageio.mimsave(directory+ '_overlay_image.tif', np.array(im, dtype=int))

#Saving to a pandas file.
Network_Snow1 = op.io.Pandas.to_dataframe(network = pn) 

#Saving to a .CSV file
op.io.CSV.save(network = pn, filename = directory)

print('--- SNOW-extraction took: %s seconds ---' % (time.time() - start_time))