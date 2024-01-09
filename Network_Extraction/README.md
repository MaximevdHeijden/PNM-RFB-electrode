# Network Extraction 
## Using OpenPNM v3.0.0 for OpenPNM v2.2.0 extracted modules

This document will discuss the extraction of pore networks from tif images and the used physics by importing OpenPNM networks created in OpenPNM v2.2.0 into OpenPNM v3.0.0.\
Networks are extracted in the script “Network_Extraction_SNOW1” which uses the SNOW 1 function to extract the network instead of the more recent SNOW 2 function. The files pertaining to OpenPNM v2.2.0 are stored in the folder “output”, while the files pertaining to OpenPNM v3.0.0 are stored in the folder “output_version_3_0”. The extraction consists of the following steps:
1.  Place the to be extracted network (.tif) in the input folder.

2.	Set the OpenPNM to v2.2.0 using GitKraken and set Porespy to v1.3.1 (using pip install porespy==1.3.1 in Anaconda prompt). We use Porespy 1.3.1. as it is the last Porespy version which still contains the SNOW 1 function. Any higher Porespy version automatically employs SNOW 2, which gives slightly different extracted networks than SNOW 1. SNOW 1 operates with OpenPNM v2.2.0 but NOT with OpenPNM v3.0.0.\
In case you switch between OpenPNM v3.0.0 to OpenPNM v2.2.0 you often run into package competency issues including:\
    •	numba.errors. This error cannot be solved by installing packages, but by disabling the module. Follow the path openpnm->algorithms->InvasionPercolation and comment out the following three lines: "from numba.errors import NumbaPendingDeprecationWarning", "import warnings", and "warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)".\
    •	Install the following packages in Anaconda prompt using the "pip install" command: unyt, flatdict.\
    •   Comment the following lines in the mentioned scripts:\
        •	In porespy - tools - _funcs_: "from skimage.measure import marching_cubes_lewiner".\
        •	In porespy - metrics - _regionprops_: "from from skimage.measure import mesh_surface_area, marching_cubes_lewiner".\
        •	In openpnm - openpnm - materials - VoronoiFibers: "from scipy.stats import itemfreq".\
    •   Change in porespy – tools - _funcs_ in function randomize_colors "new_vals = sp.random.permutation(im_vals)" to "new_vals = np.random.permutation(im_vals)".

3.	Initialize the “Network_Extraction_SNOW1” script by setting the correct file name, voxel_size, type of diameter used for the extraction (inscribed/equivalent), and the name of the exported network.

4.	Run the extraction script. All output data is saved in the “output” folder. These include:\
    •   A .csv file containing the properties of the network (pore and throat diameter, throat connections, surface area, cross-sectional area etc.).\
    •   A .pnm file to be used solely for OpenPNM version 2.2.0 simulations.\
    •   A .vtp file for visualization in Paraview.\
    •   An .tif overlay image for visualization in Paraview.

5.	Changing the file data type. The extracted .pnm network created in step 4 can be used ONLY for OpenPNM v2.2.0. OpenPNM v3.0.0. cannot load in the .pnm network created in step 4 due to missing geometry and physics modules. A network which can be used in OpenPNM v3.0.0 is created from the .csv file created in step 4. This .csv is loaded into the script the script “Netwerk_Properties_V3_with_V2Data”, which outputs a .pnm file which can be used in OpenPNM v3.0.0.\
    The .csv file, however, needs to be altered before loading it into the script “Netwerk_Properties_V3_with_V2Data”. Only the first row of the .csv file needs to be altered (where the strings are stored). Press ctrl+F, select the replace functionality and enter  “network | net_01 | “ (**NOTE:** the space at the end is NOT a typo). Leave the “Replace with” bar blanc and select “Replace all”.

6.	Creating the network for OpenPNM v3.0.0 from the altered .csv file. First, change to OpenPNM v3.0.0 in the Gitkraken app, but leave the Porespy version unaltered. Then open the script “Netwerk_Properties_V3_with_V2Data” and change the input file name of the altered .csv file, the name of the original .tif file, and the voxel size.

7.	Optional: rotation of the network (only for interdigitated flow fields). Depending on the stitching of two networks in ImageJ, the network might require a rotation along the x-axis. Rotation is “turned” on in the script Network_Properties_V3_wth_V2Data by assigning the option “Shift” a value of 1. Note that with a rotation you also need to shift the network such that the origin of the network aligns with the origin of the “original” network (e.g., 1.65e-6; 1.65e-6; 1.65e-6) in case boundary pores are included into the network. The function “Rotating_network” performs: (1) A rotation of network for a specified angle around x-axes and re-assigns the correct labels from a positive z-direction reference. The standard is pi/2 (90 degrees). By visualization of the .vtp file in Paraview, you can check if your rotation was successful. (2) Shifting of the rotated network. This shifting is performed to align the network coordinates with the origin. Here the minimum X, Y and Z coordinates of the network should overlay with the origin (e.g., 1.65e-6; 1.65e-6; 1.65e-6) in case boundary pores are included into the network. The shifting currently only takes place for the y-direction. You find the correct value of the shift by trial and error, and checking the minimum coordinates of the network in all three directions (X,Y,Z). **NOTE:** if you don’t want to rotate the network during the extraction step, you should alter the assignment of boundary conditions in the main script for running the electrochemical interdigitated flow field simulations.

8.	Run the script (make sure that you do so in a new console/kernel to allow switching of the OpenPNM version in Spyder). All data is stored in the network path utput_version_3_0, and the following is outputted: a .csv file containing network properties (porosity, permeability, PSD etc), a .pnm network that can be used only in OpenPNM v3.0.0, and a .vtp file which can be used for visualization.

9.	The .pnm network created in step 8 can be used for simulations within OpenPNM v3.0.0. Important to note is that we still use the physics of OpenPNM v2.2.0. in OpenPNM v3.0.0 simulations. The reason hereof is that throat diameters exceeding the pore diameters crashes the spheres and cylinder physics module of OpenPNM v3.0.0 (an error warning is captured, and if the warning is commented out one obtains negative conduit lengths and thus negative transport physics, which are unphysical). These physics function are stored in the module “Costum_functions”.
