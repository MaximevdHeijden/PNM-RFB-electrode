# Pore network modeling scripts
## About
This repository contains pore network modeling approaches for the simulation of porous electrode microstructures for redox flow batteries which is for example used for the bottom-up design of electrode structures using a [genetic algorithm](https://github.com/MaximevdHeijden/GA-RFB-electrode). PNM-RFB-electrode is an open-source code written in Python using the open-source software OpenPNM.

## Installation
The scripts and functions in this folder use OpenPNM version 3.0.0. which can be installed using [OpenPNM documentation](https://openpnm.org/installation.html) and are written for Windows. To change the current version to OpenPNM version 3.0.0, [Gitkraken](https://www.gitkraken.com/) can be used. Before running the code, minor adjustments need to be made to the OpenPNM documentation, which can be found in the “READ Me – OpenPNM changes” file.
For the installation of OpenPNM from scratch:

1.	Download [Anaconda](https://www.anaconda.com/download/) and follow the installment procedure. The open-source Anaconda Distribution is the easiest way to perform Python/R data science and machine learning on Linux, Windows, and Mac OS X. Python is used as the programming language in this project and Spyder can be used as editor.

2.    Check pip installation (pip is the package manager of python) and install pipenv. First, open 'Anaconda prompt' via the Windows start menu. Check if pip is installed using:

      ```
      python -m pip --version
      ```

      If pip is installed, install pipenv using:

      ```
      python -m pip install --user pipenv
      ```

      You can get a warning that the installed packages are not on your PATH. Copy this path, go to your windows start menu and type in “Edit the system environment variables”, select “environment variables”, click on Path and select “Edit”. Then click “new” and paste the path you just copied.

3.    Install the latest version of [OpenPNM](https://openpnm.org/installation.html), following the “hard (but correct way)” will allow you to install the developing code:\
      •	First make a repository (folder) in which you want to save the code and copy the directory adress.\
      •	Open Anaconda Prompt and go to the folder by using the "cd" command.\
      •	Clone the repo via the command:
 
            ```
            git clone https://github.com/PMEAL/OpenPNM
            ```

      •	Enter the root folder of openPNM by using the "cd OpenPNM" command.\
      •	Enter the following command and note the space between the dot and e is not a typo:
       
            ```
            pip install --no-deps -e .
            ```

      •	Enter the following commands to install the openpnm dependencies:
       
            ```
            conda install --file requirements/conda.txt -c conda-forge
            pip install -e .
            ```

      •	Install [Gitkraken](https://www.gitkraken.com/). Gitkraken allows you to switch between different versions of OpenPNM. In Gitkraken open (select) the installed OpenPNM folder. Search for version 3.0.0 and select “Check out this commit”. Your editing program will automatically refer to this version.

4.    You can now run the OpenPNM code files. It could be that you get errors due to missing software packages, the most common ones are discussed below and the other errors can be solved using Google search:\
      •	You probably need to install the following packages: docrep, chemicals, pyamg, rich, thermo, transforms3d, pypardiso, and lmfit, using pip by entering the following commands in AnacondaPrompt:
 
            ```
            pip install docrep
            pip install chemicals
            pip install pyamg
            pip install thermo
            pip install transforms3d
            pip install pypardiso
            pip install lmfit
            ```

      •     It could be that the installed numpy version is not compatible with the code. This can be solved by changing the numpy version using:
       
            ```
            pip install numpy==1.23.5
            ```

      •	numba.errors. This error cannot be solved by installing packages, but by disabling the module. Follow the path openpnm->algorithms->InvasionPercolation and comment out "from numba.errors import NumbaPendingDeprecationWarning".\
      •     Boundary conditions issue. When the program can is ran, the following error can pop up: "Another boundary condition was detected in some of the locations recieved". The solution is to comment out an if-loop. Follow the path 
OpenPNM->Algorithm->generic_transport and go the function “def_set_BC” and comment out the following if-loop "if np.intersect1d(pores, BC_locs).size".

5.    Now the codes should work.

## Documentation
This repository contains several scripts that will be used when extracting networks, running the code, or for post-processing, including:
1.    Function scripts that are used but don't have to be altered to run your simulations:\
      •	**Costum_functions_phase:**\
            Additional phase functions used are defined in this script as OpenPNM automatically assigns water models to update the anolyte and catholyte properties.\
      •	**Costum_functions_transport:**\
            Additional functions that are retrieved from older OpenPNM versions but are required for the simulations to run.\
      •	**Costum_network:**\
            Functions required to create a new network to correctly re-assign the properties and the labels of the network.\
      •	**Costum_functions_km:**\
            Additional functions defined to run the model including the definition of the boundary pores, formulation of the Butler-Volmer equation, the convergence criteria, overpotentials, pressure drop, flowrate, and conductivity.\
      •	**Costum_functions_pressure_fitting:**\
            Additional functions including the definitions of the throat flowrate and hydraulic conductance, pore and throat velocities, and mass transfer coefficients.

2.    Input dictionary defining your operating conditions and electrolyte properties:\
**inputDictTEMPO:**\
The input parameters used in the algorithm are defined in this script. The input parameters can easily be adjusted based on the desired chemistry (we advise to experimentally determine most of the parameters) and operating conditions. Examples for input parameters can be found in the citations mentioned under the 'Cite' section (TEMPO and Iron) or in the [GA-RFB-electrode](https://github.com/MaximevdHeijden/GA-RFB-electrode) repository (vanadium).

3.    Folder "Network_Extraction" with the network extraction scripts:\
      In this folder, the scripts required for the network extraction are provided, including:\
      •	**Costum_functions:**\
            Additional functions that are retrieved from older OpenPNM versions but are required for the simulations to run.\
      •	**Costum_network:**\
            Functions required to create a new network to correctly re-assign the properties and the labels of the network.\
      •	**Costum_functions_pressure_fitting:**\
            Additional functions including the definitions of the throat flowrate and hydraulic conductance, pore and throat velocities, and mass transfer coefficients.\
      •     **Network_Extraction_SNOW1:**\
            This script is used to extract the pore network from tif images. **NOTE:** in order to run the network extraction scripts with SNOW1, several steps needs to be taken, which are explained in README-NetworkExtraction.\
      •	**Netwerk_Properties_V3_with_V2Data:**\
            This script computes the pore size distribution, porosity profiles, and permeability of the extracted network. Furthermore, .pnm and .vtk files are created which can be used for the pore network model simulations and for visualization in Paraview, respectively.

      Before running the code, the following folders need to be created:\
      •	**input:** this folder contains the networks that will be extracted and the sdjusted .csv file required for the network properties script.\
      •	**output:** this folder contains the output of the network extraction including .csv, .pnm, .vtk, and .tif files of the extracted network.\
      •	**output_version_3_0:** this folder contains the output of the network properties script including a .pnm, .vtk, and .xlsx files which can be used for pore network simulations, visualization, and contain network properties, respectively.

4.    Scripts for the pressure drop correction:\
      These scripts perform the role of hydraulic calibration and obtaining the fluid field (pressure field, velocity field) inside the electrode:\
      •	**Pressure_fitting_script:**\
            Pressure drop fitting script using the conduit conductance model to obtain the fitting values (contraction curvature, expansion curvature, and throat effective aspect ratio parameters) over a range of superficial inlet velocities that can be used in your simulations. This script performs the hydraulic callibration of the pore network model where the overall pressure drop is fitted to the experimentally obtained pressure drop using the least squares scheme.\
            **NOTE:** sometimes there is a convergence issue of the stokesflow algorithm. This often arises because a few pores are not converging (e.g., three pores) in the fluid field computation. These pores can be trimmed from the network, but this is not necessary as the results from the “unconverged” stokesflow were found to be the same as those where these few pores were trimmed from the network. The index of these pores can be located based on the change of pore pressure between consecutive iterations of the fluid field computation. This difference can be found by running the script in DEBUG MODE and set a debug point in the "OpenPNM – algorithms – transport function" script (altered version) in "_run_special" on the line "phase['pore.pressure'] = self.x.copy()". First run for a few iterations for the computation to stabilize (e.g., 5 times), then run the below commands in the console:

            ```
            net = self.project['net']                          
            internal_pores = net.pores('internal')
            Pres_old_int = x0[internal_pores]
            Pres_new_int = x_new[internal_pores]
            rel_pres_diff = np.abs((Pres_new_int - Pres_old_int)/Pres_new_int)
            print(np.where(rel_pres_diff > rel_pres_diff.max()/4))
            ```

        Here we compute the relative pore pressure difference between two consecutive iterations. We then print the index of those pores with a relative difference greater than a quarter of the maximum hereof (though this might require refinement on your case). Often the same pore index is found between consecutive iterations. Trimming of these few pores in the script “Flow_field_main” after loading in the original network then avoids any non-convergence warnings.\
      •	**Flow_field_main:**\
            This script obtains the fluid field of the pore network model. This script performs three roles: altering the network through changing the throat diameter with the fitted throat effective aspect ratio, outputting the altered network into the "Input_networks"folder, and outputting the fluid field at a specified inlet velocity in the folder "Input_Fluid_field".\
            The altered network and the fluid field can be subsequently loaded into the scripts used for fitting of the local mass transfer coefficient and used in the polarization runs.

5.    Scripts for the mass transfer coefficient correction:\
      These scripts fit the pre-exponential and exponential factors of a local mass transfer correlation to a global mass transfer correlation (obtained from e.g., literature or limiting current density experiments).\
      •	**Local_km_fitting_script_initialize_fluid_field:**\
            This fitting script allows you to load in a network and compute the fluid field within the fitting script with as only input a newtork from the "input" folder. Then, the pressure field is computed using OpenPNM physics (i.e., hydraulic conductance based on shape factors).\
      •	**Local_km_fitting_script_load_in_fluid_field:**\
            This fitting script allows you to load in a pre-ran fluid field in combination with a certain network (see the scripts for the pressure drop correction). This script requires you to fit the throat diameter and make an “altered network” as explained in the scripts for the pressure drop correction. To run this fitting script, you need to place the altered network in the folder “Input_networks” and the fluid field (at all velocities you want to run the mass transfer fitting at) in the older “Input_Fluid_field”.

6.    Polarization curve fitting scripts:\
      These scripts can be used to fit the pore network model to experimental polarization data, using a velocity independent or velocity dependent mass transfer correlation (see the scripts for the mass coefficient correlation). The pre-exponential and exponential factor of a possible LOCAL mass transfer correlation can be placed in the input dictionary. Fitting of the experimental data can occur through fitting the following parameters:
      1.	The conductivity
      2.	The mass transfer coefficient. **NOTE:** This is simply a factor multiplied by the calculated km. Application hereof to the velocity dependent mass transfer coefficient is then equal to fitting the pre-exponential factor of the correlation.
      3.	The membrane resistance. **NOTE:** The membrane resistivity is computed through Pouillet’s law. Fitting the conductivity then also (partially) fits the membrane resistance.

• **Fitting_script_OpenPNM:**\
This fitting script allows you to load in a network and compute the fluid field within the fitting script with as only input a newtork from the "input" folder. Then, the pressure field is computed using OpenPNM physics (i.e, hydraulic conductance based on shape factors).\
• **Fitting_script:**\
This fitting script allows you to load in a pre-ran fluid field in combination with a certain network (see the scripts for the pressure drop correction). This script requires you to fit the throat diameter and make an “altered network” as explained in the scripts for the pressure drop correction. To run this fitting script, you need to place the altered network in the folder “Input_networks” and the fluid field (at all velocities you want to run the mass transfer fitting at) in the older “Input_Fluid_field”.

7.    Main Pore Network Model script for simulating the redox flow battery electrodes:\
      These scripts can be used to create a polarization plot and to obtain the contributions of the individual overpotentials, using a velocity independent or velocity dependent mass transfer correlation (see the scripts for the mass coefficient correlation). The pre-exponential and exponential factor of a possible LOCAL mass transfer correlation can be placed in the input dictionary.\
      •	**PNM_main_V3_Final_OpenPNM:**\
            This script allows you to load in a network and compute the fluid field with as only input a newtork from the "input" folder. Then, the pressure field is computed using OpenPNM physics (i.e, hydraulic conductance based on shape factors).\
      •	**PNM_main_V3_Final:**\
            This script allows you to load in a pre-ran fluid field in combination with a certain network (see the scripts for the pressure drop correction). This script requires you to fit the throat diameter and make an “altered network” as explained in the scripts for the pressure drop correction. To run this script, you need to place the altered network in the folder “Input_networks” and the fluid field (at all velocities you want to run the mass transfer fitting at) in the older “Input_Fluid_field”.

8.    Visualization scripts:\
      The following scripts are used for postprocessing purposes. **NOTE:** Before running these scripts, several changes have to be made to the cathode/anode.csv files:
      1.    Search for "net." and replace with nothing (empty bar). Do the same for "catholyte_vtk." or "anolyte_vtk.", depending if it is the cathode or anode electrode.
      2.    Use the 'Text to Columns' function in Excel with separated by commas.
      3.    Delete the last columns (DT - EA, with repetitions of pore.velocity_components_OpenPNM_odd_network[0,1,2], pore.velocity_components_odd_network[0,1,2], pore.velocity_magnitude_OpenPNM_odd_network[0], and pore.velocity_magnitude_odd_network[0]) and save the new version of the .csv file.
      4.    Now the scripts should work.\
      •     **Costum_network_plotting:**\
            Functions required to create a new network to correctly re-assign the properties and the labels of the network.\
      •	**Current_plot:**\
            Script to visualize the current profiles through the electrode.\
      •	**Overpotential_plots:**\
            Scripts to visualize the overpotential contributions over the electrode thickness.\
      •	**Plotting_masstransfer:**\
            Visualization of the mass transfer coefficients of the internal pores for the cases with and without velocity-dependent mass transfer coefficients.\
      In addition to these visualization scripts, the outputted VTK files can be visualized in [Paraview](https://www.paraview.org/).

## Getting started
After installing OpenPNM version 3.0.0. and making the required adjustments, the pore network model can be run for the desired operating conditions and network properties, which need to be specified in the scripts.\
The following properties can be adjusted:\
•	Electrolyte chemistry: in the inputDict, the chemistry dependent parameters can be changed (e.g., diffusion coefficient, viscosity, concentration).\
•	Operation conditions: in inputDict, the operating conditions can be changes (e.g., total electrode lenght, membrane properties, cell potential settings).\
•	Fitting parameters: in the inputDict, the mass transfer coefficients (C1 and C2) can be defined, which can be obtained from literature or experiments.
•	Flow field type: in the scripts, the flow field design can be changed where Flow_field == 0 simulates a flow-through flow field, and Flow_field == 1 an interdigitated flow field.\
•	Fitting parameters: in the main scripts, the mass transfer coefficient, membrane resistance, and conductivity factors can be altered.\
•	Pressure drop physics: in the main scripts, you can select wether you want to use OpenPNM physics or physics based on the mechanical energy balance.\
•	Velocity dependent mass transfer coefficient: in the main manuscript you can select wether you want to use a velocity dependent or independent mass transfer coefficient.\
•	File name and folder: the file name and folder can be changed in the main scripts.\
•	Electrode type: the electrode you want to investigate with the pore network model can be loaded in the main scripts using the "input" folder.\
•	Additional operating conditions and electrolyte parameters: in the main script, the exchange current density and velocity should be changed accordingly.\
•	Network in series: in the main scripts you can change if you want to simulate a small electrode (the size of the loaded electrode) or if you want to use the network-in-series approach to simulate the size of longer electrode (e.g., as used in the experiments).\
•	Fitting patameters: for the Pressure_fitting_script and Flow_field_main scripts additional fitting parameters are required as input.

More information regarding these parameters can be found in the publications mentioned under the ‘Cite’ section. **Care must be taken when changing the electrolyte type, as different inputDicts must be used.**\
**If other extracted networks are used, make sure you double check the xyz coordinates, as they affect the boundary conditions and thus the simulations.**

Before running the code, the following folders need to be created:\
•	**input:** this folder contains the networks that will be simulated for the following scripts: Fitting_script_OpenPNM, PNM_main_V3_Final_OpenPNM, Flow_field_main, Local_km_fitting_script_initialize_fluid_field, and Pressure_fitting_script. This folder should contain the following files: the extracted network in .pnm format, .xlsx files with the experimental data necessary for the fitting of the pressure drop, mass transfer coefficient, and polarization curve.\
•	**Input_networks:** this folder contains the adjusted networks with scaled throat diameters (which is the output of the Flow_field_main script (.pnm files with scaled throat diameters)) that will be simulated for the following scripts: Fitting_script, PNM_main_V3_Final, and Local_km_fitting_script_load_in_fluid_field.\
•	**Input_Fluid_field:** for loading in the fluid field (which is the output of the Flow_field_main script (.pnm files)) for the scripts: Fitting_script, PNM_main_V3_Final, and Local_km_fitting_script_load_in_fluid_field.\
•	**output:** for saving the output data from the scripts.

**NOTE:** Before starting, carefully check wether you have all the information to run the scripts (input data, networks, fitting parameters required etc.).

## Contribution
PNM-RFB-electrode is developed using a wide range of open-source tools, including OpenPNM. The code has initially been developed as part of a PhD research project, but further contributions are very welcome.  

## Publications and citations
The code has been developed and described in the following two publications and PhD thesis. Please cite them when referring to the algorithm: 

The latest additions to the code (more accurate pressure drop and mass transfer coefficient descriptions) are discussed in Chapter 9 of the following PhD dissertation:
```bash
@phdthesis{heijden2023engineering,
      title = {Engineering porous electrodes for redox flow batteries - modeling, diagnostics, and manufacturing},
      author = {Maxime van der Heijden},
      isbn = {978-90-386-5878-0},
      year = {2023},
      publisher = {Eindhoven University of Technology},
}
```

The frame of the code is discussed and used in the following publication (using OpenPNM v2.2.0 and without new features):
```bash
@article{heijden2022assessing,
      title = {Assessing the Versatility and Robustness of Pore Network Modeling to Simulate Redox Flow Battery Electrode Performance},
      author = {Maxime van der Heijden and Rik van Gorp and Mohammad Amin Sadeghi and Jeffrey Gostick and Antoni Forner-Cuenca},
      journal = {Journal of the Electrochemical Society},
      volume = {169},
      number = {4},
      pages = {040505},
      year = {2022},
      doi = {10.1149/1945-7111/ac5e46},
}
```

The addition of the interdigitated flow field is discussed and used in the following publication (using OpenPNM v2.6.0):
```bash
@article{munoz2023versatile,
      title = {Understanding the role of electrode thickness on redox flow cell performance},
      author = {Vanesa Muñoz-Perales and Maxime van der Heijden and Victor de Haas and Jacky Olinga and Marcos Vera and Antoni Forner-Cuenca},
      journal = {ChemElectroChem},
      year = {2023},
      doi = {10.1002/celc.202300380},
}
```
