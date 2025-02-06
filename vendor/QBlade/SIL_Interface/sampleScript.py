import os
import sys
from ctypes import *
from QBladeLibrary import QBladeLibrary

# Define the directory where the QBlade library is located
dll_directory = "../"

# On Windows systems, we update the PATH environment variable to include the QBlade directory. 
# This ensures that the required SSL libraries (e.g., libssl and libcrypto) are properly located and loaded.
# If experiencing issues with this DLL in a Windows Python environment see:
# https://docs.qblade.org/src/license/license_files.html#resolving-openssl-issues-on-windows
if os.name == 'nt':  # 'nt' indicates Windows
    os.environ["PATH"] = os.path.abspath(dll_directory) + ";" + os.environ.get("PATH", "")

# Search the directory below for library files matching the pattern QBlade*.dll or QBlade*.so
dll_files = [f for f in os.listdir(dll_directory) if 'QBlade' in f and ('.dll' in f or '.so' in f)]

# Check if any matching files are found
if not dll_files:
    print('No matching QBlade*.dll or QBlade*.so files found in the specified directory:',os.path.abspath(dll_directory))
    sys.exit(1)  # Exit the script with a non-zero status to indicate an error

# Use the first matching file
dll_file_path = os.path.join(dll_directory, dll_files[0])

# Display the selected shared library file
print(f'Using shared library file: {dll_file_path}')

# Create an object of the class 'QBladeLibrary' that contains the API
QBLADE = QBladeLibrary(dll_file_path)    

# Creation of a QBlade instance from the library
QBLADE.createInstance(1,32)

# Loading a project or sim-file, in this case the DTU_10MW_Demo project or simulation definition file
#QBLADE.loadSimDefinition(b"./DTU_10MW_Demo.sim") #uncomment this line to load a simulation definition file
QBLADE.loadProject(b"./NREL_5MW_Sample.qpr") 

# Initializing the sim and ramp-up phase, call before starting the simulation loop
QBLADE.initializeSimulation()

# We will run the simulation for 500 steps before storing the results
number_of_timesteps = 500

# Start of the simulation loop
for i in range(number_of_timesteps):

    #advance the simulation
    success = QBLADE.advanceTurbineSimulation() 
    
    # Check if the simulation step was successful
    if not success:  # If success is False, exit the loop
        print(f"Simulation failed at timestep {i}. Exiting loop.")
        break
    
    # Assign the c-type double array 'loads' with length [6], initialized with zeros
    loads = (c_double * 6)(0,0,0,0,0,0) 
    # Retrieve the tower loads and store the in the array 'loads' by calling the function getTowerBottomLoads_at_num()
    QBLADE.getTowerBottomLoads_at_num(loads,0)
    
    # Uncomment the next line to try changing the position of the turbine dynamically
    #QBLADE.setTurbinePosition_at_num(-0.2*i,0,0,0,i*0.1,i*0.1,0) 
    
    # Example how to extract a variable by name from the simulation, call as often as needed with different variable names, extracting rpm and time in the lines below
    rpm = QBLADE.getCustomData_at_num(b"Rotational Speed [rpm]",0,0) 
    time = QBLADE.getCustomData_at_num(b"Time [s]",0,0) #example how to extract the variable 'Time' by name from the simulation
    AoA = QBLADE.getCustomData_at_num(b"Angle of Attack at 0.25c (at section) Blade 1 [deg]",0.85,0) #example how to extract the variable 'Angle of Attack' by name at 85% blade length from the simulation 
    
    # Example how to extract a 3 length double array with the x,y,z windspeed components at a global position of x=-50,Y=0,Z=100m from the simulation
    windspeed = (c_double * 3)(0,0,0) 
    QBLADE.getWindspeed(-50,0,100,windspeed)
    
    # Assign the c-type double array 'ctr_vars' with length [5], initialized with zeros
    ctr_vars = (c_double * 5)(0); 
    # Advance the turbine controller and store the controller signals in the array 'ctr_vars'
    QBLADE.advanceController_at_num(ctr_vars,0)
    
    # Pass the controller signals in 'ctr_vars' to the turbine by calling setControlVars_at_num(ctr_vars,0) 
    QBLADE.setControlVars_at_num(ctr_vars,0) 
    
    # Print out a few of the recorded data, in this case torque, tower bottom force along z (weight force) and rpm
    print("Time:","{:3.2f}".format(time),"   Windspeed:","{:2.2f}".format(windspeed[0]),"  Torque:","{:1.4e}".format(ctr_vars[0]),"    RPM:","{:2.2f}".format(rpm),"   Pitch:","{:2.2f}".format(ctr_vars[2]),"   AoA at 85%:","{:2.2f}".format(AoA))

# The simulation loop ends here after all 'number_of_timesteps have been evaluated
	
# Storing the finished simulation in a project as NREL_5MW_Sample_completed, you can open this file to view the results of the simulation inside QBlade's GUI
QBLADE.storeProject(b"./NREL_5MW_Sample_completed.qpr")

# Storing the simulation results in QBlade ASCII format in the file NREL_5MW_Sample_results.txt
QBLADE.exportResults(0,b"./",b"NREL_5MW_Sample_results",b"")

# Unloading the qblade library
QBLADE.unload() 