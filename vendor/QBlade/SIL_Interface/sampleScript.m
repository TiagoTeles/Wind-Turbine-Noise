%%
% This is a minimum working example on how to interface with the QBlade
% library. A project called 'NREL_5MW_Sample.qpr' has been provided to work with. 
% This script establishes a connection to the QBlade library, calls the library
% QBLADE and then initializes QBlade and the project. It runs the
% simulation for 500 simulation steps, and outputs some values. The time length
% of the simulation equals the number of simulation steps times the time
% step with which the simulation has been set up. At the end, the QBlade
% instance is closed again.

%%
clear all
close all 
clc

% Search the directory below for library files matching the pattern QBlade*.dll or QBlade*.so
libSearchDirectory = '../';
sharedLibFiles = dir(fullfile(libSearchDirectory, '*QBlade*'));
sharedLibFiles = sharedLibFiles(contains({sharedLibFiles.name}, {'.dll', '.so'}));

if isempty(sharedLibFiles)
    fprintf('No matching QBlade*.dll files or QBlade*.so found in the specified directory.');
    return;
end
    
% Use the first matching dll file and print out its name
sharedLibFilePath = fullfile(libSearchDirectory, sharedLibFiles(1).name);
fprintf('Using DLL file: %s\n', sharedLibFilePath);

% Create an object of the class 'QBladeLibrary' that contains the API
QBLADE = QBladeLibrary(sharedLibFilePath);

QBLADE.createInstance(1,32);

% Since matlab is unable to display the console output from the library, we
% store the output in a log file
QBLADE.setLogFile(fullfile('.', 'LogFile.txt'))

QBLADE.loadProject('NREL_5MW_Sample.qpr')

QBLADE.initializeSimulation()

number_of_timesteps = 500; 

f = waitbar(0,'Initializing Simulation') ;

for i = 1:1:number_of_timesteps
    
    % Advance the simulation
    success = QBLADE.advanceTurbineSimulation();
    
    % Check if the simulation step was successful
    if ~success
        fprintf('Simulation failed at timestep %d. Exiting loop.\n', i);
        break; % Exit the loop
    end
    
    % Assign the c-type double array 'loads' with length [6], initialized with zeros
    loads = libpointer('doublePtr',zeros(6,1));
    % Retrieve the tower loads and store the in the array 'loads' by calling the function getTowerBottomLoads_at_num()
    QBLADE.getTowerBottomLoads_at_num(loads,0);
    % De-referencing the 'loads' pointer and accessing its first value
    loads.Value(1);
    
    % Uncomment the next line to try changing the position of the turbine dynamically
    %QBLADE.setTurbinePosition_at_num(-0.2*i,0,0,0,i*0.1,i*0.1,0)
    
    % Example how to extract a variable by name from the simulation, call as often as needed with different variable names, extracting rpm and time in the lines below
    rpm = QBLADE.getCustomData_at_num('Rotational Speed [rpm]',0,0);
    t = QBLADE.getCustomData_at_num('Time [s]',0,0);  %example how to extract the variable 'Time' by name from the simulation
    AoA = QBLADE.getCustomData_at_num('Angle of Attack at 0.25c (at section) Blade 1 [deg]',0.85,0); %example how to extract the variable 'Angle of Attack' by name at 85% blade length from the simulation 
    
    % Example how to extract a 3 length double array with the x,y,z windspeed components at a global position of x=-50,Y=0,Z=100m from the simulation
    windspeed = libpointer('doublePtr',zeros(3,1)); 
    QBLADE.getWindspeed(-50,0,100,windspeed);
    
    % Assign the c-type double array 'ctr_vars' with length [5], initialized with zeros
    ctr_vars = libpointer('doublePtr',zeros(5,1));
    % Advance the turbine controller and store the controller signals in the array 'ctr_vars'
    QBLADE.advanceController_at_num(ctr_vars,0)
    
    % Pass the controller signals in 'ctr_vars' to the turbine by calling setControlVars_at_num(ctr_vars,0) 
    QBLADE.setControlVars_at_num(ctr_vars,0)
    
    fprintf('Time: %3.2f	Windspeed: %2.2f    Torque: %1.4e	RPM: %2.2f	Pitch: %2.2f    AoA at 85%%: %2.2f\n',t,windspeed.Value(1),ctr_vars.Value(1),rpm,ctr_vars.Value(3),AoA);
    
    waitbar(i/number_of_timesteps,f,'QBlade Simulation Running')

end

close(f)

% Storing the finished simulation in a project as NREL_5MW_Sample_completed, you can open this file to view the results of the simulation inside QBlade's GUI
QBLADE.storeProject('./NREL_5MW_Sample_completed.qpr')

% Storing the simulation results in QBlade ASCII format in the file NREL_5MW_Sample_results.txt
QBLADE.exportResults(0,'./','NREL_5MW_Sample_Results','')

% Closing the instance of the shared library, if this fail it can lead to unexpected behavior
QBLADE.closeInstance()

% Unloading the shared library
QBLADE.unload()

