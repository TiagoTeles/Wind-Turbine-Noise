"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-18
License:  GNU GPL 3.0

Settings for the various scripts.

Classes:
    None

Functions:
    None

Exceptions:
    None
"""

# QBlade configuration
QBLADE_PATH = "bin\\QBlade\\QBladeCE_2.0.8.5.dll"                   # QBlade.DLL path
QBLADE_CL_DEVICE = 1                                                # OpenCL device
QBLADE_GROUP_SIZE = 32                                              # OpenCL work-group size

# XFOIL configuration
XFOIL_PATH = "bin\\XFoil\\xfoil.exe"                                # XFOIL.exe path

# Simulation specification
TURBINE_PATH = "data\\turbines\\DTU_10MW_v1.3\\DTU_10MW_RWT.trb"    # Turbine definition path
SIMULATION_PATH = "data\\simulations\\DTU_10MW_RWT.sim"             # Simulation definition path
WORKING_DIR = "temp"                                                # Working directory
N_TIMESTEP = 500                                                    # Number of timesteps, [-]
BLADE = 1                                                           # Blade ID TODO: Does not work for unsteady case!
TIMESTEP = -1                                                       # Timestep index to use

# Aerodynamic settings
C_0 = 340.0                                                         # Speed of sound, [m/s] TODO: Base on temperature?
RHO_0 = 1.225                                                       # Air density, [kg/m^3] TODO: Use QBalde value?
Z_0 = 0.01                                                          # Surface roughness length, [m] TODO: Base on what?

# Acoustic settings
BASE_10 = True                                                      # Use base 10?
F_MIN = 20                                                          # Minimum frequency, [Hz]
F_MAX = 20000                                                       # Maximum frequency, [Hz]
F_REF = 1000                                                        # Reference frequency, [Hz]
P_REF = 2E-5                                                        # Reference pressure, [Pa]

# Aeroacosutic settings
SPL_CORRECTION = True                                               # Apply inflow noise 10dB SPL correction?
PROBE_TOP = 0.975                                                   # Probe location at the top surface, [-]
PROBE_BOT = 0.950                                                   # Probe location at the bottom surface, [-]
THICKNESS_CUTOFF = 0.3                                              # Thickness cutoff for the airfoil, [-]