"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-18
License:  GNU GPL 3.0
"""

# QBlade settings
QBLADE_PATH = "bin\\QBlade\\QBladeCE_2.0.8.5.dll"                   # QBlade.dll path
QBLADE_CL_DEVICE = 0                                                # OpenCL device
QBLADE_GROUP_SIZE = 32                                              # OpenCL work-group size

# XFOIL settings
XFOIL_PATH = "bin\\XFoil\\xfoil.exe"                                # XFOIL.exe path
MAX_ITER = 100                                                      # Iteration limit, [-]

# Simulation settings
TURBINE_PATH = "data\\turbines\\DTU_10MW_v1.3\\DTU_10MW_RWT.trb"    # Turbine path
SIMULATION_PATH = "data\\simulations\\DTU_10MW_RWT.sim"             # Simulation path
WORKING_DIR = "temp"                                                # Working directory
N_TIMESTEP = 500                                                    # Number of timesteps, [-]
TIMESTEP = -1                                                       # Timestep index
BLADE = 1                                                           # Blade index

# Atmospheric settings
C_0 = 340.0                                                         # Speed of sound, [m/s]
RHO_0 = 1.225                                                       # Air density, [kg/m^3]
Z_0 = 0.01                                                          # Surface roughness length, [m]

# Acoustic settings
P_REF = 2E-5                                                        # Reference pressure, [Pa]
F_MIN = 20                                                          # Minimum frequency, [Hz]
F_MAX = 20000                                                       # Maximum frequency, [Hz]
F_REF = 1000                                                        # Reference frequency, [Hz]
BASE_10 = True                                                      # Use base 10?

# Aeroacoustic settings
RADIAL_CUTOFF = 0.4                                                 # Smallest radial position section, [-]
SPL_CORRECTION = True                                               # Use inflow noise 10dB SPL correction?
PROBE_TOP = 0.975                                                   # Probe location at the top surface, [-]
PROBE_BOT = 0.950                                                   # Probe location at the bottom surface, [-]
