"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-21
License:  GNU GPL 3.0
"""

# QBlade settings
QBLADE_PATH = "bin\\QBlade\\QBladeCE_2.0.8.5.dll"           # QBlade.dll path
QBLADE_CL_DEVICE = 0                                        # OpenCL device
QBLADE_GROUP_SIZE = 32                                      # OpenCL work-group size

# XFOIL settings
XFOIL_PATH = "bin\\XFOIL\\xfoil.exe"                        # xfoil.exe path
MAX_ITER = 100                                              # Iteration limit, [-]

# Simulation settings
SIMULATION_PATH = "data\\DTU_10MW_RWT\\DTU_10MW_RWT.sim"    # Simulation path
N_TIMESTEP = 500                                            # Number of timesteps, [-]
TIMESTEP = -1                                               # Timestep index
BLADE = 1                                                   # Blade index

# Universal constants
G = 9.0665                                                  # Gravitational acceleration, [m/s^2]
KAPPA = 0.41                                                # Von Karman constant, [-]

# Environment properties
C_0 = 340.0                                                 # Speed of sound, [m/s]
RHO_0 = 1.225                                               # Air density, [kg/m^3]
NU_0 = 1.48E-5                                              # Kinematic viscosity, [m^2/s]

# Acoustic settings
P_REF = 2E-5                                                # Reference pressure, [Pa]
F_MIN = 20                                                  # Minimum frequency, [Hz]
F_MAX = 20000                                               # Maximum frequency, [Hz]
F_REF = 1000                                                # Reference frequency, [Hz]
BASE_10 = True                                              # Use base 10?

# Aeroacoustic settings
SPL_CORRECTION = True                                       # Use inflow noise 10dB SPL correction?
RADIAL_CUTOFF = 0.4                                         # TBLTE noise radial cutoff, [-]
PROBE_TOP = 0.975                                           # Suction side probe location, [-]
PROBE_BOT = 0.950                                           # Pressure side probe location, [-]
