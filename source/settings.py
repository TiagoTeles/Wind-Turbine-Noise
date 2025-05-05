"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-05
License:  GNU GPL 3.0
"""

# QBlade settings
QBLADE_PATH = "bin\\QBlade\\QBladeCE_2.0.8.6.dll"           # QBlade.dll path
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
G = 9.81                                                    # Gravitational acceleration, [m/s^2]
KAPPA = 0.41                                                # Von Karman constant, [-]

# Environment properties
C_0 = 340.3                                                 # Speed of sound, [m/s]
RHO_0 = 1.225                                               # Air density, [kg/m^3]
NU_0 = 1.46E-5                                              # Kinematic viscosity, [m^2/s]

# Acoustic settings
P_REF_AIR = 2E-5                                            # Reference pressure, [Pa]
P_REF_WATER = 1E-6                                          # Reference pressure, [Pa]
F_MIN = 20                                                  # Minimum frequency, [Hz]
F_MAX = 20000                                               # Maximum frequency, [Hz]
F_REF = 1000                                                # Reference frequency, [Hz]
BASE_10 = True                                              # Use base 10?

# Aeroacoustic settings
SPL_CORRECTION = False                                      # Use inflow noise 10dB SPL correction?
RADIAL_CUTOFF = 0.4                                         # TBLTE noise radial cutoff, [-]
PROBE_TOP = 0.975                                           # Suction side probe location, [-]
PROBE_BOT = 0.950                                           # Pressure side probe location, [-]
SPEED_RATIO = 1.0/0.7                                       # U / U_c, [-]
CORRELATION_COEFFICIENT = 1.5                               # Correlation coefficient, [-]
