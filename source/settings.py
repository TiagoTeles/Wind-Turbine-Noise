"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-11
License:  GNU GPL 3.0
"""

import numpy as np


# Atmosphere settings
P_0 = 101325.0  # Pressure, [Pa]
T_0 = 288.15    # Temperature, [K]
H_R = 0.80      # Relative humidity, [-]

# Wind turbine settings
LAT = np.radians(41.6865)   # Latitude, [rad]
LON = np.radians(-9.0574)   # Longitude, [rad]

# Ocean settings
BATHYMETRY_PATH = "data\\study_domain\\windfloat_atlantic\\bathymetry.csv"      # Bathymetry path
TEMPERATURE_PATH = "data\\study_domain\\windfloat_atlantic\\temperature.csv"    # Temperature path
SALINITY_PATH = "data\\study_domain\\windfloat_atlantic\\salinity.csv"          # Salinity path
SEABED_PATH = "data\\study_domain\\windfloat_atlantic\\seabed.csv"              # Seabed path

# QBlade settings
QBLADE_DLL_PATH = "bin\\QBlade\\QBladeCE_2.0.9.4.dll"   # QBlade.dll path
QBLADE_OPENCL_DEVICE = 0                                # OpenCL device
QBLADE_OPENCL_GROUP_SIZE = 32                           # OpenCL work-group size

# QBlade simulation settings
QBLADE_SIMULATION_PATH = "data\\wind_turbine\\IEA_22MW_RWT\\IEA-22-280-RWT-Monopile.sim"    # Simulation path
QBLADE_RESULTS_PATH = "results\\wind_turbine\\IEA_22MW_RWT\\qblade.txt"                     # Results path

# One-Third octave settings
BASE_10 = True      # Use base-10 formulation?
F_MIN = 20.0        # Minimum frequency, [Hz]
F_MAX = 20000.0     # Maximum frequency, [Hz]

# Aeroacoustic discretisation settings
N_AZIMUTH = 12      # Number of azimuthal positions, [-]
ASPECT_RATIO = 2.0  # Panel aspect ratio, [-]
RADIUS_CUTOFF = 0.4 # Radial cutoff, [-]

# XFOIL settings
XFOIL_PATH = "bin\\XFOIL\\xfoil.exe"    # xfoil.exe path
XFOIL_MAX_ITERATIONS = 100              # Iteration limit, [-]
XFOIL_TRANSITION_UPPER = 0.065          # Upper transition position, [-]
XFOIL_TRANSITION_LOWER = 0.200          # Lower transition position, [-]
XFOIL_AMPLIFICATION = 9.0               # Critical amplification factor, [-]

# Boundary layer settings
PROBE_UPPER = 0.975 # Upper probe position, [-]
PROBE_LOWER = 0.950 # Lower probe position, [-]

# # NM80 settings
# P_0 = 103236.2  # Pressure, [Pa]
# T_0 = 292.15    # Temperature, [K]
# H_R = 0.70      # Relative humidity, [-]

# QBLADE_SIMULATION_PATH = "data\\wind_turbine\\NM_2MW_RWT\\NM-2-80-RWT.sim"                  # Simulation path
# QBLADE_RESULTS_PATH = "results\\wind_turbine\\NM_2MW_RWT\\qblade.txt"                       # Results path
