"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-10-21
License:  GNU GPL 3.0
"""

import numpy as np


# Atmosphere properties
P_0 = 101325.0  # Pressure, [Pa]
T_0 = 288.15    # Temperature, [K]
H_R = 0.80      # Relative humidity, [-]

# Study domain properties
LATITUDE = np.radians(41.6865)  # Study domain latitude, [rad]
LONGITUDE = np.radians(-9.0574) # Study domain longitude, [rad]

# Ocean properties
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
# QBLADE_SIMULATION_PATH = "data\\wind_turbine\\NM_2MW_RWT\\NM-2-80-RWT.sim"                  # Simulation path
# QBLADE_RESULTS_PATH = "results\\wind_turbine\\NM_2MW_RWT\\qblade.txt"                       # Results path

# Acoustic settings
F_MIN = 20.0        # Minimum frequency, [Hz]
F_MAX = 20000.0     # Maximum frequency, [Hz]
F_REF = 1000.0      # Reference frequency, [Hz]
BASE_10 = True      # Use base-10 formulation?

# Aeroacoustic settings
N_AZIMUTH = 12  # Number of azimuthal positions, [-]
AR = 2.0        # Panel aspect ratio, [-]
