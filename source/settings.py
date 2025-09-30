"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-23
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
BATHYMETRY_PATH = "data\\environments\\windfloat_atlantic\\bathymetry.csv"      # Bathymetry path
TEMPERATURE_PATH = "data\\environments\\windfloat_atlantic\\temperature.csv"    # Temperature path
SALINITY_PATH = "data\\environments\\windfloat_atlantic\\salinity.csv"          # Salinity path
SEABED_PATH = "data\\environments\\windfloat_atlantic\\seabed.csv"              # Seabed path

# QBlade settings
QBLADE_PATH = "bin\\QBlade\\QBladeCE_2.0.9.3.dll"   # QBlade.dll path
OPENCL_DEVICE = 0                                   # OpenCL device
OPENCL_GROUP_SIZE = 32                              # OpenCL work-group size

# QBlade simulation settings
SIMULATION_PATH = "data\\turbines\\IEA_22MW_RWT\\IEA-22-280-RWT-Monopile.sim"   # Simulation path
RESULTS_PATH = "results\\IEA_22MW_RWT\\windfloat_atlantic\\qblade.txt"          # Results path

# Acoustic settings
F_MIN = 20.0        # Minimum frequency, [Hz]
F_MAX = 20000.0     # Maximum frequency, [Hz]
F_REF = 1000.0      # Reference frequency, [Hz]
BASE_10 = True      # Use base-10 formulation?

# Aeroacoustic settings
# N_AZIMUTH = 12      # Number of azimuthal positions, [-]
# BLADE_ID = 0        # QBlade blade index
AR = 2.0            # Panel aspect ratio, [-]
