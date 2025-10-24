"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-10-06
License:  GNU GPL 3.0
"""

# ISA model constants
R = 287.1               # Specific gas constant, [J/(k*Kg)]
GAMMA = 1.4             # Specific heat capacity ratio, [-]

# Sutherland's model constants
S = 110.4               # Sutherland's empirical coefficient, [K]
BETA_S = 1.458E-6       # Sutherland's empirical coefficient, [kg/(m*s*K^0.5)]

# Acoustic constants
P_REF_AIR = 2.0E-5      # Reference acoustic pressure in air, [Pa]
P_REF_WATER = 1.0E-6    # Reference acoustic pressure in water, [Pa]

# Hersbach's model constants
ALPHA_CH = 0.018        # Charnock's empirical coefficient, [-]
ALPHA_M = 0.11          # Hersbach's empirical coefficient, [-]
P = -12.0               # Hersbach's blending factor, [-]

# Universal constants
G = 9.80665             # Gravitational acceleration, [m/s^2]
KAPPA = 0.41            # von Karman constant, [-]

# Aeroacoustic settings
ALPHA_C = 1.0 / 0.7     # Freestream-Convection velocity ratio, [-]
B_C = 1.5               # Spanwise correlation coefficient, [-]

# Atmospheric absorption constants
P_REF = 101325.0        # Reference pressure, [Pa]
T_REF = 293.15          # Reference temperature, [K]
T_01 = 273.16           # Triple-point temperature of water, [K]
X_N = 0.7808            # Fractional molar concentration of nitrogen, [-]
X_O = 0.2095            # Fractional molar concentration of oxygen, [-]
THETA_N = 3352.0        # Characteristic vibrational temperature of nitrogen, [K]
THETA_O = 2239.1        # Characteristic vibrational temperature of oxygen, [K]
