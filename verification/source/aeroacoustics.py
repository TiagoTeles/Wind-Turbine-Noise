"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2026-05-05
License:  GNU GPL 3.0
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from source.settings import AEROACOUSTIC_RESULTS_PATH, F_MAX, F_MIN
from verification.source.settings import CAA_3DS_PATH, CAA_DLR_PATH, CAA_DTU_PATH
from verification.source.settings import CAA_IAG_PATH, CAA_NREL_PATH, CAA_TNO_PATH


# Read the simulation results
data = np.load(AEROACOUSTIC_RESULTS_PATH)

# Read the verification data
caa_3ds = pd.read_csv(CAA_3DS_PATH)
caa_dlr = pd.read_csv(CAA_DLR_PATH)
caa_dtu = pd.read_csv(CAA_DTU_PATH)
caa_iag = pd.read_csv(CAA_IAG_PATH)
caa_nrel = pd.read_csv(CAA_NREL_PATH)
caa_tno = pd.read_csv(CAA_TNO_PATH)

# Extract the data
azimuth = data["azimuth"]
frequency = data["frequency"]
radius = data["radius"]
spl_le = data["spl_le"]
spl_te = data["spl_te"]
x_o = data["x_o"]

# Offset the inflow turbulence noise SPL by 10.0 [dB]
spl_le -= 10.0

# Determine the blade SPLs
spl_le_blade = 10.0 * np.log10(np.sum(np.pow(10.0, spl_le / 10.0), axis=2))
spl_te_blade = 10.0 * np.log10(np.sum(np.pow(10.0, spl_te / 10.0), axis=2))

# Deterine the azimuthally averaged SPLs
spl_le_turbine = 10.0 * np.log10(np.mean(np.pow(10.0, spl_le_blade / 10.0), axis=0))
spl_te_turbine = 10.0 * np.log10(np.mean(np.pow(10.0, spl_te_blade / 10.0), axis=0))

# Plot the inflow turbulence noise spectra
plt.semilogx(frequency, spl_le_turbine[:, 0], label="T. Teles")
plt.semilogx(caa_3ds["frequency"], caa_3ds["spl_le"], ls="--", label="3DS")
plt.semilogx(caa_dlr["frequency"], caa_dlr["spl_le"], ls="--", label="DLR")
plt.semilogx(caa_dtu["frequency"], caa_dtu["spl_le"], ls="--", label="DTU")
plt.semilogx(caa_iag["frequency"], caa_iag["spl_le"], ls="--", label="IAG")
plt.semilogx(caa_nrel["frequency"], caa_nrel["spl_le"], ls="--", label="NREL")
plt.semilogx(caa_tno["frequency"], caa_tno["spl_le"], ls="--", label="TNO")
plt.xlabel("Frequency, [Hz]")
plt.ylabel("SPL, [dB]")
plt.xlim(F_MIN, F_MAX)
plt.ylim(15.0, 55.0)
plt.grid(which="both")
plt.legend()
plt.show()

# Plot the TBLTE noise spectra
plt.semilogx(frequency, spl_te_turbine[:, 0], label="T. Teles")
plt.semilogx(caa_3ds["frequency"], caa_3ds["spl_te"], ls="--", label="3DS")
plt.semilogx(caa_dlr["frequency"], caa_dlr["spl_te"], ls="--", label="DLR")
plt.semilogx(caa_dtu["frequency"], caa_dtu["spl_te"], ls="--", label="DTU")
plt.semilogx(caa_iag["frequency"], caa_iag["spl_te"], ls="--", label="IAG")
plt.semilogx(caa_nrel["frequency"], caa_nrel["spl_te"], ls="--", label="NREL")
plt.semilogx(caa_tno["frequency"], caa_tno["spl_te"], ls="--", label="TNO")
plt.xlabel("Frequency, [Hz]")
plt.ylabel("SPL, [dB]")
plt.xlim(F_MIN, F_MAX)
plt.ylim(15.0, 55.0)
plt.grid(which="both")
plt.legend()
plt.show()
