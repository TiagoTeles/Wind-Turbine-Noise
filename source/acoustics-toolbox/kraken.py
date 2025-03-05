"""
TODO
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from readshd import readshd

filepath = "data/acoustics-toolbox/uniform"
filename = "uniform"

rmax = 20.0
Dmax = 1500.0
freq = 25.0

# Run KRAKEN
print("Running KRAKEN...")
os.system(f"./bin/krakenc.exe {filepath}/{filename}")
os.system(f"mv {filepath}/{filename}.prt {filepath}/kraken.prt")
os.system(f"./bin/field.exe {filepath}/{filename}")
os.system(f"mv field.prt {filepath}/field.prt")

# Read results
print("Reading output data...")

pressure, geometry = readshd(f"{filepath}/{filename}.shd", np.nan, np.nan, freq)

zs     = geometry["zs"]
rarray = geometry["rarray"]
zarray = geometry["zarray"]

p = np.squeeze(pressure, axis=(0, 1))
tl = 20*np.log10(abs(p))

plt.imshow(tl,extent=[0,rmax,-Dmax,0], aspect='auto',cmap='jet',vmin=-80,vmax=-30)

plt.colorbar()
plt.xlabel('Range (km)')
plt.ylabel('Depth (m)')
plt.xlim(0,rmax)
plt.ylim(-Dmax,0)
plt.clim(-95, -75)

plt.show()