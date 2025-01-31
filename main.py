import matplotlib.pyplot as plt
import numpy as np
import os

filepath = "data/acoustics-toolbox/uniform/uniform"

# Run KRAKEN
print( "Running KRAKEN..." )
os.system(f"./bin/krakenc.exe {filepath}")
os.system(f"./bin/field.exe {filepath} > field.prt")

# print( "Reading output data..." )

# filename = 'elba.shd'
# xs = nan
# ys = nan
# pressure,geometry = readshd(filename,xs,ys)

# Modes,bc = readmod('elba.mod')

# zs     = geometry["zs"]
# rarray = geometry["rarray"]
# zarray = geometry["zarray"]
# print(pressure)
# p = squeeze( pressure, axis=(0,1) )
# tl = 20*log10( abs( p ) )

# figure(1)
# imshow(tl,extent=[0,rmax,-Dmax,0], aspect='auto',cmap='jet',vmin=-80,vmax=-30)
# colorbar()
# plot([0,rmax],[-98  , -98],'k',linewidth=2)
# plot([0,rmax],[-103, -103],'k',linewidth=2)
# plot([0,rmax],[-128, -128],'k',linewidth=2)
# plot(rs,-zs,marker="<",markersize=16,color="k")
# xlabel('Range (m)')
# ylabel('Depth (m)')
# title('KRAKEN - Elba waveguide')
# ylim(-Dmax,0)

# phi = Modes["phi"]
# z   = Modes["z"]
# k   = Modes["k"]

# figure(2)
# for i in range(4):  
#    rphi = real( phi[ : , i ] ) 
#    iphi = imag( phi[ : , i ] )
#    thetitle = 'Z_'  + str(i+1) + '(z)' 
#    subplot(1,4,i+1)
#    plot(rphi,-z,iphi,-z,'r--')
#    title( thetitle )
#    grid(True)
# subplot(141)
# ylabel('Depth (m)')

# show()

# print("done.")
