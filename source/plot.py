""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-18
License:  GNU GPL 3.0

Plot the results of the aeroacoustic calculations.

Classes:
    None

Functions:
    spectrum
    directivity
    map

Exceptions:
    None
"""

import matplotlib.pyplot as plt
import numpy as np

from settings import *


def spectrum(f, spl, n_blades, x):
    """
    Plot the SPL spectrum at the rotor axis.

    Parameters:
        f : np.array -- frequency array
        spl : np.array -- SPL array
        n_blades : int -- number of blades, [-]
        x : float -- observer distance, [m]

    Returns:
        None
    """

    # Determine the SPL spectrum of each panel
    spl_panel = 10 * np.log10(np.mean(np.pow(10, spl/10), axis=2))

    # Determine the SPL spectrum of each blade
    spl_blade = 10 * np.log10(np.sum(np.pow(10, spl_panel/10), axis=1))

    # Determine the SPL spectrum of the turbine
    spl_turbine = 10 * np.log10(n_blades * np.pow(10, spl_blade/10))

    # Determine the total SPL of the turbine
    spl_total = 10 * np.log10(np.sum(np.pow(10, spl_blade/10), axis=0))

    # Plot the panel spectra
    for i in range(spl.shape[1]):
        plt.plot(f, spl[:, i], label=f"Panel {i}")

    # Plot the blade spectra
    plt.plot(f, spl_blade, label="Blade", lw=2, ls='--')

    # Plot the turbine spectra
    plt.plot(f, spl_turbine, label="Turbine", lw=2, ls='-.')

    # Set the plot title
    plt.title(f"SPL spectra of the turbine at a downwind distance of {x:.0f} [m] from the rotor axis. Total SPL: {spl_total:.0f} dB")

    # Set axis labels
    plt.xlabel("Frequency, [Hz]")
    plt.ylabel("SPL, [dB]")

    # Set axis scales
    plt.xscale("log")
    plt.yscale("linear")

    # Set axis limits
    plt.xlim(F_MIN, F_MAX)
    plt.ylim(0, spl_total)

    # Show plot
    plt.grid(which="both")
    plt.legend()
    plt.show()


def directivity(spl, n_blades, n_azimuthal, n_angles, h, r):
    """
    Plot the directivity pattern about the z-axis.

    Parameters:
        spl : np.array -- SPL array
        n_blades : int -- number of blades, [-]
        n_azimuthal : int -- number of azimuthal angles, [-]
        n_angles : int -- number of observer angles, [-]
        h : float -- observer height, [m]
        r : float -- observer distance, [m]

    Returns:
        None
    """

    # Determine the SPL of each panel
    spl_panel = 10 * np.log10(np.sum(np.pow(10, spl/10), axis=0))

    # Determine the SPL of each blade
    spl_blade = 10 * np.log10(np.sum(np.pow(10, spl_panel/10), axis=0))

    # Seperate the azimuthal and observer angle axis
    spl_blade = spl_blade.reshape(n_angles, n_azimuthal)

    # Determine the SPL of the turbine
    p2_blade = np.pow(10, spl_blade/10) * np.square(P_REF)

    p2_turbine = np.zeros(p2_blade.shape)

    for i in range(n_blades):

        # Determine the azimuthal index
        index = i * n_azimuthal // n_blades
        
        # Add the contribution of each blade
        p2_turbine += np.roll(p2_blade, index, axis=1)
    
    spl_turbine = 10 * np.log10(p2_turbine / np.square(P_REF))

    # Average the SPL over all azimuthal angles
    spl_total = 10 * np.log10(np.mean(np.pow(10, spl_turbine/10), axis=1))

    # Plot the directivity pattern
    plt.polar(np.linspace(0, 2*np.pi, n_angles), spl_total)

    # Set the plot title
    plt.title(f"Directivity pattern of the turbine at an observer distance of {r:.0f} [m] and height of {h:.2f} [m].")

    # Set axis labels
    # plt.xlabel(r"$\gamma$, [deg]")
    # plt.ylabel("SPL, [dB]")

    # Set axis scales
    plt.xscale("linear")
    plt.yscale("linear")

    # Set axis limits
    plt.xlim(0, 2*np.pi)
    plt.ylim(40, 80)

    # Show plot
    plt.show()


def map(spl, n_panels, n_azimuthal, x):
    """
    Plot the SPL of each panel at each azimuthal angle.

    Parameters:
        spl : np.array -- SPL array
        n_panels : int -- number of panels, [-]
        n_azimuthal : int -- number of azimuthal angles
        x : float -- observer distance, [m]
    
    Returns:
        None
    """

    # Determine the SPL of each panel
    spl_panel = 10 * np.log10(np.sum(np.pow(10, spl/10), axis=0))

    # Determine R and THETA for datapoint
    r = np.linspace(0, 1, n_panels)
    gamma = np.linspace(0, 2*np.pi, n_azimuthal)

    R, THETA = np.meshgrid(r, gamma)

    # Plot the heatmap
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    cs = ax.contourf(THETA, R, spl_panel.T, levels=np.linspace(0, np.max(spl_panel), 21), cmap="rainbow")

    # Set the colorbar
    fig.colorbar(cs, label="SPL, [dB]")

    # Set the plot title
    plt.title(f"SPL of one blade along an entire CCW rotation \
              at a downwind distance of {x:0f} [m].")

    # Set axis labels
    # plt.xlabel(r"$\gamma$, [deg]")
    # plt.ylabel("r/R, [-]")

    # Set axis scales
    plt.xscale("linear")
    plt.yscale("linear")

    # Set axis limits
    plt.xlim(0, 2*np.pi)
    plt.ylim(0, 1)

    # Show plot
    plt.show()
