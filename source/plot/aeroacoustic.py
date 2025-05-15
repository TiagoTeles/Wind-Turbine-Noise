""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-12
License:  GNU GPL 3.0

Plot the results of the aeroacoustics calculations.

Classes:
    None

Functions:
    plot_spectrum
    plot_directivity

Exceptions:
    None
"""

import matplotlib.pyplot as plt
import numpy as np


def plot_spectrum(f, spl, f_min, f_max, name):
    """
    Plot the SPL spectrum.

    Parameters:
        f : np.array -- frequency, [Hz]
        spl : np.array -- SPL, [dB]
        f_min : float -- minimum frequency, [Hz]
        f_max : float -- maximum frequency, [Hz]
        name : str -- filename

    Returns:
        None
    """

    # Plot the spectrum
    plt.plot(f, spl)

    # Set the axis labels
    plt.xlabel("Frequency, [Hz]")
    plt.ylabel("SPL, [dB]")

    # Set the axis scales
    plt.xscale("log")
    plt.yscale("linear")

    # Set the axis limits
    plt.xlim(f_min, f_max)
    plt.ylim(0, 10 * np.ceil(np.max(spl) / 10))

    # Set the grid lines
    plt.grid(which="both")

    # Save the plot
    plt.savefig(name + ".pdf")

    # Show the plot
    plt.show()


def plot_directivity(yaw, spl, name):
    """
    Plot the directivity pattern.

    Parameters:
        yaw : np.array -- yaw angle, [rad]
        spl : np.array -- SPL, [dB]
        name : str -- filename

    Returns:
        None
    """

    # Plot the directivity pattern
    plt.polar(yaw, spl)

    # Set the axis limits
    plt.xlim(0, 2*np.pi)
    plt.ylim(0, 10 * np.ceil(np.max(spl) / 10))

    # Save the plot
    plt.savefig(name + ".pdf")

    # Show the plot
    plt.show()


def plot_contour(azimuth, r, spl, name):
    """
    Plot the SPL map.

    Parameters:
        azimuth : np.array -- azimuth angle, [rad]
        r : np.array -- r/R, [m]
        spl : np.array -- SPL, [dB]
        name : str -- filename

    Returns:
        None
    """

    # Create the meshgrid
    R, AZIMUTH = np.meshgrid(r, azimuth)

    # Plot the SPL map
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    cs = ax.contourf(AZIMUTH, R, spl, cmap="rainbow")

    # Set the colorbar
    fig.colorbar(cs, label="SPL, [dB]")

    # Save the plot
    plt.savefig(name + ".pdf")

    # Show plot
    plt.show()
