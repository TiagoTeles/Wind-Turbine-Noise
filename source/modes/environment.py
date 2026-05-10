""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-30
License:  GNU GPL 3.0

Write data to .env files.

Classes:
    None

Functions:
    write_env

Exceptions:
    None
"""

def write_env(path, title, frequency, n_media, interpolation, bc_top, attenuation, ssp, bc_bot, \
              c_low, c_high, r_r_max, s_z, r_z):
    """
    Write the .env file.

    Arguments:
        path : str -- path to the .env file
        title : str -- title of the environment file
        frequency : float -- frequency, [Hz]
        interpolation : str -- SSP interpolation type
        bc_top : dict -- top boundary conditions
        ssp : dict -- sound speed profile
        attenuation : str -- attenuation units
        bc_bot : dict -- bottom boundary conditions
        c_low : float -- lower phase speed limit, [m/s]
        c_high : float -- upper phase speed limit, [m/s]
        r_r_max : float -- maximum receiver range, [m]
        s_z : np.ndarray -- source z coordinates, [m]
        r_z : np.ndarray -- receiver z coordinates, [m]
    """

    # Convert densities from [kg/m^3] to [g/cm^3]
    if bc_top["type"] == "A":
        bc_top["rho"] /= 1.0E3

    for i in range(len(ssp["profile"])):
        ssp["profile"][i]["rho"] /= 1.0E3

    if bc_bot["type"] == "A":
        bc_bot["rho"] /= 1.0E3

    # Convert positions from [m] to [km]
    r_r_max /= 1E3

    # Create the environment file
    f = open(path, "w", encoding="utf-8")

    # Write the title
    f.write(f"{title}\n")

    # Write the frequency
    f.write(f"{frequency}\n")

    # Write the number of media
    f.write(f"{n_media}\n")

    # Write the top options
    f.write(f"{interpolation}{bc_top["type"]}{attenuation}\n")

    if bc_top["type"] == "A":
        f.write(f"{bc_top["z"]} {bc_top["c_p"]} {bc_top["c_s"]} {bc_top["rho"]} \
                  {bc_top["alpha_p"]} {bc_top["alpha_s"]}\n")

    # Write the sound speed profile
    f.write(f"{ssp["n_mesh"]} {ssp["sigma"]} {ssp["depth"]}\n")

    for j in range(len(ssp[i]["profile"])):
        z = ssp[i]["profile"][j]["z"]
        c_p = ssp[i]["profile"][j]["c_p"]
        c_s = ssp[i]["profile"][j]["c_s"]
        rho = ssp[i]["profile"][j]["rho"]
        alpha_p = ssp[i]["profile"][j]["alpha_p"]
        alpha_s = ssp[i]["profile"][j]["alpha_s"]

        f.write(f"{z} {c_p} {c_s} {rho} {alpha_p} {alpha_s}\n")

    # Write the bottom properties
    f.write(f"{bc_bot["type"]} {bc_bot["sigma"]}\n")

    if bc_bot["type"] == "A":
        f.write(f"{bc_bot["z"]} {bc_bot["c_p"]} {bc_bot["c_s"]} {bc_bot["rho"]} {bc_bot["alpha_p"]} {bc_bot["alpha_s"]}\n")

    # Write the phase speed limits
    f.write(f"{c_low} {c_high}\n")

    # Write the maximum range
    f.write(f"{r_r_max}\n")

    # Write the source z coordinates
    f.write(f"{len(s_z)}\n")
    f.write(f"{s_z[0]} {s_z[-1]} /\n")

    # Write the receiver z coordinates
    f.write(f"{len(r_z)}\n")
    f.write(f"{r_z[0]} {r_z[-1]} /\n")

    # Close the file
    f.close()
