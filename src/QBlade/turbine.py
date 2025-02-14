""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-14
License:  GNU GPL 3.0

Store turbine data.

Classes:
    Turbine

Functions:
    None

Exceptions:
    None
"""

import os
import sys
import numpy as np

from aero import Aero
from misc import parse


class Turbine:
    """ 
    A class to store the turbine data. 
    
    Methods:
        __init__ -- parse the turbine file
        
    Attributes:
        aerodynamic_panel_loads : bool -- enable distributed aero panel loads and gradients?
        am_gb : float -- Am constant for the GORMONT-BERG dynamic stall model
        bem_speedup : float -- initial BEM convergence acceleration time [s]
        bem_tip_loss : bool -- use BEM tip loss factor?
        blade_path : str -- path of the blade file
        bound_core_radius : float -- the fixed core radius of the bound blade vortex (fraction of local chord) [-]
        controller_path : str -- controller path (leave blank if unused)
        controller_type : int -- type of turbine controller (0 = none, 1 = BLADED, 2 = DTU, 3 = TUB)
        convection_type : int -- wake convection type (0 = BL, 1 = HH, 2 = LOC)
        conversion_length : float -- wake conversion length (to particles) [-]
        discretisation_type : int -- type of rotor discretization (0 = from Bladetable, 1 = linear, 2 = cosine)
        dynamic_stall_type : int -- dynamic stall model (0 = none; 1 = OYE; ; 2 = IAG; 3 = GORMONT-BERG; 4 = ATEFLAP)
        first_wake_row_length : float -- first wake row length [-]
        gamma_epsilon : float -- relative gamma convergence criteria [m^2/s]
        gamma_iterations : int -- maximum number of gamma iterations [-]
        gamma_relaxation : float -- relaxation factor used in the gamma iteration [-]
        geometric_stiffness : bool -- enable geometric stiffness?
        himmelskamp_stall_delay : bool -- include the Himmelskamp Stall delay? (HAWT only)
        iag_parameters : list -- the parameters for the IAG DS model (not implemented)
        interpolation_type : int -- type of rotor interpolation (0 = linear, 1 = Akima splines)
        lift_drag_correction : bool -- include the 2 point lift drag correction?
        max_strain : int -- maximum element strain, before elements are removed from the wake [-]
        maximum_wake_distance : float -- maximum wake distance from the rotor plane (normalized by dia) [-]
        maximum_wake_size : int -- maximum number of wake elements [-]
        n_blades : int -- number of blades (a Structural Model overrides this value)
        n_panels : int -- number of aerodynamic panels per blade (unused if DISCTYPE = 0)
        name : str -- name of the turbine object
        near_wake_length : float -- near wake length [-]
        parameter_path : str -- controller parameter path (leave blank if unused)
        path : str -- path to the .trb file
        polar_discretisation : int -- the polar discretization for the unsteady BEM [-]
        rotational_direction : int -- direction of rotor rotation (0 = STANDARD or 1 = REVERSED)
        rotor_clearance : float -- rotor clearance to ground [m] (VAWT only)
        rotor_cone_angle : float -- rotor cone angle [rad] (HAWT only)
        rotor_configuration : int -- rotor configuration (0 = UPWIND or 1 = DOWNWIND)
        rotor_overhang : float -- rotor overhang [m] (HAWT only)
        rotor_shaft_tilt : float -- shaft tilt angle [rad] (HAWT only)
        rotor_x_tilt : float -- rotor x-tilt angle [rad] (VAWT only)
        rotor_y_tilt : float -- rotor y-tilt angle [rad] (VAWT only)
        shed_vortex : bool -- include shed vortex elements?
        structure_path : str -- path to the structural file (leave blank if unused)
        tf_ate : float -- Tf constant for the ATEFLAP dynamic stall model
        tf_oye : float -- Tf constant for the OYE dynamic stall model
        tower_bot_radius : float -- tower bottom radius [m]
        tower_drag_coefficient : float -- tower drag coefficient [-] (if a Structural Model is used the tower drag is defined in the tower input file)
        tower_height : float -- tower height [m]
        tower_shadow : bool -- include the tower shadow effect
        tower_top_radius : float -- tower top radius [m]
        tp_ate : float -- Tp constant for the ATEFLAP dynamic stall model
        trailing_vortex : bool -- include trailing vortex elements?
        turbine_type : int -- turbine type (0 = HAWT or 1 = VAWT)
        unsteady_aerodynamics : bool -- include unsteady non-circulatory aerodynamics?
        vortex_strain : bool -- calculate vortex strain?
        vortex_viscosity : float -- turbulent vortex viscosity factor [-]
        wake_core_radius : float -- initial core radius of the free wake vortex (fraction of local chord) [-]
        wake_integration_type : int -- wake integration type (0 = EF; 1 = ET, 2 = PC; 3 = PC2B)
        wake_length_type : int -- wake length type (0 = counted in rotor revolutions, 1 = counted in time steps)
        wake_reduction_factor : float -- wake reduction factor [-]
        wake_relaxation_factor : float -- wake relaxation factor [0-1]
        wake_self_induction : bool -- calculate wake self-induction?
        wake_type : int -- wake type (0 = free vortex wake; 1 = unsteady BEM) (unsteady BEM is only available for HAWT)
        zone_1_factor : int -- wake zone 1 streamwise factor [-]
        zone_1_factor_s : int -- wake zone 1 spanwise factor [-]
        zone_1_length : float -- wake zone 1 length [-]
        zone_2_factor : int -- wake zone 2 streamwise factor [-]
        zone_2_factor_s : int -- wake zone 2 spanwise factor [-]
        zone_2_length : float -- wake zone 2 length [-]
        zone_3_factor : int -- wake zone 3 streamwise factor [-]
        zone_3_factor_s : int -- wake zone 3 spanwise factor [-]
        zone_3_length : float -- wake zone 3 length [-]
    """

    def __init__(self, path):
        """
        Parse the Turbine object.
        
        Arguments:
            path: str -- The path to the .trb file.

        Returns:
            None
        """

        self.path = path

        # Open file
        if os.path.isfile(path):
            f = open(path, "r", encoding="utf-8")
        else:
            print(f"No file found at {path}!")
            sys.exit(1)

        # Parse data in file
        # Object name
        self.name = parse(f, "OBJECTNAME", 0, str)

        # Rotor definition
        self.blade_path           = parse(f,     "BLADEFILE", 0, str)
        self.turbine_type         = parse(f,      "TURBTYPE", 0, int)
        self.n_blades             = parse(f,     "NUMBLADES", 0, int)
        self.rotor_configuration  = parse(f,   "ROTORCONFIG", 0, int)
        self.rotational_direction = parse(f, "ROTATIONALDIR", 0, int)
        self.discretisation_type  = parse(f,      "DISCTYPE", 0, int)
        self.interpolation_type   = parse(f,      "INTPTYPE", 0, int)
        self.n_panels             = parse(f,     "NUMPANELS", 0, int)

        # Turbine geometry pParameters
        self.rotor_overhang   = parse(f,    "OVERHANG", 0, float)
        self.rotor_shaft_tilt = parse(f,   "SHAFTTILT", 0, float)
        self.rotor_cone_angle = parse(f,   "ROTORCONE", 0, float)
        self.rotor_clearance  = parse(f,   "CLEARANCE", 0, float)
        self.rotor_x_tilt     = parse(f,       "XTILT", 0, float)
        self.rotor_y_tilt     = parse(f,       "YTILT", 0, float)
        self.tower_height     = parse(f, "TOWERHEIGHT", 0, float)
        self.tower_top_radius = parse(f, "TOWERTOPRAD", 0, float)
        self.tower_bot_radius = parse(f, "TOWERBOTRAD", 0, float)

        # Dynamic stall models
        self.dynamic_stall_type = parse(f, "DYNSTALLTYPE", 0,   int)
        self.tf_oye             = parse(f,       "TF_OYE", 0, float)
        self.am_gb              = parse(f,        "AM_GB", 0, float)
        self.tf_ate             = parse(f,       "TF_ATE", 0, float)
        self.tp_ate             = parse(f,       "TP_ATE", 0, float)
        self.iag_parameters     = None  # Not implemented

        # Aerodynamic models
        self.unsteady_aerodynamics   = parse(f, "UNSTEADYAERO", 0,  bool)
        self.lift_drag_correction    = parse(f,   "2PLIFTDRAG", 0,  bool)
        self.himmelskamp_stall_delay = parse(f,  "HIMMELSKAMP", 0,  bool)
        self.tower_shadow            = parse(f,  "TOWERSHADOW", 0,  bool)
        self.tower_drag_coefficient  = parse(f,    "TOWERDRAG", 0, float)

        # Wake type
        self.wake_type = parse(f, "WAKETYPE", 0, int)

        # Vortex wake parameters (only used if waketype = 0)
        self.wake_integration_type  = parse(f,      "WAKEINTTYPE", 0,   int)
        self.wake_self_induction    = parse(f,       "WAKEROLLUP", 0,  bool)
        self.trailing_vortex        = parse(f,     "TRAILINGVORT", 0,  bool)
        self.shed_vortex            = parse(f,         "SHEDVORT", 0,  bool)
        self.convection_type        = parse(f,   "CONVECTIONTYPE", 0,   int)
        self.wake_relaxation_factor = parse(f,   "WAKERELAXATION", 0, float)
        self.first_wake_row_length  = parse(f,     "FIRSTWAKEROW", 0, float)
        self.maximum_wake_size      = parse(f,      "MAXWAKESIZE", 0,   int)
        self.maximum_wake_distance  = parse(f,      "MAXWAKEDIST", 0, float)
        self.wake_reduction_factor  = parse(f,    "WAKEREDUCTION", 0, float)
        self.wake_length_type       = parse(f,   "WAKELENGTHTYPE", 0,   int)
        self.conversion_length      = parse(f, "CONVERSIONLENGTH", 0, float)
        self.near_wake_length       = parse(f,   "NEARWAKELENGTH", 0, float)
        self.zone_1_length          = parse(f,      "ZONE1LENGTH", 0, float)
        self.zone_2_length          = parse(f,      "ZONE2LENGTH", 0, float)
        self.zone_3_length          = parse(f,      "ZONE3LENGTH", 0, float)
        self.zone_1_factor          = parse(f,      "ZONE1FACTOR", 0,   int)
        self.zone_2_factor          = parse(f,      "ZONE2FACTOR", 0,   int)
        self.zone_3_factor          = parse(f,      "ZONE3FACTOR", 0,   int)
        self.zone_1_factor_s        = parse(f,    "ZONE1FACTOR_S", 0,   int)
        self.zone_2_factor_s        = parse(f,    "ZONE2FACTOR_S", 0,   int)
        self.zone_3_factor_s        = parse(f,    "ZONE3FACTOR_S", 0,   int)

        # Vortex core parameters (only used if waketype = 0)
        self.bound_core_radius = parse(f, "BOUNDCORERADIUS", 0, float)
        self.wake_core_radius  = parse(f,  "WAKECORERADIUS", 0, float)
        self.vortex_viscosity  = parse(f, "VORTEXVISCOSITY", 0, float)
        self.vortex_strain     = parse(f,    "VORTEXSTRAIN", 0,  bool)
        self.max_strain        = parse(f,       "MAXSTRAIN", 0,   int)

        # Gamma iteration parameters (only used if waketype = 0)
        self.gamma_relaxation = parse(f, "GAMMARELAXATION", 0, float)
        self.gamma_epsilon    = parse(f,    "GAMMAEPSILON", 0, float)
        self.gamma_iterations = parse(f, "GAMMAITERATIONS", 0,   int)

        # Unsteady BEM parameters
        self.polar_discretisation = parse(f,  "POLARDISC", 0,   int)
        self.bem_tip_loss         = parse(f, "BEMTIPLOSS", 0,  bool)
        self.bem_speedup          = parse(f, "BEMSPEEDUP", 0, float)

        # Structural model
        self.structure_path          = parse(f, "STRUCTURALFILE", 0,  str)
        self.geometric_stiffness     = parse(f,  "GEOMSTIFFNESS", 0, bool)
        self.aerodynamic_panel_loads = parse(f, "AEROPANELLOADS", 0, bool)

        # Turbine Controller
        self.controller_type = parse(f, "CONTROLLERTYPE", 0, int)
        self.controller_path = parse(f, "CONTROLLERFILE", 0, str)
        self.parameter_path  = parse(f,  "PARAMETERFILE", 0, str)

        # Format parsed data
        self.blade_path      = os.path.join(os.path.dirname(path), self.blade_path.replace("/", "\\"))
        self.structure_path  = os.path.join(os.path.dirname(path), self.structure_path.replace("/", "\\"))
        self.controller_path = os.path.join(os.path.dirname(path), self.controller_path.replace("/", "\\"))
        self.parameter_path  = os.path.join(os.path.dirname(path), self.parameter_path.replace("/", "\\"))

        self.rotor_shaft_tilt = np.radians(self.rotor_shaft_tilt)
        self.rotor_cone_angle = np.radians(self.rotor_cone_angle)
        self.rotor_x_tilt     = np.radians(self.rotor_x_tilt)
        self.rotor_y_tilt     = np.radians(self.rotor_y_tilt)

        # Aero, Structure, and Control objects
        self.aero = Aero(self.blade_path)
        self.structure = None
        self.control = None

        # Close file
        f.close()

if __name__ == "__main__":

    # Parse turbine file
    turbine = Turbine("data\\turbines\\DTU_10MW\\DTU_10MW_RWT.trb")

    # Print turbine data
    print("Turbine Name:", turbine.name)
    print("Blade File Path:", turbine.blade_path)
    print("Turbine Type:", turbine.turbine_type)
    print("Number of Blades:", turbine.n_blades, "[-]")
    print("Rotor Configuration:", turbine.rotor_configuration)
    print("Rotational Direction:", turbine.rotational_direction)
    print("Discretisation Type:", turbine.discretisation_type)
    print("Interpolation Type:", turbine.interpolation_type)
    print("Number of Panels:", turbine.n_panels, "[-]")
    print("Rotor Overhang:", turbine.rotor_overhang, "[m]")
    print("Rotor Shaft Tilt:", turbine.rotor_shaft_tilt, "[rad]")
    print("Rotor Cone Angle:", turbine.rotor_cone_angle, "[rad]")
    print("Rotor Clearance:", turbine.rotor_clearance, "[m]")
    print("Rotor X-Tilt:", turbine.rotor_x_tilt, "[rad]")
    print("Rotor Y-Tilt:", turbine.rotor_y_tilt, "[rad]")
    print("Tower Height:", turbine.tower_height, "[m]")
    print("Tower Top Radius:", turbine.tower_top_radius, "[m]")
    print("Tower Bottom Radius:", turbine.tower_bot_radius, "[m]")
    print("Dynamic Stall Type:", turbine.dynamic_stall_type)
    print("Tf (OYE dynamic stall model):", turbine.tf_oye, "[-]")
    print("Am (GORMONT-BERG dynamic stall model):", turbine.am_gb, "[-]")
    print("Tf (ATEFLAP dynamic stall model):", turbine.tf_ate, "[-]")
    print("Tp :", turbine.tp_ate, "[-]")
    print("IAG Parameters (IAG dynamic stall model):", turbine.iag_parameters)
    print("Unsteady Aerodynamics:", turbine.unsteady_aerodynamics)
    print("Lift Drag Correction:", turbine.lift_drag_correction)
    print("Himmelskamp Stall Delay:", turbine.himmelskamp_stall_delay)
    print("Tower Shadow:", turbine.tower_shadow)
    print("Tower Drag Coefficient:", turbine.tower_drag_coefficient, "[-]")
    print("Wake Type:", turbine.wake_type)
    print("Wake Integration Type:", turbine.wake_integration_type)
    print("Wake Self Induction:", turbine.wake_self_induction)
    print("Trailing Vortex:", turbine.trailing_vortex)
    print("Shed Vortex:", turbine.shed_vortex)
    print("Convection Type:", turbine.convection_type)
    print("Wake Relaxation Factor:", turbine.wake_relaxation_factor, "[-]")
    print("First Wake Row Length:", turbine.first_wake_row_length, "[-]")
    print("Maximum Wake Size:", turbine.maximum_wake_size, "[-]")
    print("Maximum Wake Distance:", turbine.maximum_wake_distance, "[-]")
    print("Wake Reduction Factor:", turbine.wake_reduction_factor, "[-]")
    print("Wake Length Type:", turbine.wake_length_type)
    print("Conversion Length:", turbine.conversion_length, "[-]")
    print("Near Wake Length:", turbine.near_wake_length, "[-]")
    print("Zone 1 Length:", turbine.zone_1_length, "[-]")
    print("Zone 2 Length:", turbine.zone_2_length, "[-]")
    print("Zone 3 Length:", turbine.zone_3_length, "[-]")
    print("Zone 1 Factor:", turbine.zone_1_factor, "[-]")
    print("Zone 2 Factor:", turbine.zone_2_factor, "[-]")
    print("Zone 3 Factor:", turbine.zone_3_factor, "[-]")
    print("Zone 1 Factor S:", turbine.zone_1_factor_s, "[-]")
    print("Zone 2 Factor S:", turbine.zone_2_factor_s, "[-]")
    print("Zone 3 Factor S:", turbine.zone_3_factor_s, "[-]")
    print("Bound Core Radius:", turbine.bound_core_radius, "[-]")
    print("Wake Core Radius:", turbine.wake_core_radius, "[-]")
    print("Vortex Viscosity:", turbine.vortex_viscosity, "[-]")
    print("Vortex Strain:", turbine.vortex_strain)
    print("Max Strain:", turbine.max_strain, "[-]")
    print("Gamma Relaxation:", turbine.gamma_relaxation, "[-]")
    print("Gamma Epsilon:", turbine.gamma_epsilon, "[-]")
    print("Gamma Iterations:", turbine.gamma_iterations, "[-]")
    print("Polar Discretisation:", turbine.polar_discretisation, "[-]")
    print("BEM Tip Loss:", turbine.bem_tip_loss)
    print("BEM Speedup:", turbine.bem_speedup, "[s]")
    print("Structure File Path:", turbine.structure_path)
    print("Geometric Stiffness:", turbine.geometric_stiffness)
    print("Aerodynamic Panel Loads:", turbine.aerodynamic_panel_loads)
    print("Controller Type:", turbine.controller_type)
    print("Controller Path:", turbine.controller_path)
    print("Parameter Path:", turbine.parameter_path)
