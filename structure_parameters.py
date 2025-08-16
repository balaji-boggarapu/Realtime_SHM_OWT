import numpy as np
class StructureParameters:
    def __init__(self):
        # Parameters of the tower
        self.r_base = 2  # Radius at the base of the tower [in m]
        self.r_top = 0.75  # Radius at the top of the tower [in m]
        self.d_base = 2 * self.r_base  # Base diameter of the tower [in m]
        self.d_top = 2 * self.r_top  # Top diameter of the tower [in m]
        self.L_tower = 65  # Height of the tower [in m]
        self.c_d = 0.47  # Drag co-efficient for sphere [dimensionless]
        self.c_aero = 0.5  # Aerodynamic co-efficient [dimensionless]
        self.rho_tower = 78500  # Density of material of tower [in kg/m^3]
        self.E_tower = 2.1 * 10 * 11  # Young's modulus of tower [in N/m^2]
        self.c = 0.005  # Damping co-efficient of system [dimensionless]
        # X = 0 # Position of tower in co-ordinate system

        # Parameter of water
        self.H = 1.80  # Height of water [in m]
        self.d = 20  # Depth of water [in m]
        self.T = 5  # Time [in s]
        self.rho_water = 1025  # Density of sea water [in kg/m^3]
        self.c_d = 2.4  # Drag force co-efficient [dimensionless]
        self.c_m = 0.7  # Inertial force co-efficient [dimensionless]
        self.L_wave = 38.90  # Wave length of water [in m]
        self.T_wave = 5  # Time [in s]

        # Parameters of air
        self.rho_air = 1.22  # Density of air [in kg/m^3]
        self.v_ref = 3.58  # Mean velocity of air [in m/s]
        self.z_r = 10  # Reference height [in m] (From 31)
        self.k_r = 0.17  # Terrian factor (from paper)
        self.z_0 = 0.01  # SurLength of roughness [in m]

        # v_air = 18 # Mean velocity of air [in m/s]

        # Parameters of RNC assembly
        self.m_rna = 83000  # Mass of the RNC assembly [in kg]
        self.g = 9.81  # Aceeleration due to gravity [in m/s^2]
        self.r_hub = 3  # Radius of concentrated mass [in m]
        self.E_blade = 65 * 10 ** 9  # Young's modulus of blade material [in N/m^2]
        self.rho_rnc = 2.14 * 10 * 3  # Density of RNC [in kg/m^3]
        self.L_blade = 24  # Blade length [in m]
        self.d_hub = 2 * self.r_hub  # Diameter of the hub [in m]
        self.A_swept = 624  # Swept area of blades [in m^2]

        # Parameters for equation of motion
        self.c = 0.005  # Damping coefficient [dimensionless]
        self.E_tower = 2.1 * 10 ** 11  # Youngs modulus of tower [in N/m^2]
        self.rho_tower = 8 * 10 ** 3  # Density of tower [in kg/m^3]

        # Derived Parameters
        self.k = (2 * np.pi) / self.L_wave  # Wave number [in m^-1]
        self.omega = (2 * np.pi) / self.T  # Angular frequency [in Hz]
        self.a = self.H / 2  # Wave amplitude [in m]

    # Mass of the tower
    def tower_mass(self):
        V_tower = (np.pi * self.L_tower * ((self.r_base ** 2) + (self.r_base * self.r_top) + (self.r_top ** 2))) / 3
        return self.rho_tower * V_tower

    # Total mass of the system
    def total_mass(self):
        return self.tower_mass() + self.m_rna

    # Axial force due to concentrated mass
    def dead_load(self):
        F_axial = self.m_rna * self.g
        return F_axial
