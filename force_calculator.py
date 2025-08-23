import numpy as np
from scipy.integrate import quad
from structure_parameters import StructureParameters
import random

class ForceCalculator(StructureParameters):
    def __init__(self):
        super().__init__()
        self.val = random.uniform(3,6)


    # Calculating the force on concentrated mass due to wind
    def force_wind(self):
        # Calculate the reference area
        A = np.pi * (self.r_hub ** 2)
        v_air = self.v_ref * (np.log(self.L_tower / self.z_0) / np.log(self.z_r / self.z_0))  # v_ref * k_r * np.log(L_tower / z_r)
        F_wind = 0.5 * self.c_d * self.rho_air * A * (v_air ** 2) * self.val
        return F_wind

    # Calculating the aerodynamic force
    def force_aero(self):
        # Defining the integrad function
        def integrand_wind(y):
            r_y = self.r_base + ((self.r_top - self.r_base) * y / self.L_tower)  # Radius at height y (Linear interpolation)
            v_air = self.v_ref * (np.log(y / self.z_0) / np.log(
                self.z_r / self.z_0))  # From [31] # v_ref * k_r * np.log(y / z_r) (This is what authors have used)
            A_y = 2 * np.pi * r_y  # Circumference at height y (upon integrating gives the exposed area)
            return 0.5 * self.c_aero * self.rho_air * A_y * (v_air ** 2) * self.val

        # Calculating the wind force by integrating the integrand from d to L
        F_aero = quad(integrand_wind, self.d, self.L_tower)
        return F_aero[0]

    # Wave elevation function
    def wave_elevation(self, x, t):
        return self.a * np.cos(self.k * x - self.omega * t)

    # Calculating the force due to waves
    def force_wave(self, t):
        # Wave height eveluation
        nta = self.wave_elevation(0, t)

        # Calculating the drag force
        def integrand_drag(y, t):
            d_y = self.d_base + ((self.d_top - self.d_base) * y / self.L_tower)  # Diameter at height y
            u = ((np.pi * self.H) / self.T) * np.cosh(self.k * (self.d + y)) * np.cos(self.k * y - self.omega * t) / np.sinh(
                self.k * self.d)  # Particle velocity
            return 0.5 * self.rho_water * self.c_d * d_y * u * abs(u)

        # Calculating the drag force by integrating the integrand over -d to nta
        F_d = quad(integrand_drag, -self.d, nta, args=(t,))

        def integrand_inertia(y, t):
            d_y = self.d_base + ((self.d_top - self.d_base) * y / self.L_tower)  # Diameter at the height y
            u_dot = (-(2 * np.pi ** 2 * self.H) / self.T ** 2) * np.cosh(self.k * (self.d + y)) * np.cos(self.k * y - self.omega * t) / np.sinh(
                self.k * self.d)  # Particle acceleration from linear wave theory
            return 0.25 * np.pi * self.rho_water * self.c_m * d_y ** 2 * u_dot

            # Calculating the inertial force by integrating the integrand over -d to nta

        F_i = quad(integrand_inertia, -self.d, nta, args=(t,))

        # From Morisson's equation

        F_wave = F_d[0] + F_i[0]

        return F_wave