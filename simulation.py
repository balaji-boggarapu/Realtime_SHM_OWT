import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from structure_parameters import StructureParameters
from force_calculator import ForceCalculator


class TowerSimulation(ForceCalculator):
    def __init__(self):
        super().__init__()
        self.t_span = (0, 15)
        self.time_steps = 101
        self.t_eval = np.linspace(*self.t_span, self.time_steps)
        self.C = 3.563 * (self.c * self.E_tower / self.L_tower ** 3)  # In final equation C = 0.02 (as per the authors)
        self.K = (3.563 * self.E_tower / self.L_tower ** 3) - (
                    1.285 * self.dead_load() / self.L_tower)  # In final equation K = 4.04 (as per the authors)

        self.force_wave_values = np.array([self.force_wave(t) for t in self.t_eval])

        self.time_independent_forces = self.force_wind() + self.force_aero() + self.dead_load()

    # Shape function definition (Two shape functions are considered)
    def shape_function_underwater(self,y):
        return 1.62 * (y / self.L_tower) ** 3 - 3.11 * (y / self.L_tower) ** 2 + 2.49 * (y / self.L_tower)

    def shape_function_abovewater(self,y):
        return 1.62 * (y / self.L_tower) ** 3 - 3.11 * (y / self.L_tower) ** 2 + 2.49 * (y / self.L_tower)

    def shape_function(self,y_vals):
        phi_vals = np.zeros_like(y_vals, dtype=float)

        for i, y in enumerate(y_vals):

            if y < self.d:
                phi_vals[i] = self.shape_function_underwater(y)  # Section under water
            else:
                phi_vals[i] = self.shape_function_abovewater(y)  # Section above water

        return phi_vals

        # Displacement calculation at different heights

    def displacement(self,Z_t, y_vals):
        phi_vals = self.shape_function(y_vals)  # Shape function at heights
        return np.outer(phi_vals, Z_t)

        # Given acceleration equation from authors

    def acceleration_underwater(self,t, y):
        Z, dZ = y
        a = ((self.time_independent_forces + np.interp(t, self.t_eval, self.force_wave_values) - (self.C * dZ) - (self.K * Z)) / self.total_mass())
        return [dZ, a]

    def acceleration_abovewater(self,t, y):
        Z, dZ = y
        a = ((self.time_independent_forces - (self.C * dZ) - (self.K * Z)) / self.total_mass())
        return [dZ, a]

    # Solve the ODE system for Z(t) and Z_dot(t)
    def solve_motion(self,t_span, t_eval, Z0, dZ0, section):
        y0 = [Z0, dZ0]
        if section == "underwater":
            sol = solve_ivp(self.acceleration_underwater, t_span, y0, t_eval=t_eval, method='BDF', rtol=1e-4, atol=1e-6)
        else:
            sol = solve_ivp(self.acceleration_abovewater, t_span, y0, t_eval=t_eval, method='BDF', rtol=1e-4, atol=1e-6)
        return sol.t, sol.y[0], sol.y[1]

    # Calculate the average time period of displacements
    def time_period(self,signal, t, decimals):
        peaks, _ = self.find_peaks(signal)
        peaks_time = t[peaks]
        period = np.diff(peaks_time)
        mean_period = round(np.mean(period), decimals)
        return mean_period

    def run(self):
        y_vals = np.array([20, 35, 50, 65])  # Heights at which displacement is desired

        # Creating gradiant colour maps
        colors = [
            "#2ca02c",
            "#DAA520",
            "#ff7f0e",
            "#d62728"
        ]

        IC = [(0, 0)]  # Set of initial conditions

        # Splitting the y_vals array to above and below sea level
        y_vals_underwater = y_vals[y_vals <= self.d]
        y_vals_abovewater = y_vals[y_vals > self.d]

        plt.figure(figsize=(8, 5))

        for Z0, dZ0 in IC:
            for i, y_b in enumerate(y_vals_underwater):
                t_vals_1, Z_t_1, Z_dot_1 = self.solve_motion(self.t_span, self.t_eval, Z0, dZ0,
                                                        "underwater")  # Soling for section underwater

                displacements_underwater = self.displacement(Z_t_1,
                                                        y_vals_underwater)  # Calculate displacements at different heights

                plt.plot(t_vals_1, displacements_underwater[i], label=f"y = {y_b:.2f} m",
                         color=colors[i])  # Plotting the curves

                column_labels = [f"Z = {Z0}, dZ0 = {dZ0}" for Z0, dZ0 in IC]

            for j, y_a in enumerate(y_vals_abovewater):
                t_vals_2, Z_t_2, Z_dot_2 = self.solve_motion(self.t_span, self.t_eval, Z_t_1[0], Z_dot_1[0],
                                                        "abovewater")  # Solving for section above water

                displacements_abovewater = self.displacement(Z_t_2,
                                                        y_vals_abovewater)  # Calculate displacement at different heights

                plt.plot(t_vals_1, displacements_abovewater[j], label=f"y = {y_a:.2f} m",
                         color=colors[i + j + 1])  # Plotting the curves

        plt.xlabel(r"$t \,\mathrm{[s]}$")
        plt.ylabel(r"$s(y,t) \,\mathrm{[m]}$")
        plt.legend(loc='upper right')
        plt.grid()
        plt.show()