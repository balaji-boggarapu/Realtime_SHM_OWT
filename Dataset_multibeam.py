import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.integrate import solve_ivp
from scipy.integrate import quad
from scipy.signal import find_peaks
import pandas as pd

# Parameters of the tower
r_base  = 2 # Radius at the base of the tower [in m]
r_top = 0.75 # Radius at the top of the tower [in m]
d_base = 2 * r_base # Base diameter of the tower [in m]
d_top = 2 * r_top # Top diameter of the tower [in m]
L_tower = 65 # Height of the tower [in m]
c_d = 0.47 # Drag co-efficient for sphere [dimensionless]
c_aero = 0.5 # Aerodynamic co-efficient [dimensionless]
rho_tower = 78500 # Density of material of tower [in kg/m^3]
E_tower = 2.1 * 10*11 # Young's modulus of tower [in N/m^2]
c = 0.005 # Damping co-efficient of system [dimensionless]
# X = 0 # Position of tower in co-ordinate system 

# Parameter of water 
H = 1.80 # Height of water [in m]
d = 20 # Depth of water [in m]
T = 5 # Time [in s]
rho_water = 1025 # Density of sea water [in kg/m^3]
c_d = 2.4 # Drag force co-efficient [dimensionless]
c_m = 0.7 # Inertial force co-efficient [dimensionless]
L_wave = 38.90 # Wave length of water [in m]
T_wave = 5 # Time [in s]

# Parameters of air
rho_air = 1.22 # Density of air [in kg/m^3]
v_ref = 3.58 # Mean velocity of air [in m/s]
z_r = 10 # Reference height [in m] (From 31)
k_r = 0.17 # Terrian factor (from paper)
z_0 = 0.01 # SurLength of roughness [in m]

# v_air = 18 # Mean velocity of air [in m/s]

# Parameters of RNC assembly
m_rna = 83000 # Mass of the RNC assembly [in kg]
g = 9.81 # Aceeleration due to gravity [in m/s^2]
r_hub = 3 # Radius of concentrated mass [in m]
E_blade = 65 * 10**9 # Young's modulus of blade material [in N/m^2]
rho_rnc = 2.14 * 10*3 # Density of RNC [in kg/m^3]
L_blade = 24 # Blade length [in m]
d_hub = 2 * r_hub # Diameter of the hub [in m]
A_swept = 624 # Swept area of blades [in m^2]

# Parameters for equation of motion
c = 0.005 # Damping coefficient [dimensionless]
E_tower = 2.1 * 10**11 # Youngs modulus of tower [in N/m^2]
rho_tower = 8 * 10**3 # Density of tower [in kg/m^3]

# Derived Parameters
k = (2 * np.pi) / L_wave  # Wave number [in m^-1]
omega = (2 * np.pi) / T  # Angular frequency [in Hz]
a = H / 2  # Wave amplitude [in m]

# Mass of the tower
def tower_mass():
    V_tower = (np.pi * L_tower * ((r_base**2) + (r_base*r_top) + (r_top**2))) / 3
    return (rho_tower * V_tower)

# Total mass of the system
def total_mass():
    return tower_mass() + m_rna

# Axial force due to concentrated mass
def dead_load():
    F_axial = m_rna * g
    return F_axial

# Calculating the force on concentrated mass due to wind
def force_wind():
    
    # Calculate the reference area
    A = np.pi * (r_hub**2)
    v_air = v_ref * (np.log(L_tower / z_0) / np.log(z_r / z_0)) #v_ref * k_r * np.log(L_tower / z_r)
    F_wind = 0.5 * c_d * rho_air * A * (v_air**2)
    return F_wind

# Calculating the aerodynamic force
def force_aero (): 
    
    # Defining the integrad function
    def integrand_wind(y):
        r_y = r_base + ((r_top - r_base) * y / L_tower)  # Radius at height y (Linear interpolation)
        v_air = v_ref * (np.log(y / z_0) / np.log(z_r / z_0)) # From [31] # v_ref * k_r * np.log(y / z_r) (This is what authors have used)
        A_y = 2 * np.pi * r_y # Circumference at height y (upon integrating gives the exposed area)
        return 0.5 * c_aero * rho_air * A_y * (v_air**2)
    
    # Calculating the wind force by integrating the integrand from d to L
    F_aero = quad (integrand_wind, d, L_tower)
    return F_aero[0]

# Wave elevation function
def wave_elevation(x, t):
    return a * np.cos(k * x - omega * t)

# Calculating the force due to waves
def force_wave(t):
    
    # Wave height eveluation
    nta = wave_elevation(0,t)
    
    # Calculating the drag force
    def integrand_drag(y, t):
        d_y = d_base + ((d_top - d_base) * y / L_tower)  # Diameter at height y
        u = ((np.pi * H) / T) * np.cosh(k * (d + y)) * np.cos(k * y - omega * t) / np.sinh(k * d)  # Particle velocity
        return 0.5 * rho_water * c_d * d_y * u * abs(u)
    
    # Calculating the drag force by integrating the integrand over -d to nta
    F_d = quad(integrand_drag, -d, nta, args=(t,))
    
    def integrand_inertia(y,t):
        d_y = d_base + ((d_top - d_base) * y / L_tower) # Diameter at the height y
        u_dot = (-(2 * np.pi**2 * H)/ T**2) * np.cosh(k * (d+y)) * np.cos(k * y - omega * t) / np.sinh(k * d) # Particle acceleration from linear wave theory
        return 0.25 * np.pi * rho_water *  c_m * d_y**2 * u_dot 
    
    # Calculating the inertial force by integrating the integrand over -d to nta
    F_i = quad(integrand_inertia, -d, nta, args=(t,))
    
    # From Morisson's equation
    F_wave = F_d[0] + F_i[0]
     
    return F_wave

# Equation of motion
t_span = (0, 15)
time_steps = 101
t_eval = np.linspace(*t_span, time_steps)

C = 3.563 * (c * E_tower/ L_tower**3) #In final equation C = 0.02 (as per the authors)
K = (3.563 * E_tower / L_tower**3) - (1.285 * dead_load()/L_tower) #In final equation K = 4.04 (as per the authors)

force_wave_values = np.array([force_wave(t) for t in t_eval])

time_independent_forces = force_wind() + force_aero() + dead_load()

# Shape function definition (Two shape functions are considered)
def shape_function_underwater(y):
    return 1.62 * (y / L_tower)**3 - 3.11 * (y / L_tower)**2 + 2.49 * (y / L_tower)

def shape_function_abovewater(y):
    return 1.62 * (y / L_tower)**3 - 3.11 * (y / L_tower)**2 + 2.49 * (y / L_tower)

def shape_function(y_vals):
    phi_vals = np.zeros_like(y_vals, dtype = float)
    
    for i, y in enumerate(y_vals):
        
        if y < d:
            phi_vals[i] = shape_function_underwater(y) # Section under water
        else:
            phi_vals[i] = shape_function_abovewater(y) # Section above water
    
    return phi_vals        

   
# Displacement calculation at different heights
def displacement(Z_t, y_vals):
    phi_vals = shape_function(y_vals)  # Shape function at heights
    return np.outer(phi_vals, Z_t)    

# Given acceleration equation from authors
def acceleration_underwater(t, y):
    Z, dZ = y
    a = ((time_independent_forces + np.interp(t, t_eval, force_wave_values) - (C * dZ) - (K * Z)) / total_mass())
    return [dZ, a]

def acceleration_abovewater(t, y):
    Z, dZ = y
    a = ((time_independent_forces - (C * dZ) - (K * Z)) / total_mass())
    return [dZ, a]

# Solve the ODE system for Z(t) and Z_dot(t)
def solve_motion(t_span, t_eval, Z0, dZ0, section):
    y0 = [Z0, dZ0]
    if section == "underwater":
        sol = solve_ivp(acceleration_underwater, t_span, y0, t_eval=t_eval, method = 'BDF', rtol = 1e-4, atol = 1e-6)
    else:
        sol = solve_ivp(acceleration_abovewater, t_span, y0, t_eval=t_eval, method = 'BDF', rtol = 1e-4, atol = 1e-6)
    return sol.t, sol.y[0], sol.y[1]

# Calculate the average time period of displacements
def time_period(signal, t, decimals):
    peaks, _ = find_peaks(signal)
    peaks_time = t[peaks]
    period = np.diff(peaks_time)
    mean_period = round(np.mean(period), decimals)
    return mean_period

def IC(M):
    np.random.seed(0)
    initial_conditions = np.random.rand(M, 2)
    return [tuple(pair) for pair in initial_conditions]

y_vals = np.array([20, 35, 50, 65]) # Heights at which displacement is desired

IC_vals = IC(1000) # Set of initial conditions

# Splitting the y_vals array to above and below sea level
y_vals_underwater = y_vals[y_vals <= d]
y_vals_abovewater = y_vals[y_vals > d]

system_matrix_underwater = []
system_matrix_abovewater = []
IC_counter = 0

for Z0, dZ0 in IC_vals:  
    IC_counter+=1
    print(IC_counter)
    t_vals_1, Z_t_1, Z_dot_1 = solve_motion(t_span, t_eval, Z0, dZ0, "underwater") # Solving for section underwater
    displacements_underwater = displacement(Z_t_1, y_vals_underwater) # Calculate displacements at different heights
    displacements_underwater_flattened = displacements_underwater.flatten()
    row_displacements_underwater = displacements_underwater_flattened.transpose()
    system_matrix_underwater.append(row_displacements_underwater)
        
    t_vals_2, Z_t_2, Z_dot_2 = solve_motion(t_span, t_eval,  Z_t_1[0] ,Z_dot_1[0] , "abovewater") # Solving for section above water
    displacements_abovewater = displacement(Z_t_2, y_vals_abovewater) # Calculate displacement at different heights
    displacements_abovewater_flattened = displacements_abovewater.flatten()
    row_displacements_abovewater = displacements_abovewater_flattened.transpose()
    system_matrix_abovewater.append(row_displacements_abovewater)
        
system_matrix_underwater = np.column_stack(system_matrix_underwater)
system_matrix_abovewater = np.column_stack(system_matrix_abovewater)

system_matrix = np.vstack((system_matrix_underwater, system_matrix_abovewater))
# np.save("Displacement data from OWT.nby", system_matrix)