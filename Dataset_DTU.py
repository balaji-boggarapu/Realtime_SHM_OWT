import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

file_path = "C:\Projects\DTU Dataset\load_results\load_results\constrained_box_shearShifted_wind_speed_8.5_event_200504081210_seed_1180.bin"
file_size = os.path.getsize(file_path)

# print(f"File size = {file_size}")

data = np.fromfile(file_path, dtype = np.float32).reshape(29, 60000)

channels = [
    "Tower base moment about the x axis (kNm)",
    "Tower base moment about the y axis (kNm)",
    "Tower base moment about the z axis (kNm)",
    "Tower top moment about the x axis (kNm)",
    "Tower top moment about the y axis (kNm)",
    "Tower top moment about the z axis (kNm)",
    "Main shaft main bearing moment about the x axis (kNm)",
    "Main shaft main bearing moment about the y axis (kNm)",
    "Main shaft main bearing moment about the z axis (kNm)",
    "Blade 1 root moment about the x axis (kNm)",
    "Blade 1 root moment about the y axis (kNm)",
    "Blade 1 root moment about the z axis (kNm)",
    "Blade 2 root moment about the x axis (kNm)",
    "Blade 2 root moment about the y axis (kNm)",
    "Blade 2 root moment about the z axis (kNm)",
    "Blade 3 root moment about the x axis (kNm)",
    "Blade 3 root moment about the y axis (kNm)",
    "Blade 3 root moment about the z axis (kNm)",
    "Tower top acceleration in the x direction (m/s^2)",
    "Tower top acceleration in the y direction (m/s^2)",
    "Blade 1 tip deflection in the x direction (m)",
    "Blade 1 tip deflection in the y direction (m)",
    "Blade 1 tip deflection in the z direction (m)",
    "Blade 2 tip deflection in the x direction (m)",
    "Blade 2 tip deflection in the y direction (m)",
    "Blade 2 tip deflection in the z direction (m)",
    "Blade 3 tip deflection in the x direction (m)",
    "Blade 3 tip deflection in the y direction (m)",
    "Blade 3 tip deflection in the z direction (m)"
]

df = pd.DataFrame(data.T, columns = channels)

plt.plot(df.iloc[:,-3])
plt.axhline(df.iloc[:,-3].mean(), color = 'red')
plt.show()
