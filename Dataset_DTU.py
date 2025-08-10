import numpy as np
import os
import pandas as pd
from sqlalchemy  import create_engine
import urllib.parse

# Fetching the load data of a randon timestamp from the dataset

file_path_load = "C:\Projects\DTU Dataset\load_results\load_results\constrained_box_shearShifted_wind_speed_8.5_event_200504081210_seed_1180.bin"


data = np.fromfile(file_path_load, dtype = np.float32).reshape(29, 60000)

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



# Fetching the load data of a randon timestamp from the dataset

file_path_wind = r"C:\Projects\DTU Dataset\Hov_au100max_ts_uv\wind_data_from_Mark\Timeseries_filtered_O2butter_0.1Hz\Hov_au100max_ts_uv\Hov_au100max_uv100,160_U8_200504081210.csv"

data_wind_raw = pd.read_csv(file_path_wind, header=None)
data_wind_df = data_wind_raw.T
columns_wind = [
    "Stream-wise velocity in m/s at 100m",
    "Stream-wise velocity in m/s at 160m",
    "Transverse velocity in m/s at 100m",
    "Transverse velocity in m/s at 160m"
]

data_wind_df.columns = columns_wind

# Downsampling the load data from 100Hz to 1Hz to match the wind speed data
load_df = pd.DataFrame(data.T, columns=channels)
load_df_downsampled = load_df.iloc[::100]
print(data_wind_df.shape)
print(load_df_downsampled.shape)

# Joining the dataframes
main_df = pd.concat([load_df_downsampled.reset_index(drop=True),data_wind_df.reset_index(drop=True)], axis=1)
print(load_df_downsampled.shape)

'''# Pushing the code to MySQL
password = urllib.parse.quote_plus("Balaji@1998") # To encode the @ character
engine = create_engine(f"mysql+mysqlconnector://balaji_user:{password}@127.0.0.1:3306/wind_db") # Pushing both the dataframes to MySQL server

# Pushing the dataframes to MySQL server
data_wind.to_sql(name='wind_measurements', con=engine, if_exists='replace', index=False)
print(engine.url)'''