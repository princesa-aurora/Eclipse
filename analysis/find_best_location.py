import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import PIL

from utils import *



def time_visible(occultation_data, times):
    # for every point compute the total timespan for which the eclipse is visible

    start_idcs = np.argmax(occultation_data > 0, axis=0)
    end_idcs = len(occultation_data)-1 - np.argmax(occultation_data[::-1] > 0, axis=0)

    start_times = times[start_idcs]
    end_times = times[end_idcs]

    time = end_times - start_times

    # the above calculation returns nonsense at locations where the eclipse isn't visible at any time
    # so we mask those locations and set the time manually to 0
    time[np.all(occultation_data == 0, axis=0)] = 0

    return time



def time_in_totality(occultation_data, times):
    # for every point compute the timespan for which the eclipse is total

    start_idcs = np.argmax(occultation_data == 1, axis=0)
    end_idcs = len(occultation_data)-1 - np.argmax(occultation_data[::-1] == 1, axis=0)

    start_times = times[start_idcs]
    end_times = times[end_idcs]

    time = end_times - start_times

    # the above calculation returns nonsense at locations where the eclipse isn't total at any time
    # so we mask those locations and set the time manually to 0
    time[np.all(occultation_data < 1, axis=0)] = 0

    return time



file_path = "/home/aurora/eclipse_data/22.05.2026_20:12:50/solar_eclipse_of_5.2.2000.nc"


dataset = xr.open_dataset(file_path)

times = dataset["time"].values
sun_lon = dataset["lon_sun"].values
sun_lat = dataset["lat_sun"].values
num_steps = len(times)

lon_grid = dataset["lon_grid"].values
lat_grid = dataset["lat_grid"].values
grid_size = len(lat_grid)

occultation_data = dataset["occultation_data"].values

dataset.close()

land_mask = is_on_land(lon_grid, lat_grid)


visible_times = time_visible(occultation_data, times)
totality_times = time_in_totality(occultation_data, times)

if (visible_times > 0).any():
    argmax_vis = np.argmax(visible_times)
    lon_opt_vis = lon_grid[argmax_vis]
    lat_opt_vis = lat_grid[argmax_vis]
else:
    lon_opt_vis = np.nan
    lat_opt_vis = np.nan

if (visible_times[land_mask] > 0).any():
    argmax_vis_land = np.argmax(visible_times[land_mask])
    lon_opt_vis_land = lon_grid[land_mask][argmax_vis_land]
    lat_opt_vis_land = lat_grid[land_mask][argmax_vis_land]
else:
    lon_opt_vis_land = np.nan
    lat_opt_vis_land = np.nan

if (totality_times > 0).any():
    argmax_tot = np.argmax(totality_times)
    lon_opt_tot = lon_grid[argmax_tot]
    lat_opt_tot = lat_grid[argmax_tot]
else:
    lon_opt_tot = np.nan
    lat_opt_tot = np.nan

if (totality_times[land_mask] > 0).any():
    argmax_tot_land = np.argmax(totality_times[land_mask])
    lon_opt_tot_land = lon_grid[land_mask][argmax_tot_land]
    lat_opt_tot_land = lat_grid[land_mask][argmax_tot_land]
else:
    lon_opt_tot_land = np.nan
    lat_opt_tot_land = np.nan


if np.isnan(lat_opt_vis):
    print("There is no eclipse to be found in the data.")
    exit(0)
else:
    print(f"The eclipse is visible for the longest at lon({lon_opt_vis*180/np.pi}°), lat({lat_opt_vis*180/np.pi}°) anywhere on earth.")
    print(f"It can be seen there for {visible_times[argmax_vis] /60} minutes.")
    print()

if np.isnan(lat_opt_vis_land):
    print("However it isn't visible anywhere on land.")
    exit(0)
else:
    print(f"On land it is visible for the longest at lon({lon_opt_vis_land*180/np.pi}°), lat({lat_opt_vis_land*180/np.pi}°).")
    print(f"It can be seen there for {visible_times[land_mask][argmax_vis_land] /60} minutes.")
    print()

if np.isnan(lat_opt_tot):
    print("The eclipse isn't total anywhere.")
    exit(0)
else:
    print(f"Totality lasts the longest at: lon({lon_opt_tot*180/np.pi}°), lat({lat_opt_tot*180/np.pi}°) anywhere on earth.")
    print(f"It lasts there for {totality_times[argmax_tot] /60} minutes.")
    print()

if np.isnan(lat_opt_tot_land):
    print("However it isn't total anywhere on land.")
    exit(0)
else:
    print(f"On land totality lasts the longest at lon({lon_opt_tot_land*180/np.pi}°), lat({lat_opt_tot_land*180/np.pi}°).")
    print(f"It lasts there for {totality_times[land_mask][argmax_tot_land] /60} minutes.")
    print()




