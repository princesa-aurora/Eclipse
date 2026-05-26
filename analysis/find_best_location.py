import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import PIL

from utils import *


file_path = "/home/aurora/eclipse_data/25.05.2026_04:24:13/solar_eclipse_of_5.2.2000.nc"


dataset = xr.open_dataset(file_path)

times = dataset["time"].values
sun_lon = dataset["lon_sun"].values
sun_lat = dataset["lat_sun"].values

lon_grid = dataset["lon_grid"].values
lat_grid = dataset["lat_grid"].values

classification_data = dataset["classification_data"].values
topology_data = np.mod(classification_data, 4)
angle_data = (classification_data - topology_data) /4

dataset.close()

day_mask = (np.cos(lat_grid[None,:])*np.cos(sun_lat[:,None])*np.cos(lon_grid[None,:]-sun_lon[:,None]) + np.sin(lat_grid[None,:])*np.sin(sun_lat[:,None]) > 0)
land_mask = is_on_land(lon_grid, lat_grid)


# for every point compute the total timespan for which the eclipse is visible
visibility_mask = (topology_data > 0)
visibility_mask = np.logical_and(visibility_mask, day_mask)
start_idcs = np.argmax(visibility_mask, axis=0)
end_idcs = len(visibility_mask)-1 - np.argmax(visibility_mask[::-1], axis=0)
visibility_times = times[end_idcs] - times[start_idcs]
# the above calculation returns nonsense at locations where the eclipse isn't visible at any time
# so we mask those locations and set the time manually to 0
visibility_times[np.all(~visibility_mask, axis=0)] = 0

# for every point compute the timespan for which the eclipse is annular
annularity_mask = (topology_data == 2)
annularity_mask = np.logical_and(annularity_mask, day_mask)
start_idcs = np.argmax(annularity_mask, axis=0)
end_idcs = len(annularity_mask)-1 - np.argmax(annularity_mask[::-1], axis=0)
annularity_times = times[end_idcs] - times[start_idcs]
# the above calculation returns nonsense at locations where the eclipse isn't annular at any time
# so we mask those locations and set the time manually to 0
annularity_times[np.all(~annularity_mask, axis=0)] = 0

# for every point compute the timespan for which the eclipse is total
totality_mask = (topology_data == 3)
totality_mask = np.logical_and(totality_mask, day_mask)
start_idcs = np.argmax(totality_mask, axis=0)
end_idcs = len(totality_mask)-1 - np.argmax(totality_mask[::-1], axis=0)
totality_times = times[end_idcs] - times[start_idcs]
# the above calculation returns nonsense at locations where the eclipse isn't total at any time
# so we mask those locations and set the time manually to 0
totality_times[np.all(~totality_mask, axis=0)] = 0


if (visibility_times == 0).all():
    print(f"{BOLD}There is no eclipse to be found in the data.{RESET}")
    exit(0)

print(f"{BOLD}Best eclipse viewing locations:{RESET}")
print(f"{BOLD}---------------------------------{RESET}")

print(f"{BOLD}Visibility anywhere on earth:{RESET}")
argmax_vis = np.argmax(visibility_times)
lon_opt_vis = lon_grid[argmax_vis]
lat_opt_vis = lat_grid[argmax_vis]
print(f"lon({lon_opt_vis*180/np.pi}°), lat({lat_opt_vis*180/np.pi}°) ({get_location_name(lon_opt_vis, lat_opt_vis)})")
print(f"It is visible there for {visibility_times[argmax_vis] /60} minutes.")
print()

print(f"{BOLD}Annularity anywhere on earth:{RESET}")
if (annularity_times > 0).any():
    argmax_ann = np.argmax(annularity_times)
    lon_opt_ann = lon_grid[argmax_ann]
    lat_opt_ann = lat_grid[argmax_ann]
    print(f"lon({lon_opt_ann*180/np.pi}°), lat({lat_opt_ann*180/np.pi}°) ({get_location_name(lon_opt_ann, lat_opt_ann)})")
    print(f"It is annular there for {annularity_times[argmax_ann] /60} minutes.")
else:
    print("The eclipse isn't annular anywhere on earth.")
print()

print(f"{BOLD}Totality anywhere on earth:{RESET}")
if (totality_times > 0).any():
    argmax_tot = np.argmax(totality_times)
    lon_opt_tot = lon_grid[argmax_tot]
    lat_opt_tot = lat_grid[argmax_tot]
    print(f"lon({lon_opt_tot*180/np.pi}°), lat({lat_opt_tot*180/np.pi}°) ({get_location_name(lon_opt_tot, lat_opt_tot)})")
    print(f"It is total there for {totality_times[argmax_tot] /60} minutes.")
else:
    print("The eclipse isn't total anywhere on earth.")
print()


if (visibility_times[land_mask] == 0).all():
    print(f"{BOLD}The eclipse isn't visible anywhere on land.{RESET}")
    exit(0)

print(f"{BOLD}Visibility anwhere on land:{RESET}")
argmax_vis_land = np.argmax(visibility_times[land_mask])
lon_opt_vis_land = lon_grid[land_mask][argmax_vis_land]
lat_opt_vis_land = lat_grid[land_mask][argmax_vis_land]
print(f"lon({lon_opt_vis_land*180/np.pi}°), lat({lat_opt_vis_land*180/np.pi}°) ({get_location_name(lon_opt_vis_land, lat_opt_vis_land)})")
print(f"It is visible there for {visibility_times[land_mask][argmax_vis_land] /60} minutes.")
print()

print(f"{BOLD}Annularity anywhere on land:{RESET}")
if (annularity_times[land_mask] > 0).any():
    argmax_ann_land = np.argmax(annularity_times[land_mask])
    lon_opt_ann_land = lon_grid[land_mask][argmax_ann_land]
    lat_opt_ann_land = lat_grid[land_mask][argmax_ann_land]
    print(f"lon({lon_opt_ann_land*180/np.pi}°), lat({lat_opt_ann_land*180/np.pi}°) ({get_location_name(lon_opt_ann_land, lat_opt_ann_land)})")
    print(f"It is annular there for {annularity_times[land_mask][argmax_ann_land] /60} minutes.")
else:
    print("The eclipse isn't annular anywhere on land.")
print()

print(f"{BOLD}Totality anywhere on land:{RESET}")
if (totality_times[land_mask] > 0).any():
    argmax_tot_land = np.argmax(totality_times[land_mask])
    lon_opt_tot_land = lon_grid[land_mask][argmax_tot_land]
    lat_opt_tot_land = lat_grid[land_mask][argmax_tot_land]
    print(f"lon({lon_opt_tot_land*180/np.pi}°), lat({lat_opt_tot_land*180/np.pi}°) ({get_location_name(lon_opt_tot_land, lat_opt_tot_land)})")
    print(f"It is total there for {totality_times[land_mask][argmax_tot_land] /60} minutes.")
else:
    print("The eclipse isn't total anywhere on land.")
print()




