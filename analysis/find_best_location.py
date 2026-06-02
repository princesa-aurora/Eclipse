import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from tqdm import tqdm

from utils import *


file_path = "/data/eclipse_data/28.05.2026_18-36-25/solar_eclipse_of_05.02.2000.nc"


dataset = xr.open_dataset(file_path, chunks={"steps": 1, "grid": -1})

times = dataset["time"].values
sun_lon = dataset["lon_sun"].values
sun_lat = dataset["lat_sun"].values
num_steps = times.shape[0]
dt = (times[-1]-times[0])/(num_steps-1)

lon_grid = dataset["lon_grid"].values
lat_grid = dataset["lat_grid"].values
grid_size = lat_grid.shape[0]

classification_data = dataset["classification_data"].data

dataset.close()


visibility_start_idcs = np.full(lat_grid.shape, num_steps, dtype=int)
visibility_end_idcs = np.full(lat_grid.shape, num_steps, dtype=int)
annularity_start_idcs = np.full(lat_grid.shape, num_steps, dtype=int)
annularity_end_idcs = np.full(lat_grid.shape, num_steps, dtype=int)
totality_start_idcs = np.full(lat_grid.shape, num_steps, dtype=int)
totality_end_idcs = np.full(lat_grid.shape, num_steps, dtype=int)

for i in tqdm(range(num_steps), desc="Processing steps"):
    classification_chunk = classification_data[i].compute()
    topology_chunk = np.mod(classification_chunk, 4)

    # compute day mask
    day_mask = (np.cos(lat_grid)*np.cos(sun_lat[i])*np.cos(lon_grid-sun_lon[i]) + np.sin(lat_grid)*np.sin(sun_lat[i]) > 0)

    # visibility:
    # compute visibility mask (both eclipse topology must be nonzero and it must be day)
    visibility_mask = (topology_chunk > 0)
    visibility_mask = np.logical_and(visibility_mask, day_mask)
    # find start and end indices:
    # end_idcs gets set until the last True (per location), so ultimately its value is the index of the last True
    # start_idcs gets set until the first True (per location), so ultimately its value is the index of the first True
    # (this is since then end_idcs gets set for the first time and is not num_steps any more in the next step)
    visibility_start_idcs[visibility_end_idcs==num_steps] = i
    visibility_end_idcs[visibility_mask] = i

    # annularity
    annularity_mask = (topology_chunk == 2)
    annularity_mask = np.logical_and(annularity_mask, day_mask)
    annularity_start_idcs[annularity_end_idcs==num_steps] = i
    annularity_end_idcs[annularity_mask] = i

    # totality
    totality_mask = (topology_chunk == 2)
    totality_mask = np.logical_and(totality_mask, day_mask)
    totality_start_idcs[totality_end_idcs==num_steps] = i
    totality_end_idcs[totality_mask] = i


# Some comments on edge cases:
# when the eclipse is visible/annular/total already at the first step (at some location):
# start_idx = 0, end_idx = last => all good
# when the eclipse is visible/annular/total until the last step (at some location):
# start_idx = first, end_idx = num_steps-1 => all good
# when the eclipse is never visible/annular/total (at some location):
# start_idx = num_steps-1, end_idx = num_steps => end_idx out ouf range
# => need to manually mask those positions and set the times to 0

# note the +dt is since end_idx is the last one where the eclipse is visible/annular/total,
# not the first one after that as is usual in programming
visibility_times = np.zeros(lat_grid.shape)
visibility_times[visibility_end_idcs!=num_steps] = times[visibility_end_idcs[visibility_end_idcs!=num_steps]] - times[visibility_start_idcs[visibility_end_idcs!=num_steps]] +dt

annularity_times = np.zeros(lat_grid.shape)
annularity_times[annularity_end_idcs!=num_steps] = times[annularity_end_idcs[annularity_end_idcs!=num_steps]] - times[annularity_start_idcs[annularity_end_idcs!=num_steps]] +dt

totality_times = np.zeros(lat_grid.shape)
totality_times[totality_end_idcs!=num_steps] = times[totality_end_idcs[totality_end_idcs!=num_steps]] - times[totality_start_idcs[totality_end_idcs!=num_steps]] +dt


land_mask = is_on_land(lon_grid, lat_grid)


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




