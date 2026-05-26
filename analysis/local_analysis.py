import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import PIL

from utils import *


file_path = "/home/aurora/eclipse_data/25.05.2026_04:24:13/solar_eclipse_of_5.2.2000.nc"

lon_observ_deg = 0
lat_observ_deg = -90


dataset = xr.open_dataset(file_path)

times = dataset["time"].values
sun_lon = dataset["lon_sun"].values
sun_lat = dataset["lat_sun"].values
num_steps = len(times)
minutes_after_start = (times-times[0])/60

lon_grid = dataset["lon_grid"].values
lat_grid = dataset["lat_grid"].values
grid_size = len(lat_grid)

occultation_data = dataset["occultation_data"].values
classification_data = dataset["classification_data"].values
topology_data = np.mod(classification_data, 4)
angle_data = (classification_data - topology_data) /4 *6

dataset.close()


lon_observ = lon_observ_deg * np.pi/180
lat_observ = lat_observ_deg * np.pi/180
location_name = get_location_name(lon_observ, lat_observ)
# check for day (since if night the eclipse isn't visible)
local_day_mask = (np.cos(lat_observ)*np.cos(sun_lat)*np.cos(lon_observ-sun_lon) + np.sin(lat_observ)*np.sin(sun_lat) > 0)

# get the index of the grid point closest to the chosen location
# i.e. an interpolation of type "nearest" will be done
nearest_idx = get_nearest_idcs(lon_grid, lat_grid,
                                lon_observ, lat_observ, k=1)

local_occults = occultation_data[:, nearest_idx] *100 # local occultation rates in percent
local_topologies = topology_data[:, nearest_idx]
local_angles = angle_data[:, nearest_idx]


# plot the local occultation rate
y_min = -1
y_max = 101
fig, ax = plt.subplots(dpi=200)
ax.set_xlabel(f"time[minutes after {j2000_to_utc_datetime(times[0])}]")
ax.set_ylabel("occultation[%]", color='navy')
ax.tick_params(axis='y', labelcolor='navy')
ax.set_title(f"location: lon({lon_observ_deg}°), lat({lat_observ_deg}°) ({location_name})")
ax.set_xlim(0, minutes_after_start[-1])
ax.set_ylim(y_min, y_max)

# plot the day mask (i.e. make background grey when it is night)
for i in range(num_steps-1):
    color = 'white' if local_day_mask[i] else 'grey'
    ax.axvspan(minutes_after_start[i], minutes_after_start[i+1], y_min, y_max, color=color)

# plot the entire line in turquoise
ax.plot(minutes_after_start, local_occults, color='turquoise', label="no eclipse")
# overlay the partial segment in navy
ax.plot(minutes_after_start[local_topologies > 0], local_occults[local_topologies > 0], color='navy', label="partial")
# overlay the annular segment in orange
ax.plot(minutes_after_start[local_topologies == 2], local_occults[local_topologies == 2], color='orange', label="annular")
# overlay the total segment in crimson
ax.plot(minutes_after_start[local_topologies == 3], local_occults[local_topologies == 3], color='crimson', label="total")

# add a second y-axis on the right-hand side
ax2 = ax.twinx()
ax2.set_ylabel('clockwise moon angle from north[deg]', color='mediumseagreen')
ax2.tick_params(axis='y', labelcolor='mediumseagreen')
ax2.set_ylim(0, 360)

# plote the binned angle data
ax2.fill_between(minutes_after_start, local_angles, local_angles+6, alpha=0.7, color='mediumseagreen')


# finish the plot
ax.grid()
ax.legend(loc=('upper left' if local_angles[0] < local_angles[-1] else 'upper right'))
plt.tight_layout()
plt.show()
plt.close(fig)


# apply day mask
local_occults = local_occults[local_day_mask]
local_topologies = local_topologies[local_day_mask]
times = times[local_day_mask]

print(f"Chosen location: lon({lon_observ_deg}°), lat({lat_observ_deg}°) ({location_name})")

if (local_topologies > 0).any():
    overall_start_time = times[local_topologies > 0][0]
    overall_end_time = times[local_topologies > 0][-1]
else:
    print("At this location the eclipse is not visible at all.")
    exit(0)

if (local_topologies == 3).any(): # the eclipse undergoes totality at the chosen location
    totality_start_time = times[local_topologies == 3][0]
    totality_end_time = times[local_topologies == 3][-1]

    print("At this location the eclipse undergoes totality.")
    print(f"It is visible there overall from {j2000_to_utc_datetime(overall_start_time)} to {j2000_to_utc_datetime(overall_end_time)}.")
    print(f"Totality lasts from {j2000_to_utc_datetime(totality_start_time)} to {j2000_to_utc_datetime(totality_end_time)}.")

elif (local_topologies == 2).any(): # the eclipse undergoes annularity at the chosen location
    annularity_start_time = times[local_topologies == 2][0]
    annularity_end_time = times[local_topologies == 2][-1]

    annular_occult = np.mean(local_occults[local_topologies == 2])

    print("At this location the eclipse undergoes annularity.")
    print(f"It is visible there overall from {j2000_to_utc_datetime(overall_start_time)} to {j2000_to_utc_datetime(overall_end_time)}.")
    print(f"Annularity lasts from {j2000_to_utc_datetime(annularity_start_time)} to {j2000_to_utc_datetime(annularity_end_time)}.")
    print(f"The occultation during annularity is {annular_occult}%.")

else: # the eclipse is only partial at the chosen location
    occult_argmax = np.argmax(local_occults)
    max_occult = local_occults[occult_argmax]
    max_time = times[occult_argmax]

    print("At this location the eclipse is only partial.")
    print(f"It is visible there overall from {j2000_to_utc_datetime(overall_start_time)} to {j2000_to_utc_datetime(overall_end_time)}.")
    print(f"The maximum occultation is {max_occult}% at {j2000_to_utc_datetime(max_time)}.")




