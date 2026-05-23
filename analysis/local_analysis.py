import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import PIL

from utils import *


file_path = "/home/aurora/eclipse_data/22.05.2026_20:12:50/solar_eclipse_of_5.2.2000.nc"

lon_observ_deg = 0
lat_observ_deg = -90


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


lon_observ = lon_observ_deg * np.pi/180
lat_observ = lat_observ_deg * np.pi/180

# setup type "nearest" interpolation
nearest_idcs = get_nearest_idcs(lon_grid, lat_grid,
                                       lon_observ, lat_observ, k=3)

# pre-compute the interpolation weights
interp_weights = get_interp_weights(lon_grid[nearest_idcs], lat_grid[nearest_idcs],
                                           lon_observ, lat_observ)


# interpolate the occultation at all times
local_occults = []
for idx in range(num_steps):
    sun_lon_rad = sun_lon[idx]
    sun_lat_rad = sun_lat[idx]

    # compute the relative intensity of the sun
    #sun_cos_zen = np.cos(lat_observ)*np.cos(sun_lat_rad)*np.cos(lon_observ-sun_lon_rad) + np.sin(lat_observ)*np.sin(sun_lat_rad)
    #sun_rel_intensity = np.maximum(0, sun_cos_zen)**0.7

    # only need the occultations at the nearest points for interpolation
    occults = occultation_data[idx][nearest_idcs]
    # interpolate the occultation rate at the observation point
    interp_occult = interp_weights.dot(occults)

    local_occults.append(interp_occult)


local_occults = np.array(local_occults)
# convert from fraction to percent
local_occults *= 100

minutes_after_start = (times-times[0])/60


# plot the local occultation rate
fig = plt.figure(dpi=200)
plt.grid()
plt.xlabel(f"time[minutes after {j2000_to_utc_datetime(times[0])}]")
plt.ylabel("occultation[%]")
plt.title(f"Location: lon({lon_observ_deg}°), lat({lat_observ_deg}°)")

plt.plot(minutes_after_start, local_occults, color="navy")
plt.grid()
plt.show()
plt.close(fig)



print(f"Chosen location: lon({lon_observ_deg}°), lat({lat_observ_deg}°)")


# compute the times of the eclipse: start and stop of partiality and totality (if any)
if (local_occults > 0).any():
    first_idx = np.argmax(local_occults > 0)
    last_idx = len(local_occults)-1 - np.argmax(local_occults[::-1] > 0)

    t_start = interpolate_linearly(0,
                local_occults[first_idx], times[first_idx],
                local_occults[first_idx+1], times[first_idx+1])

    t_end = interpolate_linearly(0,
                local_occults[last_idx], times[last_idx],
                local_occults[last_idx-1], times[last_idx-1])

    print(f"The eclipse is visible between {j2000_to_utc_datetime(t_start)} and {j2000_to_utc_datetime(t_end)} at the chosen location.")

else:
    print("The eclipse cannot be seen at the chosen location at all.")


if ((local_occults) == 1).any():
    first_idx = np.argmax(local_occults == 1)
    last_idx = len(local_occults)-1 - np.argmax(local_occults[::-1] == 1)

    t_start = interpolate_linearly(1,
                local_occults[first_idx-1], times[first_idx-1],
                local_occults[first_idx-2], times[first_idx-2])

    t_end = interpolate_linearly(1,
                local_occults[last_idx+1], times[last_idx+1],
                local_occults[last_idx+2], times[last_idx+2])

    print(f"Totality lasts from {j2000_to_utc_datetime(t_start)} to {j2000_to_utc_datetime(t_end)} at the chosen location.")

else:
    print("The eclipse isn't total during any time at the chosen loaction.")




