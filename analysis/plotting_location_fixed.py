import os
import numpy as np
import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
import xarray as xr
import PIL
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from tqdm import tqdm
import imageio

from utils import *


# specify location to be kept fixed
center_lon_deg = 0
center_lat_deg = -90

file_path = "/data/eclipse_data/28.05.2026_18-36-25/solar_eclipse_of_05.02.2000.nc"

anim_dir = os.path.join("/home/aurora/eclipse_data", file_path.split('/')[-2])
if not os.path.exists(anim_dir): # create the directory if it doesn't exist yet
    os.mkdir(anim_dir)
anim_path = os.path.join(anim_dir, file_path.split('/')[-1][:-3] + f"_location_fixed({center_lon_deg}, {center_lat_deg}).mp4")
version = 0
while os.path.exists(anim_path): # add a version number if the file already exists
    version += 1
    anim_path = anim_path[:-4] + f"({version})" + ".mp4"


dataset = xr.open_dataset(file_path, chunks={"steps": 1, "grid": -1})

times = dataset["time"].values
sun_lon = dataset["lon_sun"].values
sun_lat = dataset["lat_sun"].values
num_steps = times.shape[0]

lon_grid = dataset["lon_grid"].values
lat_grid = dataset["lat_grid"].values
grid_size = lat_grid.shape[0]

occultation_data = dataset["occultation_data"].data

dataset.close()

center_lon_rad = center_lon_deg * np.pi/180
center_lat_rad = center_lat_deg * np.pi/180

lon_mesh, lat_mesh = np.meshgrid(np.linspace(-np.pi, np.pi, 2000), np.linspace(np.pi/2, -np.pi/2, 1000))

# setup type "nearest" interpolation
nearest_idcs = get_nearest_idcs(lon_grid, lat_grid,
                                lon_mesh, lat_mesh, k=1)

# load background image of earth
earth_image_path = "/home/aurora/Desktop/NE1_50M_SR_W.tif"
earth_image = PIL.Image.open(earth_image_path)


def draw_shadow(idx, ax, shadow_plot, under_sun_point_plot):

    sun_lon_rad = sun_lon[idx]
    sun_lat_rad = sun_lat[idx]

    # compute the relative intensity of the sun
    sun_cos_zen = np.cos(lat_mesh)*np.cos(sun_lat_rad)*np.cos(lon_mesh-sun_lon_rad) + np.sin(lat_mesh)*np.sin(sun_lat_rad)
    sun_rel_intensity = np.maximum(0, sun_cos_zen)**0.7

    # interpolate the occultation due to the eclipse
    occults = occultation_data[idx].compute()
    occults_interp = occults[nearest_idcs]

    # plot the combined shadow (night and eclipse)
    shadow_image = np.zeros(shape=(*lat_mesh.shape, 4))
    shadow_image[:,:,3] = 1 - sun_rel_intensity * (1-occults_interp)
    shadow_plot.remove()
    shadow_plot = ax.imshow(shadow_image, origin="upper", transform=ccrs.PlateCarree(), extent=[-180,180,-90,90], zorder=2)

    # plot the under-sun point
    under_sun_point_plot.remove()
    under_sun_point_plot = ax.scatter(sun_lon_rad*180/np.pi, sun_lat_rad*180/np.pi, color="orange", marker="*", transform=ccrs.PlateCarree(), zorder=3)

    # for debugging: print the original Fibonacci sphere
    # ax.scatter(lon_grid*180/np.pi, lat_grid*180/np.pi, color="red", marker=".", s=1, transform=ccrs.PlateCarree())

    # add the time as title
    ax.set_title(
        j2000_to_utc_datetime(times[idx]),
        color='white',
        fontsize=14,
        pad=20,
        fontname='sans-serif'
    )

    return shadow_plot, under_sun_point_plot



# initialize an empty figure
fig = plt.figure(figsize=(10, 10), dpi=112)
plt.suptitle(os.path.split(file_path)[1][:-3].replace('_', ' ')
            + f", centered at lon({center_lon_deg})°, lat({center_lat_deg})° ({get_location_name(center_lon_rad, center_lat_rad)}).")

# make the background black
fig.patch.set_facecolor('black')

# set up the Orthographic (spherical globe) projection
projection = ccrs.Orthographic(central_longitude=center_lon_deg, central_latitude=center_lat_deg)

# create the axes and explicitly strip out the default rectangular box border
ax = plt.axes(projection=projection, frameon=False)
ax.set_global()

# plot earth
ax.stock_img()
#ax.imshow(earth_image, origin="upper", transform=ccrs.PlateCarree(), extent=[-180,180,-90,90])

# add country borders and coastlines
ax.add_feature(
    cfeature.COASTLINE.with_scale('50m'),
    edgecolor='white',
    linewidth=0.6,
    alpha=0.75
)
ax.add_feature(
    cfeature.BORDERS.with_scale('50m'),
    edgecolor='white',
    linewidth=0.4,
    alpha=0.9,
    linestyle='-'
)
# add subtle geographic gridlines
gridlines = ax.gridlines(
    crs=ccrs.PlateCarree(),
    draw_labels=False,
    linewidth=0.5,
    color='white',
    alpha=0.5,
    linestyle='--'
)
gridlines.xlocator = plt.MultipleLocator(15)
gridlines.ylocator = plt.MultipleLocator(15)

# plot the central point
ax.scatter(0, 0, color="crimson", marker="x", zorder=3)

# initialize the actual animated layers: the shadow and the under-sun point
shadow_plot = ax.imshow(np.zeros(shape=(*lat_mesh.shape, 4)), origin="upper", transform=ccrs.PlateCarree(), extent=[-180,180,-90,90], zorder=2)
under_sun_point_plot = ax.scatter([], [], color="orange", marker="*", transform=ccrs.PlateCarree(), zorder=3)


# use imageio to open a secure video container using its built-in ffmpeg binary
print("Initializing Video Container with ImageIO...")
video_writer = imageio.get_writer(anim_path, fps=5, codec='libx264', quality=8, pixelformat='yuv420p')

# generate all the frames
for idx in tqdm(range(num_steps), desc="Rendering"):
    shadow_plot, under_sun_point_plot = draw_shadow(idx, ax, shadow_plot, under_sun_point_plot)

    # draw the canvas, extract raw pixels as an RGBA image matrix, and write directly to video
    fig.canvas.draw()
    frame_pixels = np.asarray(fig.canvas.buffer_rgba())
    # convert RGBA matrix to an RGB image (dropping the 4th alpha layer for mp4 video storage compatibility)
    rgb_frame = frame_pixels[:, :, :3]

    # pass the frame to the video writer
    video_writer.append_data(rgb_frame)


video_writer.close()
plt.close(fig)

print(f"\n[SUCCESS] Render complete! Animation saved to: {anim_path}")



