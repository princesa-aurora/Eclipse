import os
import numpy as np
import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
import xarray as xr
import dask
import PIL
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from tqdm import tqdm
import multiprocessing as mp
import subprocess

from utils import *


file_path = "/data/eclipse_data/01.06.2026_19-22-00/solar_eclipse_of_05.02.2000.nc"

num_procs = 4 #mp.cpu_count()

anim_dir = os.path.join("/home/aurora/eclipse_data", file_path.split('/')[-2])
anim_path = os.path.join(anim_dir, file_path.split('/')[-1][:-3] + "_sun_fixed.mp4")
frame_dir = os.path.join(os.path.split(file_path)[0], "frames_sun_fixed")


dask.config.set(scheduler='single-threaded')
dataset = xr.open_dataset(file_path, chunks={"steps": 1, "grid": -1})

times = dataset["time"].values
sun_lon = dataset["lon_sun"].values
sun_lat = dataset["lat_sun"].values
num_steps = times.shape[0]
dt = (times[-1]-times[0])/(num_steps-1)

lon_grid = dataset["lon_grid"].values
lat_grid = dataset["lat_grid"].values
grid_size = lat_grid.shape[0]

occultation_data = dataset["occultation_data"].data

dataset.close()


lon_mesh, lat_mesh = np.meshgrid(np.linspace(-np.pi, np.pi, 2000), np.linspace(np.pi/2, -np.pi/2, 1000))

# setup type "nearest" interpolation
nearest_idcs = get_nearest_idcs(lon_grid, lat_grid,
                                lon_mesh, lat_mesh, k=1)

# load background image of earth
earth_image_path = "/home/aurora/Desktop/NE1_50M_SR_W.tif"
earth_image = PIL.Image.open(earth_image_path)



# initialize an empty figure
fig = plt.figure(figsize=(10, 10), dpi=112)


def draw_frame(idx, fig):

    sun_lon_rad = sun_lon[idx]
    sun_lat_rad = sun_lat[idx]

    sun_lon_deg = sun_lon_rad *180/np.pi
    sun_lat_deg = sun_lat_rad *180/np.pi

    # make the background black
    fig.patch.set_facecolor('black')

    # set up the Orthographic (spherical globe) projection
    projection = ccrs.Orthographic(central_longitude=sun_lon_deg, central_latitude=sun_lat_deg)

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

    # compute the relative intensity of the sun
    sun_cos_zen = np.cos(lat_mesh)*np.cos(sun_lat_rad)*np.cos(lon_mesh-sun_lon_rad) + np.sin(lat_mesh)*np.sin(sun_lat_rad)
    sun_rel_intensity = np.maximum(0, sun_cos_zen)**0.7

    # interpolate the occultation due to the eclipse
    occults = occultation_data[idx].compute()
    occults_interp = occults[nearest_idcs]

    # plot the combined shadow (night and eclipse)
    shadow_image = np.zeros(shape=(*lat_mesh.shape, 4))
    shadow_image[:,:,3] = 1 - sun_rel_intensity * (1-occults_interp)
    ax.imshow(shadow_image, origin="upper", transform=ccrs.PlateCarree(), extent=[-180,180,-90,90], zorder=2)

    # plot the under-sun point
    ax.scatter(0, 0, color="orange", marker="*")

    # for debugging: scatter plot the original Fibonacci sphere
    # ax.scatter(lon_grid*180/np.pi, lat_grid*180/np.pi, color="red", marker=".", s=1, transform=ccrs.PlateCarree())

    # add titles
    plt.suptitle(os.path.split(file_path)[1][:-3].replace('_', ' ') + ", centered below the sun")
    ax.set_title(
        j2000_to_utc_datetime(times[idx]),
        color='white',
        fontsize=14,
        pad=20,
        fontname='sans-serif'
    )


def render_frames(proc_params):
    # render and save the frames in idx_range
    proc_id, start_idx, stop_idx = proc_params

    # setup progress bar
    progress_bar = tqdm(total=stop_idx-start_idx,
                        desc=f"process #{proc_id}", 
                        position=proc_id,
                        leave=False)

    for i in range(start_idx, stop_idx):
        frame_path = os.path.join(frame_dir, f"frame_{i}.png")
        if os.path.exists(frame_path):
            progress_bar.update(1)
            continue # if the frame exists already there's no need to render it again

        # clear the frame
        fig.clear()

        # draw the new frame
        draw_frame(i, fig)

        # save the frame
        fig.savefig(frame_path, bbox_inches='tight')

        # update progress bar
        progress_bar.update(1)



if __name__ == "__main__":
    # create the directories if they don't exist yet
    if not os.path.exists(anim_dir):
        os.mkdir(anim_dir)
    if not os.path.exists(frame_dir):
        os.mkdir(frame_dir)

    # distribute frames to processes
    tasks = [(i, int(num_steps/num_procs*i), int(num_steps/num_procs*(i+1))) for i in range(num_procs)]

    print(f"Rendering across {num_procs} CPU processes...")

    # spawn the processes
    with mp.Pool(processes=num_procs) as pool:
        pool.map(render_frames, tasks)

    print("All frames successfully rendered. Starting FFmpeg video compilation...")

    # stitch the pngs together using ffmpeg in a separate terminal command
    ffmpeg_cmd = [
        "ffmpeg", "-y",                    # overwrite output file if it exists
        "-framerate", str(120/dt),         # framerate [fps] (sped up by 120 times compared to reality)
        "-i", f"{frame_dir}/frame_%d.png", # input directory
        "-c:v", "libx264",                 # industry standard h.264 video codec
        "-pix_fmt", "yuv420p",             # colorspace compatible with all video players
        "-crf", "18",                      # high quality, low compression factor
        anim_path                          # file path for the animation 
    ]
    subprocess.run(ffmpeg_cmd, check=True)

    print(f"Animation completed successfully, saved to {anim_path}.")


plt.close(fig)


