import numpy as np
import matplotlib.pyplot as plt
import sunpy.map
import astropy.units as u
import pandas as pd
from astropy.io import fits
from sunpy.coordinates import frames

npz_file = "sdo_dop_eclipse.npz"

# Directory with your FITS files
hmi_fits = "hmi.ic_45s.20231014_152700_TAI.2.continuum.fits"

# Load NPZ data
data = np.load(npz_file)
print("Arrays in file:", list(data.keys()))

# Arrays available: raw, satellite, rotation, mer_flows, corrected
fields = ["raw", "satellite", "rotation", "mer_flows", "corrected"]

# Load reference HMI map for WCS
hmi_map = sunpy.map.Map(hmi_fits)

# Ensure the HMI map is 4096×4096 to match NPZ
if hmi_map.data.shape != (4096, 4096):
    hmi_map = hmi_map.resample((4096, 4096) * u.pixel)

def plot_velocity_field(name, arr):

    # Create a new SunPy map with HMI WCS but using NPZ velocity data
    vmap = sunpy.map.Map(arr, hmi_map.meta.copy())

    # Set symmetric velocity limits for blue-red colormap
    vmax = np.nanmax(arr)
    vmin = np.nanmin(arr)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection=vmap)

    im = vmap.plot(
        axes=ax,
        cmap="RdBu_r",     # blue = negative, red = positive
        vmin=vmin,
        vmax=vmax
    )

    plt.title(f"{name} Velocity Field (HPLN/HPLT)")
    plt.colorbar(im, ax=ax, label="Velocity (m/s)")
    plt.savefig(f"{name}_velocity.png", dpi=150)
    plt.show()

# for f in fields:
#     print(f"Plotting {f}...")
#     plot_velocity_field(f, data[f])

# plot_velocity_field("rotation_satellite_corrected", data["raw"] - data["rotation"] - data["satellite"])

hmi_mapl = sunpy.map.Map(data["raw"] - data["rotation"] - data["satellite"], hmi_map.meta.copy())
# hmi_map = vmap.resample(((1024, 1024)) * u.pixel) #(198, 396)
nx, ny = hmi_map.data.shape
x_pix, y_pix = np.meshgrid(np.arange(nx), np.arange(ny))

# Convert pixel → world → HGS lat/lon
world = hmi_map.pixel_to_world(x_pix * u.pix, y_pix * u.pix).heliographic_stonyhurst
hgs = world.transform_to(frames.HeliographicCarrington(obstime="2023-10-14T16:35:39.500000"))

lon = hgs.lon.deg
lat = hgs.lat.deg

# Build DataFrame
df = pd.DataFrame({
    "lon_deg": lon.ravel(),
    "lat_deg": lat.ravel(),
    "velocity": hmi_map.data.ravel()
})
df = df.dropna()
# Save CSV
csv_name = f"velocity_field.csv"
df.to_csv(csv_name, index=False)

# # Set symmetric velocity limits for blue-red colormap
# arr = data["raw"] - data["rotation"] - data["satellite"]
# vmax = np.nanmax(arr)
# vmin = np.nanmin(arr)

# fig = plt.figure(figsize=(8, 8))
# ax = fig.add_subplot(111, projection=hmi_map)

# im = hmi_map.plot(
#     axes=ax,
#     cmap="RdBu_r",     # blue = negative, red = positive
#     vmin=vmin,
#     vmax=vmax
#     )
# name = "resampled"
# plt.title(f"{name} Velocity Field (HPLN/HPLT)")
# plt.colorbar(im, ax=ax, label="Velocity (m/s)")
# plt.savefig(f"{name}_velocity.png", dpi=150)
# plt.show()
