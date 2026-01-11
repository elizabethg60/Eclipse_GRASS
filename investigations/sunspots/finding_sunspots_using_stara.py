import glob
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import sunpy.map
from sunkit_image.stara import stara
from skimage.measure import label, regionprops_table
import matplotlib.pyplot as plt
import pandas as pd
from sunpy.coordinates import frames
from astropy.coordinates import EarthLocation, AltAz
import numpy as np
from skimage.morphology import disk, binary_dilation
from astropy.constants import R_sun

# Directory with your FITS files
fits_dir = "/storage/home/efg5335/sunpy/data/"
fits_files = sorted(glob.glob(f"{fits_dir}/*.fits"))

# --- Read all FITS headers and observation times ---
fits_times = []
for f in fits_files:
    with fits.open(f) as hdul:
        hdr = hdul[1].header  # Assuming the data is in extension 1
        fits_times.append(Time(hdr["DATE-OBS"]))

# Target timestamps (your NEID October times)
neid_october = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
target_times = Time(neid_october)

# --- Main loop over target times (example: first target only) ---
for ind, t in enumerate([target_times[50]]):
    # Find closest FITS file
    idx = min(range(len(fits_times)), key=lambda i: abs(fits_times[i]-t))
    f = fits_files[idx]

    # Load FITS into SunPy map
    hmi_map = sunpy.map.Map(f).resample((1024, 1024) * u.pixel).rotate(order=3)
    # print(hmi_map.scale) #2.01616335 arcsec / pix

    # Run STARA sunspot detection
    stara_segments = stara(hmi_map, limb_filter=10 * u.percent)

    # Label and extract region properties
    labelled = label(stara_segments)
    regions = regionprops_table(
        labelled, 
        hmi_map.data,
        properties=["label", "weighted_centroid", "area", "min_intensity", "equivalent_diameter"]
    )
    # Add observation time
    regions["obstime"] = [hmi_map.date] * len(regions["label"])
    # Convert to pandas DataFrame
    df = pd.DataFrame(regions)
    df = df[df["area"] >= 50]

    # Extract centroid pixel coordinates
    centroids_x = df["weighted_centroid-1"].values  # x / column
    centroids_y = df["weighted_centroid-0"].values  # y / row
    equivalent_diameter = df["equivalent_diameter"].values #pixels
    equivalent_diameter_km = ((equivalent_diameter * 2.01616335 * u.arcsec).to(u.rad) * 149597870).value
    equivalent_diameter_arcsec = equivalent_diameter * 2.01616335

    # Convert to world coordinates (Heliographic Stonyhurst)
    centroids_world = hmi_map.pixel_to_world(centroids_x * u.pix, centroids_y * u.pix).heliographic_stonyhurst
    # Convert HGS -> Heliographic Carrington (closest to IAU_SUN / body-fixed)
    centroids_world_hgc = centroids_world.transform_to(frames.HeliographicCarrington(obstime=hmi_map.date))

    contrasts = []
    latitude = []
    longitude = []
    for _, row in df.iterrows():
        lab = int(row["label"])
        spot_mask = (labelled == lab)

        if not np.any(spot_mask):
            continue  # skip if no pixels left for this label

        area_pixels = float(row["area"])

        # Spot intensity
        spot_vals = hmi_map.data[spot_mask]
        spot_mean = np.nanmean(spot_vals)

        # Background annulus
        r_pix = np.sqrt(area_pixels / np.pi)
        selem_out = disk(int(r_pix * 2) + 2)
        selem_in = disk(int(r_pix) + 1)
        dil_out = binary_dilation(spot_mask, selem_out)
        dil_in = binary_dilation(spot_mask, selem_in)
        annulus_mask = dil_out & (~dil_in) & (labelled == 0)
        bg_vals = hmi_map.data[annulus_mask]
        bg_median = np.nanmedian(bg_vals)

        # Contrast
        contrast = (bg_median - spot_mean) / bg_median
        contrasts.append(contrast)

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=hmi_map)
    hmi_map.plot(axes=ax)
    ax.contour(stara_segments, levels=0)
    ax.scatter(centroids_x, centroids_y, color="red", s=30)
    plt.savefig("image.png")

    df_hgs = pd.DataFrame({
        "lon": [c.lon.deg for c in centroids_world_hgc],
        "lat": [c.lat.deg for c in centroids_world_hgc],
        "diameter_km": equivalent_diameter_km,
        "diameter_arcsec": equivalent_diameter_arcsec,
        "contrast": contrasts
    })

    # Optional: save to CSV
    df_hgs.to_csv("sunspots.csv", index=False)
