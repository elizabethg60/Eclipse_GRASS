import numpy as np
from sunpy.sun import constants
from sunpy.time import parse_time
from sunpy.coordinates.frames import HeliographicStonyhurst
import astropy.units as u
from astropy.coordinates import AltAz, HADec, SkyCoord, Angle
from astropy.coordinates.representation import CartesianRepresentation, SphericalRepresentation
from astropy.modeling.rotations import SphericalRotationSequence
from sunpy.coordinates import spice
from sunpy.data import cache
import spiceypy as spiceypy
import math
from datetime import datetime
from sunpy.coordinates.spice import get_rotation_matrix
from scipy.spatial.transform import Rotation as R
import pdb

def rotation_matrix_from_vectors(vec1, vec2):
    """
    Calculates the rotation matrix that aligns vec1 to vec2.

    Args:
        vec1 (numpy.ndarray): The initial vector.
        vec2 (numpy.ndarray): The target vector.

    Returns:
        numpy.ndarray: The rotation matrix.
    """
    a = vec1 / np.linalg.norm(vec1)
    b = vec2 / np.linalg.norm(vec2)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2)) if s != 0 else np.eye(3)
    return rotation_matrix

from scipy.spatial.transform import Rotation as R

def rotation_matrix_from_state_vectors(v1, v2):
    """
    Compute the 6x6 rotation matrix to align vector v1 to v2 in 6D space.
    
    Parameters:
    v1 (numpy.ndarray): The first vector of length 6.
    v2 (numpy.ndarray): The second vector of length 6.
    
    Returns:
    numpy.ndarray: The 6x6 rotation matrix that aligns v1 to v2.
    """
    
    # Ensure vectors are numpy arrays
    v1 = np.array(v1)
    v2 = np.array(v2)

    # Compute the cross-covariance matrix between v1 and v2
    A = np.outer(v1, v2)

    # Perform Singular Value Decomposition (SVD)
    U, _, Vt = np.linalg.svd(A)

    # Rotation matrix is given by U @ Vt
    R = U @ Vt
    
    return R


kernel_urls = [
    "lsk/naif0012.tls",
    "spk/planets/de440.bsp",
    "spk/satellites/jup365.bsp",
    "pck/moon_pa_de440_200625.bpc",
    "pck/earth_latest_high_prec.bpc",
    "pck/earth_1962_240827_2124_combined.bpc",
    "fk/planets/earth_assoc_itrf93.tf",
    "pck/pck00010.tpc"
]
kernel_urls = [f"https://naif.jpl.nasa.gov/pub/naif/generic_kernels/{url}"
               for url in kernel_urls]

kernel_files = [cache.download(url) for url in kernel_urls]
spice.initialize(kernel_files)
spice.install_frame('IAU_Earth')
spice.install_frame('IAU_Sun')
spice.install_frame('J2000')

def solar_orientation(location, time="now", coordinates="altaz"):
    """
    The orientation of the sun as viewed from a particular location on earth
    with an image plane relative to the given coordinate system.

    This function is inspired by sunpy.coordinates.sun.orientation, whose
    output should be the same as `phi`.

    Parameters
    ----------
    location : `~astropy.coordinates.EarthLocation`
        Observer location on Earth.
    time : time_like
        Time in a sunpy.time.parse_time-compatible format.
    coordinates : str, optional
        Coordinate system in which the sun is viewed. The default is "altaz".

    Returns
    -------
    phi : `~astropy.coordinates.Angle`
        The solar rotation in the image plane.
    theta : `~astropy.coordinates.Angle`
        The solar rotation towards or away from the image plane.

    """
    # Define the frame where its Z axis is aligned with local zenith
    obstime = parse_time(time)
    if coordinates=="altaz":
        frame = AltAz(obstime=obstime, location=location)
    elif coordinates=="equat":
        frame = HADec(obstime=obstime, location=location)

    # Find the Sun's center in HGS at the frame's observation time(s)
    sun_center_repr = SphericalRepresentation(0*u.deg, 0*u.deg, 0*u.km)
    sun_center = SkyCoord(sun_center_repr._apply('repeat', frame.obstime.size),
                          frame=HeliographicStonyhurst, obstime=frame.obstime)
    sun_center = sun_center.transform_to(frame)

    # Find the Sun's north pole in HGS at the frame's observation time(s)
    sun_pole_repr = SphericalRepresentation(0*u.deg, 90*u.deg, constants.radius)
    sun_pole = SkyCoord(sun_pole_repr._apply('repeat', frame.obstime.size),
                         frame=HeliographicStonyhurst, obstime=frame.obstime)
    sun_pole = sun_pole.transform_to(frame)

    # Find solar north in the frame's coordinate system
    solar_north = sun_pole.data.to_cartesian() - sun_center.data.to_cartesian()
    solar_north /= solar_north.norm()

    # Find the image axes in the frame's coordinate system
    # The image plane is first defined relative to the view along the x axis of
    # the coordinate system (0/0 lon/lat) and then rotated to point at the Sun
    x_repr = CartesianRepresentation(0,-1,0) # x axis in the image plane
    y_repr = CartesianRepresentation(0,0,1) # y axis in the image plane
    z_repr = CartesianRepresentation(-1,0,0) # z axis in the image plane
    x = SkyCoord(x_repr._apply('repeat', frame.obstime.size), frame=frame)
    y = SkyCoord(y_repr._apply('repeat', frame.obstime.size), frame=frame)
    z = SkyCoord(z_repr._apply('repeat', frame.obstime.size), frame=frame)

    # Rotate such that the image plane z axis points towards the Sun
    im_x = []
    im_y = []
    im_z = []
    for (sci, xi, yi, zi) in zip(sun_center, x, y, z):
        rot = SphericalRotationSequence( # first rotate alt, then rotate az
            (sci.data.lat, -sci.data.lon), "yz") # (1,0,0) is the origin (0,0) of spherical representation

        lon_x, lat_x = rot(xi.data.lon, xi.data.lat) * u.deg
        im_x.append(SphericalRepresentation(lon_x, lat_x, 1).to_cartesian().xyz)

        lon_y, lat_y = rot(yi.data.lon, yi.data.lat) * u.deg
        im_y.append(SphericalRepresentation(lon_y, lat_y, 1).to_cartesian().xyz)

        lon_z, lat_z = rot(zi.data.lon, zi.data.lat) * u.deg
        im_z.append(SphericalRepresentation(lon_z, lat_z, 1).to_cartesian().xyz)
    im_x = CartesianRepresentation(*zip(*im_x))
    im_y = CartesianRepresentation(*zip(*im_y))
    im_z = CartesianRepresentation(*zip(*im_z))

    # Calculate Angles
    phi = np.arctan2(solar_north.dot(im_y), solar_north.dot(im_x)).to('deg')
    theta = np.arccos(solar_north.dot(im_z)).to('deg')

    # Return scalar if there is only one time
    if len(phi) == 1:
        phi = phi[0]
        theta = theta[0]
    
    #AltAz Frame
    sun_coord = SkyCoord(cartesian_rep, frame="spice_IAU_SUN", obstime=frame.obstime)
    sun_coord_frame = (sun_coord.transform_to(frame)).data.to_cartesian()
    AltAz_angle = np.arctan2(sun_coord_frame.dot(im_y), sun_coord_frame.dot(im_x)).to('deg')


    zenith = SkyCoord(lon=obs_long*u.deg, lat=obs_lat*u.deg, distance=alt*u.km, d_lon=0*u.arcmin/u.week, d_lat=0*u.arcmin/u.d, d_distance=0*u.km/u.min,
              frame='spice_IAU_EARTH', obstime=frame.obstime)
    zenith_transformed = zenith.transform_to("spice_J2000")
    #IAU_Earth Frame Spice
    sun_coord_frame = np.dot(spiceypy.pxform("IAU_SUN", "IAU_Earth", spiceypy.utc2et(time.strftime("%Y-%m-%dT%H:%M:%S"))), sun_coord.cartesian.xyz.value)
    IAU_Earth_angle_spice = sun_coord_frame - (np.dot(sun_coord_frame, zenith.cartesian.xyz.value) / (np.linalg.norm(zenith.cartesian.xyz.value)**2))*zenith.cartesian.xyz.value
    # IAU_Earth_angle_spice = np.arccos(np.dot(sun_coord_frame, zenith.cartesian.xyz.value) / (np.linalg.norm(sun_coord_frame) * np.linalg.norm(zenith.cartesian.xyz.value)))
    #IAU_Earth Frame Sunpy
    # pdb.set_trace()
    sun_coord_frame = (sun_coord.transform_to("spice_IAU_EARTH")).cartesian.xyz.value
    IAU_Earth_angle_sunpy = sun_coord_frame - (np.dot(sun_coord_frame, zenith.cartesian.xyz.value) / (np.linalg.norm(zenith.cartesian.xyz.value)**2))*zenith.cartesian.xyz.value
    # IAU_Earth_angle_sunpy = np.arccos(np.dot(sun_coord_frame, zenith.cartesian.xyz.value) / (np.linalg.norm(sun_coord_frame) * np.linalg.norm(zenith.cartesian.xyz.value)))
    
    # rot_mat = (rotation_matrix_from_vectors(sun_coord.cartesian.xyz.value, (sun_coord.transform_to("spice_J2000")).cartesian.xyz.value))
    pos_vector = list(zenith.cartesian.xyz.value)
    vel_vector = list(np.atleast_1d(zenith.velocity._values)[0])
    pos_vector_transformed = list(zenith_transformed.cartesian.xyz.value)
    vel_vector_transformed = list(np.atleast_1d(zenith_transformed.velocity._values)[0])
    rot_mat = rotation_matrix_from_state_vectors(np.array(pos_vector + vel_vector), np.array(pos_vector_transformed+vel_vector_transformed))
    
    # Return angles relative to y axis in image plane
    return Angle(phi) - 90*u.deg, Angle(theta) - 90*u.deg, Angle(AltAz_angle) - 90*u.deg, np.degrees(IAU_Earth_angle_spice), np.degrees(IAU_Earth_angle_sunpy), rot_mat

#-------------------------------------------------------------

from astropy.coordinates import EarthLocation
import astropy.units as u
obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938

earth_radius = spiceypy.bodvrd("EARTH", "RADII", 3)[1][0]	
earth_radius_pole = spiceypy.bodvrd("EARTH", "RADII", 3)[1][2]	
sun_radius = spiceypy.bodvrd("SUN","RADII", 3)[1][1]

sun_radius = 1.0 * u.R_sun
sun_axis_sun = [0.0, 0.0, sun_radius.to(u.AU).value]
flat_coeff = (earth_radius - earth_radius_pole) / earth_radius
zenith_coord = (spiceypy.georec(math.radians(obs_long), math.radians(obs_lat), alt, earth_radius, flat_coeff))
cartesian_rep = CartesianRepresentation(x=sun_axis_sun[0] * u.AU, y=sun_axis_sun[1] * u.AU, z=sun_axis_sun[2] * u.AU)
location = EarthLocation.from_geodetic(obs_long, obs_lat, alt)


neid_october = ["2023-10-14T15:26:45.500000", "2023-10-14T15:28:07.500000", "2023-10-14T15:29:30.500000", "2023-10-14T15:30:53.500000", "2023-10-14T15:32:15.500000", "2023-10-14T15:33:38.500000", "2023-10-14T15:35:01.500000", "2023-10-14T15:36:23.500000", "2023-10-14T15:37:46.500000", "2023-10-14T15:39:09.500000", "2023-10-14T15:40:31.500000", "2023-10-14T15:41:54.500000", "2023-10-14T15:43:17.500000", "2023-10-14T15:44:39.500000", "2023-10-14T15:46:02.500000", "2023-10-14T15:47:25.500000", "2023-10-14T15:48:47.500000", "2023-10-14T15:50:10.500000", "2023-10-14T15:51:33.500000", "2023-10-14T15:52:56.500000", "2023-10-14T15:54:18.500000", "2023-10-14T15:55:41.500000", "2023-10-14T15:57:04.500000", "2023-10-14T15:58:26.500000", "2023-10-14T15:59:49.500000", "2023-10-14T16:01:12.500000", "2023-10-14T16:02:34.500000", "2023-10-14T16:03:57.500000", "2023-10-14T16:05:20.500000", "2023-10-14T16:06:42.500000", "2023-10-14T16:08:05.500000", "2023-10-14T16:09:28.500000", "2023-10-14T16:10:50.500000", "2023-10-14T16:12:13.500000", "2023-10-14T16:13:36.500000", "2023-10-14T16:14:58.500000", "2023-10-14T16:16:21.500000", "2023-10-14T16:17:44.500000", "2023-10-14T16:19:06.500000", "2023-10-14T16:20:29.500000", "2023-10-14T16:21:52.500000", "2023-10-14T16:23:15.500000", "2023-10-14T16:24:37.500000", "2023-10-14T16:26:00.500000", "2023-10-14T16:27:23.500000", "2023-10-14T16:28:45.500000", "2023-10-14T16:30:08.500000", "2023-10-14T16:31:31.500000", "2023-10-14T16:32:53.500000", "2023-10-14T16:34:16.500000", "2023-10-14T16:35:39.500000", "2023-10-14T16:37:01.500000", "2023-10-14T16:38:24.500000", "2023-10-14T16:39:47.500000", "2023-10-14T16:41:09.500000", "2023-10-14T16:42:32.500000", "2023-10-14T16:43:55.500000", "2023-10-14T16:45:17.500000", "2023-10-14T16:46:40.500000", "2023-10-14T16:48:03.500000", "2023-10-14T16:49:25.500000", "2023-10-14T16:50:48.500000", "2023-10-14T16:52:11.500000", "2023-10-14T16:53:33.500000", "2023-10-14T16:54:56.500000", "2023-10-14T16:56:19.500000", "2023-10-14T16:57:42.500000", "2023-10-14T16:59:04.500000", "2023-10-14T17:00:27.500000", "2023-10-14T17:01:50.500000", "2023-10-14T17:03:12.500000", "2023-10-14T17:04:35.500000", "2023-10-14T17:05:58.500000", "2023-10-14T17:07:20.500000", "2023-10-14T17:08:43.500000", "2023-10-14T17:10:06.500000", "2023-10-14T17:11:28.500000", "2023-10-14T17:12:51.500000", "2023-10-14T17:14:14.500000", "2023-10-14T17:15:36.500000", "2023-10-14T17:16:59.500000", "2023-10-14T17:18:22.500000", "2023-10-14T17:19:44.500000", "2023-10-14T17:21:07.500000", "2023-10-14T17:22:30.500000", "2023-10-14T17:23:52.500000", "2023-10-14T17:25:15.500000", "2023-10-14T17:26:38.500000", "2023-10-14T17:28:01.500000", "2023-10-14T17:29:23.500000", "2023-10-14T17:30:46.500000", "2023-10-14T17:32:09.500000", "2023-10-14T17:33:31.500000", "2023-10-14T17:34:54.500000", "2023-10-14T17:36:17.500000", "2023-10-14T17:37:39.500000", "2023-10-14T17:39:02.500000", "2023-10-14T17:40:25.500000", "2023-10-14T17:41:47.500000", "2023-10-14T17:43:10.500000", "2023-10-14T17:44:33.500000", "2023-10-14T17:45:55.500000", "2023-10-14T17:47:18.500000", "2023-10-14T17:48:41.500000", "2023-10-14T17:50:03.500000", "2023-10-14T17:51:26.500000", "2023-10-14T17:52:49.500000", "2023-10-14T17:54:11.500000", "2023-10-14T17:55:34.500000", "2023-10-14T17:56:57.500000", "2023-10-14T17:58:20.500000", "2023-10-14T17:59:42.500000", "2023-10-14T18:01:05.500000", "2023-10-14T18:02:28.500000", "2023-10-14T18:03:50.500000", "2023-10-14T18:05:13.500000", "2023-10-14T18:06:36.500000", "2023-10-14T18:07:58.500000", "2023-10-14T18:09:21.500000", "2023-10-14T18:10:44.500000", "2023-10-14T18:12:06.500000", "2023-10-14T18:13:29.500000", "2023-10-14T18:14:52.500000", "2023-10-14T18:16:14.500000", "2023-10-14T18:17:37.500000", "2023-10-14T18:19:00.500000", "2023-10-14T18:20:22.500000", "2023-10-14T18:21:45.500000", "2023-10-14T18:23:08.500000", "2023-10-14T18:24:30.500000", "2023-10-14T18:25:53.500000", "2023-10-14T18:27:16.500000", "2023-10-14T18:28:38.500000", "2023-10-14T18:30:01.500000", "2023-10-14T18:31:24.500000", "2023-10-14T18:32:47.500000", "2023-10-14T18:34:09.500000", "2023-10-14T18:35:32.500000", "2023-10-14T18:36:55.500000", "2023-10-14T18:38:17.500000", "2023-10-14T18:39:40.500000", "2023-10-14T18:41:03.500000", "2023-10-14T18:42:25.500000", "2023-10-14T18:43:48.500000", "2023-10-14T18:45:11.500000", "2023-10-14T18:46:33.500000", "2023-10-14T18:47:56.500000", "2023-10-14T18:49:19.500000", "2023-10-14T18:50:41.500000", "2023-10-14T18:52:04.500000", "2023-10-14T18:53:27.500000", "2023-10-14T18:54:49.500000", "2023-10-14T18:56:12.500000", "2023-10-14T18:57:35.500000", "2023-10-14T18:58:57.500000", "2023-10-14T19:00:20.500000", "2023-10-14T19:01:43.500000", "2023-10-14T19:03:06.500000"]
rot = []
for i in neid_october:
    date_object = datetime.strptime(i, "%Y-%m-%dT%H:%M:%S.%f")
    rot.append(solar_orientation(location, time = date_object, coordinates="altaz")[5])
import h5py
with h5py.File('rotation_matrices_sxform.jld2', 'w') as f:
    for i, matrix in enumerate(rot):
        # Save each matrix to the file with a unique name (e.g., 'matrix0', 'matrix1', etc.)
        f.create_dataset(f'matrix{i+1}', data=matrix)

from datetime import datetime, timedelta, timezone
initial_epoch = datetime(2023,10,14,15,0,0, tzinfo=timezone.utc)
final_epoch = datetime(2023,10,14,18,10,0, tzinfo=timezone.utc)
solar_orientation_altaz = []
solar_orientation_equat = []
AltAz_arr = []
IAU_Earth_spice = []
IAU_Earth_sunpy = []
sunpy_angle = []
timestamp = []
import sunpy.coordinates
while initial_epoch < final_epoch:
    solar_orientation_altaz.append(solar_orientation(location, time = initial_epoch, coordinates="altaz")[0].degree)
    solar_orientation_equat.append(solar_orientation(location, time = initial_epoch, coordinates="equat")[0].degree)
    AltAz_arr.append(solar_orientation(location, time = initial_epoch, coordinates="altaz")[2])
    IAU_Earth_spice.append(solar_orientation(location, time = initial_epoch, coordinates="altaz")[3])
    IAU_Earth_sunpy.append(solar_orientation(location, time = initial_epoch, coordinates="altaz")[4])
    sunpy_angle.append(sunpy.coordinates.sun.orientation(location, time = initial_epoch).degree)
    timestamp.append(initial_epoch)
    initial_epoch  = initial_epoch + timedelta(minutes=5)
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
mpl = plt.matplotlib 

fig = plt.figure()
ax1 = fig.add_subplot()
ax1.scatter(timestamp, np.array(solar_orientation_altaz) + 90, s = 5, label = "Sunpy AltAz ($\phi$)")
ax1.scatter(timestamp, np.array(AltAz_arr) + 90, s = 5, label = "Spice AltAz ($\phi$)")
ax1.scatter(timestamp, 90 - np.array(IAU_Earth_sunpy), s = 5, label = "Sunpy IAU_EARTH (full)")
ax1.scatter(timestamp, 90 - np.array(IAU_Earth_spice), s = 5, label = "Spice IAU_EARTH (full)")
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
plt.xlabel("Time")
plt.ylabel("Angle (degrees)")
plt.legend()
plt.savefig("solar_north_horizon_new.png")
plt.show()

#-------------------------------------------------------------

# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.scatter(timestamp, sunpy_angle, label = "sunpy - Φ")
# ax1.scatter(timestamp, solar_orientation_altaz, label = "solar_orientation - altaz", s = 5)

# #SPICE
# spice_angle_array = [[26.130000000000003], [26.130000000000003], [26.13000000000002], [26.129999999999992], [26.130000000000003], [26.13000000000002], [26.13000000000002], [26.129999999999992], [26.129999999999992], [26.13000000000002], [26.13000000000002], [26.13000000000003], [26.13000000000002], [26.129999999999992], [26.13000000000002], [26.13000000000003], [26.129999999999992], [26.129999999999992], [26.130000000000003], [26.129999999999992], [26.130000000000003], [26.13000000000002], [26.130000000000003], [26.129999999999974], [26.129999999999992], [26.13000000000002], [26.13000000000002], [26.129999999999992], [26.129999999999992], [26.130000000000003], [26.130000000000003], [26.13000000000002], [26.13000000000002], [26.129999999999992], [26.129999999999992], [26.13000000000002], [26.13000000000002], [26.129999999999992]]
# spice_angle = []
# for i in range(0,len(spice_angle_array)):
#     spice_angle.append(-(spice_angle_array[i][0]))
# ax1.scatter(timestamp, spice_angle, label = "SPICE - Φ")
# ax1.scatter(timestamp, solar_orientation_equat, label = "solar_orientation - equat", s = 5)
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# plt.xlabel("Time")
# plt.ylabel("Angle")
# plt.legend()
# plt.savefig("tilt_confirmed.png")
# plt.show()

#-------------------------------------------------------------
# from astropy.coordinates import EarthLocation
# #EXPRES
# obs_lat_ex = 34.744444
# obs_long_ex = -111.421944 
# alt_ex = 2.359152
# location_ex = EarthLocation.from_geodetic(obs_long_ex, obs_lat_ex, alt_ex)
# #Boulder
# obs_lat_bo = 39.995380
# obs_long_bo = 360-105.262390
# alt_bo = 1.6523
# location_bo = EarthLocation.from_geodetic(obs_long_bo, obs_lat_bo, alt_bo)
# #NEID
# obs_lat_neid = 31.9583 
# obs_long_neid = 360-111.5967  
# alt_neid = 2.097938
# location_neid = EarthLocation.from_geodetic(obs_long_neid, obs_lat_neid, alt_neid)
# from datetime import datetime, timedelta, timezone

# initial_epoch_eclipse = datetime(2023,10,14,15,0,0, tzinfo=timezone.utc)
# final_epoch_eclipse = datetime(2023,10,14,18,10,0, tzinfo=timezone.utc)
# solar_orientation_phi_eclipse_ex = []
# solar_orientation_theta_eclipse_ex = []
# solar_orientation_phi_eclipse_bo = []
# solar_orientation_theta_eclipse_bo = []
# solar_orientation_phi_eclipse_neid = []
# solar_orientation_theta_eclipse_neid = []
# timestamp_eclipse = []
# while initial_epoch_eclipse < final_epoch_eclipse:
#     solar_orientation_phi_eclipse_ex.append(solar_orientation(location_ex, time = initial_epoch_eclipse)[0].degree)
#     solar_orientation_theta_eclipse_ex.append(solar_orientation(location_ex, time = initial_epoch_eclipse)[1].degree)

#     solar_orientation_phi_eclipse_bo.append(solar_orientation(location_bo, time = initial_epoch_eclipse)[0].degree)
#     solar_orientation_theta_eclipse_bo.append(solar_orientation(location_bo, time = initial_epoch_eclipse)[1].degree)

#     solar_orientation_phi_eclipse_neid.append(solar_orientation(location_neid, time = initial_epoch_eclipse)[0].degree)
#     solar_orientation_theta_eclipse_neid.append(solar_orientation(location_neid, time = initial_epoch_eclipse)[1].degree)
    
#     timestamp_eclipse.append(initial_epoch_eclipse)
#     initial_epoch_eclipse  = initial_epoch_eclipse + timedelta(minutes=5)


# initial_epoch = datetime(2023,10,14,7,0,0, tzinfo=timezone.utc)
# final_epoch = datetime(2023,10,14,20,0,0, tzinfo=timezone.utc)
# solar_orientation_phi_ex = []
# solar_orientation_theta_ex = []
# solar_orientation_phi_bo = []
# solar_orientation_theta_bo = []
# solar_orientation_phi_neid = []
# solar_orientation_theta_neid = []
# timestamp = []
# while initial_epoch < final_epoch:
#     solar_orientation_phi_ex.append(solar_orientation(location_ex, time = initial_epoch)[0].degree)
#     solar_orientation_theta_ex.append(solar_orientation(location_ex, time = initial_epoch)[1].degree)

#     solar_orientation_phi_bo.append(solar_orientation(location_bo, time = initial_epoch)[0].degree)
#     solar_orientation_theta_bo.append(solar_orientation(location_bo, time = initial_epoch)[1].degree)

#     solar_orientation_phi_neid.append(solar_orientation(location_neid, time = initial_epoch)[0].degree)
#     solar_orientation_theta_neid.append(solar_orientation(location_neid, time = initial_epoch)[1].degree)
    
#     timestamp.append(initial_epoch)
#     initial_epoch  = initial_epoch + timedelta(minutes=5)    

# import matplotlib.pyplot as plt
# import matplotlib.dates as mdates
# mpl = plt.matplotlib 
# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.plot(timestamp, solar_orientation_phi_ex, color = 'b', label = "EXPRES") 
# ax1.plot(timestamp, solar_orientation_theta_ex, color = 'grey', linestyle='dashed') 
# ax1.plot(timestamp, solar_orientation_phi_bo, color = 'r', label = "heterodyne") 
# ax1.plot(timestamp, solar_orientation_phi_neid, color = 'g', label = "NEID") 

# ax1.plot(timestamp_eclipse, solar_orientation_phi_eclipse_ex, color = 'k') 
# ax1.plot(timestamp_eclipse, solar_orientation_theta_eclipse_ex, color = 'k', linestyle='dashed') 
# ax1.plot(timestamp_eclipse, solar_orientation_phi_eclipse_bo, color = 'k') 
# ax1.plot(timestamp_eclipse, solar_orientation_phi_eclipse_neid, color = 'k') 

# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# plt.xlabel("UTC Time on 10/14/2023")
# plt.ylabel("Angle (degree)")
# plt.title("Φ - solid line, θ - dashed line")
# plt.legend()
# plt.savefig("tilt_evolution.png")
# plt.show()