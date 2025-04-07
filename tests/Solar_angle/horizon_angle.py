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
from astropy.coordinates import EarthLocation
import astropy.units as u
from datetime import datetime, timedelta, timezone
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
mpl = plt.matplotlib 

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
    sun_coord_frame = (sun_coord.transform_to(frame)).data.to_cartesian() - sun_center.data.to_cartesian()
    AltAz_angle = np.arctan2(sun_coord_frame.dot(im_y), sun_coord_frame.dot(im_x)).to('deg')

    #-------------------------------------------------------------
    
    #IAU_Earth Frame Spice
    rot = spiceypy.pxform("IAU_SUN", "IAU_Earth", spiceypy.utc2et(time.strftime("%Y-%m-%dT%H:%M:%S")))
    transformed = np.dot(rot, sun_coord.cartesian.xyz.value)
    transformed = CartesianRepresentation(x=transformed[0] * u.km, y=transformed[1] * u.km, z=transformed[2] * u.km)
    sun_coord_transformed = SkyCoord(transformed, frame="spice_IAU_EARTH", obstime=frame.obstime)
    sun_coord_frame = (sun_coord_transformed.transform_to(frame).data.to_cartesian())
    IAU_Earth_angle_spice = np.arctan2(sun_coord_transformed.data.to_cartesian().dot(im_y), sun_coord_transformed.data.to_cartesian().dot(im_x)).to('deg')
    # IAU_Earth_angle_spice = np.arctan2(sun_coord_frame.dot(im_y), sun_coord_frame.dot(im_x)).to('deg')

    #IAU_Earth Frame Sunpy
    sun_center_IAU_EARTH = sun_center.transform_to("spice_IAU_EARTH")
    sun_coord_transformed = sun_coord.transform_to("spice_IAU_EARTH")
    sun_coord_frame = sun_coord_transformed.data.to_cartesian() - sun_center_IAU_EARTH.data.to_cartesian()
    # sun_coord_frame = (sun_coord_transformed.transform_to(frame).data.to_cartesian() - sun_center.data.to_cartesian())
    IAU_Earth_angle_sunpy = np.arctan2(sun_coord_frame.dot(im_y), sun_coord_frame.dot(im_x)).to('deg')
     
    # Return angles relative to y axis in image plane
    return (Angle(phi) - 90*u.deg).degree, (Angle(AltAz_angle) - 90*u.deg).degree, (Angle(IAU_Earth_angle_spice) - 90*u.deg).degree, (Angle(IAU_Earth_angle_sunpy) - 90*u.deg).degree

obs_lat = 31.9583 
obs_long = -111.5967  
alt = 2.097938
location = EarthLocation.from_geodetic(obs_long, obs_lat, alt)

sun_radius = 1.0 * u.R_sun
sun_axis_sun = [0.0, 0.0, sun_radius.to(u.km).value]
cartesian_rep = CartesianRepresentation(x=sun_axis_sun[0] * u.km, y=sun_axis_sun[1] * u.km, z=sun_axis_sun[2] * u.km)

initial_epoch = datetime(2023,10,14,15,0,0, tzinfo=timezone.utc)
final_epoch = datetime(2023,10,14,18,10,0, tzinfo=timezone.utc)
length = 38
solar_orientation_altaz = np.zeros(length) 
AltAz_arr = np.zeros(length) 
IAU_Earth_spice = np.zeros(length) 
IAU_Earth_sunpy = np.zeros(length) 
timestamp = [] 
count = 0
while initial_epoch < final_epoch:
    solar_orientation_altaz[count], AltAz_arr[count], IAU_Earth_spice[count], IAU_Earth_sunpy[count] = solar_orientation(location, time = initial_epoch, coordinates="altaz")
    timestamp.append(initial_epoch)
    initial_epoch  = initial_epoch + timedelta(minutes=5)
    count += 1

fig = plt.figure()
ax1 = fig.add_subplot()
# ax1.plot(timestamp, np.array(solar_orientation_altaz) + 90, label = "Sunpy AltAz ($\phi$)")
# ax1.plot(timestamp, np.array(AltAz_arr) + 90, label = "Spice AltAz ($\phi$)")
ax1.scatter(timestamp, np.array(IAU_Earth_spice) + 90, s = 10, label = "Spice IAU_EARTH ($\phi$)")
ax1.scatter(timestamp, np.array(IAU_Earth_sunpy) + 90, s = 10, label = "Sunpy IAU_EARTH ($\phi$)")
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
plt.xlabel("Time")
plt.ylabel("Angle (degrees)")
plt.legend()
plt.savefig("validate_spice_rot.png")
plt.show()
plt.clf()