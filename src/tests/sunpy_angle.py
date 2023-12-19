import numpy as np
from sunpy.sun import constants
from sunpy.time import parse_time
from sunpy.coordinates.frames import HeliographicStonyhurst
import astropy.units as u
from astropy.coordinates import AltAz, HADec, SkyCoord, Angle
from astropy.coordinates.representation import CartesianRepresentation, SphericalRepresentation
from astropy.modeling.rotations import SphericalRotationSequence

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
    
    # Return angles relative to y axis in image plane
    return Angle(phi) - 90*u.deg, Angle(theta) - 90*u.deg

#-------------------------------------------------------------

# from astropy.coordinates import EarthLocation
# obs_lat = 39.995380
# obs_long = -105.262390
# alt = 1.6523
# location = EarthLocation.from_geodetic(obs_long, obs_lat, alt)
# from datetime import datetime, timedelta, timezone
# initial_epoch = datetime(2023,10,14,15,0,0, tzinfo=timezone.utc)
# final_epoch = datetime(2023,10,14,18,10,0, tzinfo=timezone.utc)
# solar_orientation_altaz = []
# solar_orientation_equat = []
# sunpy_angle = []
# timestamp = []
# import sunpy.coordinates
# while initial_epoch < final_epoch:
#     solar_orientation_altaz.append(solar_orientation(location, time = initial_epoch, coordinates="altaz")[0].degree)
#     solar_orientation_equat.append(solar_orientation(location, time = initial_epoch, coordinates="equat")[0].degree)
#     sunpy_angle.append(sunpy.coordinates.sun.orientation(location, time = initial_epoch).degree)
#     timestamp.append(initial_epoch)
#     initial_epoch  = initial_epoch + timedelta(minutes=5)
# import matplotlib.pyplot as plt
# import matplotlib.dates as mdates
# mpl = plt.matplotlib 
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
from astropy.coordinates import EarthLocation
#EXPRES
obs_lat_ex = 34.744444
obs_long_ex = -111.421944 
alt_ex = 2.359152
location_ex = EarthLocation.from_geodetic(obs_long_ex, obs_lat_ex, alt_ex)
#Boulder
obs_lat_bo = 39.995380
obs_long_bo = 360-105.262390
alt_bo = 1.6523
location_bo = EarthLocation.from_geodetic(obs_long_bo, obs_lat_bo, alt_bo)
#NEID
obs_lat_neid = 31.9583 
obs_long_neid = 360-111.5967  
alt_neid = 2.097938
location_neid = EarthLocation.from_geodetic(obs_long_neid, obs_lat_neid, alt_neid)
from datetime import datetime, timedelta, timezone

initial_epoch_eclipse = datetime(2023,10,14,15,0,0, tzinfo=timezone.utc)
final_epoch_eclipse = datetime(2023,10,14,18,10,0, tzinfo=timezone.utc)
solar_orientation_phi_eclipse_ex = []
solar_orientation_theta_eclipse_ex = []
solar_orientation_phi_eclipse_bo = []
solar_orientation_theta_eclipse_bo = []
solar_orientation_phi_eclipse_neid = []
solar_orientation_theta_eclipse_neid = []
timestamp_eclipse = []
while initial_epoch_eclipse < final_epoch_eclipse:
    solar_orientation_phi_eclipse_ex.append(solar_orientation(location_ex, time = initial_epoch_eclipse)[0].degree)
    solar_orientation_theta_eclipse_ex.append(solar_orientation(location_ex, time = initial_epoch_eclipse)[1].degree)

    solar_orientation_phi_eclipse_bo.append(solar_orientation(location_bo, time = initial_epoch_eclipse)[0].degree)
    solar_orientation_theta_eclipse_bo.append(solar_orientation(location_bo, time = initial_epoch_eclipse)[1].degree)

    solar_orientation_phi_eclipse_neid.append(solar_orientation(location_neid, time = initial_epoch_eclipse)[0].degree)
    solar_orientation_theta_eclipse_neid.append(solar_orientation(location_neid, time = initial_epoch_eclipse)[1].degree)
    
    timestamp_eclipse.append(initial_epoch_eclipse)
    initial_epoch_eclipse  = initial_epoch_eclipse + timedelta(minutes=5)


initial_epoch = datetime(2023,10,14,7,0,0, tzinfo=timezone.utc)
final_epoch = datetime(2023,10,14,20,0,0, tzinfo=timezone.utc)
solar_orientation_phi_ex = []
solar_orientation_theta_ex = []
solar_orientation_phi_bo = []
solar_orientation_theta_bo = []
solar_orientation_phi_neid = []
solar_orientation_theta_neid = []
timestamp = []
while initial_epoch < final_epoch:
    solar_orientation_phi_ex.append(solar_orientation(location_ex, time = initial_epoch)[0].degree)
    solar_orientation_theta_ex.append(solar_orientation(location_ex, time = initial_epoch)[1].degree)

    solar_orientation_phi_bo.append(solar_orientation(location_bo, time = initial_epoch)[0].degree)
    solar_orientation_theta_bo.append(solar_orientation(location_bo, time = initial_epoch)[1].degree)

    solar_orientation_phi_neid.append(solar_orientation(location_neid, time = initial_epoch)[0].degree)
    solar_orientation_theta_neid.append(solar_orientation(location_neid, time = initial_epoch)[1].degree)
    
    timestamp.append(initial_epoch)
    initial_epoch  = initial_epoch + timedelta(minutes=5)    

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
mpl = plt.matplotlib 
fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(timestamp, solar_orientation_phi_ex, color = 'b', label = "EXPRES") 
ax1.plot(timestamp, solar_orientation_theta_ex, color = 'grey', linestyle='dashed') 
ax1.plot(timestamp, solar_orientation_phi_bo, color = 'r', label = "heterodyne") 
ax1.plot(timestamp, solar_orientation_phi_neid, color = 'g', label = "NEID") 

ax1.plot(timestamp_eclipse, solar_orientation_phi_eclipse_ex, color = 'k') 
ax1.plot(timestamp_eclipse, solar_orientation_theta_eclipse_ex, color = 'k', linestyle='dashed') 
ax1.plot(timestamp_eclipse, solar_orientation_phi_eclipse_bo, color = 'k') 
ax1.plot(timestamp_eclipse, solar_orientation_phi_eclipse_neid, color = 'k') 

ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
plt.xlabel("UTC Time on 10/14/2023")
plt.ylabel("Angle (degree)")
plt.title("Φ - solid line, θ - dashed line")
plt.legend()
plt.savefig("tilt_evolution.png")
plt.show()