""" Create a rotating Earth; create the exact change in tilt, wobble, angular velocity based on how Earth is currently behaving"""

'''Create a rotating, titled Earth'''

import utils as ut
import numpy as np #generating data
import matplotlib.pyplot as plt #for plotting
import matplotlib.ticker as mtk
from mpl_toolkits import mplot3d #for 3D plots
import math as mt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

'''User input lat-long'''
target_latitude = 34.0467
target_longitude = -118.5464

'''Plot all points on Earth'''
ax = plt.figure(figsize=(10, 10)).add_subplot(111, projection="3d")
radius_earth = 1.0
latitude_data = np.linspace(-np.pi/2, np.pi/2, 30) #(pi/2)-polar angle (or elevation angle) in radians
longitude_data = np.linspace(-np.pi, np.pi, 30) #azimuth in radians
x_data = np.array([])
y_data = np.array([])
z_data = np.array([])
'''
traverse the length of each longitude, and come up with cartesian co-ords of each point on that longtitude
same can be done for each latitude instead. objective is to find the x, y, z co-ods of all points on the Earth's surface.
'''
for i in range(len(longitude_data)):
  for j in range(len(latitude_data)):
    x_data = np.append(x_data, radius_earth * mt.cos(latitude_data[j]) * mt.cos(longitude_data[i]))
    y_data = np.append(y_data, radius_earth * mt.cos(latitude_data[j]) * mt.sin(longitude_data[i]))
    z_data = np.append(z_data, radius_earth * mt.sin(latitude_data[j]))
unit_array = np.ones(np.size(longitude_data))
z_data = np.outer(unit_array, z_data)
ax.plot_wireframe(x_data, y_data, z_data, alpha=0.002, color='black', rstride=2, cstride=2)
ax.set_xlabel("X axis")
ax.set_ylabel("Y axis")
ax.set_zlabel("Z axis")

'''plot continents'''
proj = ccrs.PlateCarree() # Use Cartopy's PlateCarree projection for geographic data
continents = cfeature.NaturalEarthFeature(
    category='physical', name='land', scale='110m', facecolor='none', edgecolor='black'
) # Get continent boundaries from Cartopy

# Iterate through continent geometries
for geom in continents.geometries():
    # Handle both Polygon and MultiPolygon
    if geom.geom_type == 'Polygon':
        # Process single Polygon
        coords = np.array(geom.exterior.coords)
        lon, lat = coords[:, 0], coords[:, 1]
        # Convert latitude/longitude to spherical coordinates
        theta = np.radians(90 - lat)
        phi = np.radians(lon)
        # Convert to Cartesian coordinates for 3D plotting
        x = np.cos(phi) * np.sin(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(theta)
        # Plot the continent outline
        ax.plot(x, y, z, color='maroon', linewidth=1)
    elif geom.geom_type == 'MultiPolygon':
        # Process each Polygon in MultiPolygon
        for poly in geom.geoms:  # Use .geoms to access individual polygons
            coords = np.array(poly.exterior.coords)
            lon, lat = coords[:, 0], coords[:, 1]          
            # Convert latitude/longitude to spherical coordinates
            theta = np.radians(90 - lat)
            phi = np.radians(lon)
            # Convert to Cartesian coordinates for 3D plotting
            x = np.cos(phi) * np.sin(theta)
            y = np.sin(phi) * np.sin(theta)
            z = np.cos(theta)
            # Plot the continent outline
            ax.plot(x, y, z, color='maroon', linewidth=1)

'''A line passing through the north and south poles'''
x_poles = np.array([ut.spherical_to_cartesian(1, -90, 0)[0], ut.spherical_to_cartesian(1, 90, 0)[0]])
y_poles = np.array([ut.spherical_to_cartesian(1, -90, 0)[1], ut.spherical_to_cartesian(1, 90, 0)[1]])
z_poles = np.array([ut.spherical_to_cartesian(1, -90, 0)[2], ut.spherical_to_cartesian(1, 90, 0)[2]])
ax.plot(x_poles, y_poles, z_poles, color='darkblue', alpha=1, linestyle="--", linewidth=1.5)

'''Plot satellite orbit'''
radius_sat_orbit = 1.75 * radius_earth # radius of LEO
x_data_sat = np.array([])
y_data_sat = np.array([])
z_data_sat = np.array([])
for i in range(len(longitude_data)):
  for j in range(len(latitude_data)):
    x_data_sat = np.append(x_data_sat, radius_sat_orbit * mt.cos(latitude_data[j]) * mt.cos(longitude_data[i]))
    y_data_sat = np.append(y_data_sat, radius_sat_orbit * mt.cos(latitude_data[j]) * mt.sin(longitude_data[i]))
    z_data_sat = np.append(z_data_sat, radius_sat_orbit * mt.sin(latitude_data[j]))
z_data_sat = np.outer(unit_array, z_data_sat)
ax.plot_wireframe(x_data_sat, y_data_sat, z_data_sat, alpha=0.05, color='gray', rstride=1, cstride=1, linewidth = 0.5)

'''plot a vector from the origin to the LEO point passing through the fireprone area'''
ax.quiver(0, 0, 0, 
          ut.spherical_to_cartesian(radius_sat_orbit, target_latitude, target_longitude)[0], 
          ut.spherical_to_cartesian(radius_sat_orbit, target_latitude, target_longitude)[1], 
          ut.spherical_to_cartesian(radius_sat_orbit, target_latitude, target_longitude)[2], arrow_length_ratio=0.2, color='black')

'''plot a randomly oriented cubesat at the LEO above the fire prone area'''
# Plot a spacecraft body reference frame which is randomly oriented
scb_axis_length = 1.0 * radius_earth
x_scb_origin = ut.spherical_to_cartesian(radius_sat_orbit, target_latitude, target_longitude)[0]
y_scb_origin = ut.spherical_to_cartesian(radius_sat_orbit, target_latitude, target_longitude)[1]
z_scb_origin = ut.spherical_to_cartesian(radius_sat_orbit, target_latitude, target_longitude)[2]

ax.quiver(x_scb_origin, y_scb_origin, z_scb_origin, x_scb_origin + scb_axis_length, 0, 0, arrow_length_ratio=0.3, color='red')
ax.quiver(x_scb_origin, y_scb_origin, z_scb_origin, 0, y_scb_origin + scb_axis_length, 0, arrow_length_ratio=0.3, color='green')
ax.quiver(x_scb_origin, y_scb_origin, z_scb_origin, 0, 0, z_scb_origin + scb_axis_length, arrow_length_ratio=0.3, color='blue')


'''display chart'''
plt.show()
