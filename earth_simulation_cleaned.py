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
import random as rd
from scipy.spatial.transform import Rotation as R

'''User input lat-long''' #18.5246° N, 73.8786° E; 34.0467° N, 118.5464° W
target_latitude = 34.0467
target_longitude = -118.5464

'''Plot all points on Earth'''
ax = plt.figure(figsize=(7, 7)).add_subplot(111, projection="3d")
radius_earth = 10.0
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
ax.plot_wireframe(x_data, y_data, z_data, alpha=0.015, color='black', rstride=2, cstride=2, linewidth=0.5)
ax.set_xlabel("X axis")
ax.set_ylabel("Y axis")
ax.set_zlabel("Z axis")

'''plot ECEF reference frame'''
ecef_axis_length = 1.5 * radius_earth
ax.quiver(0, 0, 0, ecef_axis_length, 0, 0, arrow_length_ratio=0.3, color='red')
ax.quiver(0, 0, 0, 0, ecef_axis_length, 0, arrow_length_ratio=0.3, color='blue')
ax.quiver(0, 0, 0, 0, 0, ecef_axis_length, arrow_length_ratio=0.3, color='green')
ax.text(0 + ecef_axis_length, 0, 0, 'Xecef', color='red')
ax.text(0, 0 + ecef_axis_length, 0, 'Yecef', color='blue')
ax.text(0, 0, 0 + ecef_axis_length, 'Zecef', color='green')

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
        x = radius_earth * np.cos(phi) * np.sin(theta)
        y = radius_earth * np.sin(phi) * np.sin(theta)
        z = radius_earth * np.cos(theta)
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
            x = radius_earth * np.cos(phi) * np.sin(theta)
            y = radius_earth * np.sin(phi) * np.sin(theta)
            z = radius_earth * np.cos(theta)
            # Plot the continent outline
            ax.plot(x, y, z, color='maroon', linewidth=1)

'''Plot LEO sphere'''
radius_sat_orbit = 1.75 * radius_earth # radius of LEO
lat, long = np.meshgrid(latitude_data, longitude_data)
x_data_sat = radius_sat_orbit * np.cos(lat) * np.cos(long)
y_data_sat = radius_sat_orbit * np.cos(lat) * np.sin(long)
z_data_sat = radius_sat_orbit * np.sin(lat)
# Plotting the 3D sphere
ax.plot_surface(x_data_sat, y_data_sat, z_data_sat, color='white', alpha=0.15)

'''plot a vector from the origin to the LEO point passing through the fireprone area'''
ax.quiver(0, 0, 0, 
          ut.spherical_to_cartesian(radius_sat_orbit, target_latitude, target_longitude)[0], 
          ut.spherical_to_cartesian(radius_sat_orbit, target_latitude, target_longitude)[1], 
          ut.spherical_to_cartesian(radius_sat_orbit, target_latitude, target_longitude)[2], arrow_length_ratio=0.2, color='black')

# '''plot a randomly / specifically oriented cubesat at the LEO above the fire prone area'''
scb_axis_length = 0.75 * radius_earth
scb_origin = [ut.spherical_to_cartesian(radius_sat_orbit, target_latitude, target_longitude)[0], # x_scb_origin
              ut.spherical_to_cartesian(radius_sat_orbit, target_latitude, target_longitude)[1], # y_scb_origin
              ut.spherical_to_cartesian(radius_sat_orbit, target_latitude, target_longitude)[2]] # z_scb_origin
scb_ref_frame_start = np.array([scb_origin, scb_origin, scb_origin])
rotation = R.random() # random rotation object
#rotation = R.from_euler('xyz', [135, 0, 35], degrees=True) # specific rotation object
scb_ref_frame_rotated = np.array([rotation.apply([scb_axis_length, 0, 0]),
                                  rotation.apply([0, scb_axis_length, 0]),
                                  rotation.apply([0, 0, scb_axis_length])])
scb_coods = scb_ref_frame_start + scb_ref_frame_rotated
ax.quiver(scb_ref_frame_start[:, 0], scb_ref_frame_start[:, 1], scb_ref_frame_start[:, 2], 
          scb_ref_frame_rotated[:, 0], scb_ref_frame_rotated[:, 1], scb_ref_frame_rotated[:, 2],
          arrow_length_ratio=0.3, linestyle=":",
          color=['red', 'blue', 'green'])
ax.text(scb_origin[0] + rotation.apply([scb_axis_length, 0, 0])[0], 
        scb_origin[1] + rotation.apply([scb_axis_length, 0, 0])[1], 
        scb_origin[2] + rotation.apply([scb_axis_length, 0, 0])[2], 
        'Xscb', color='red')
ax.text(scb_origin[0] + rotation.apply([0, scb_axis_length, 0])[0], 
        scb_origin[1] + rotation.apply([0, scb_axis_length, 0])[1], 
        scb_origin[2] + rotation.apply([0, scb_axis_length, 0])[2], 
        'Yscb', color='blue')
ax.text(scb_origin[0] + rotation.apply([0, 0, scb_axis_length])[0], 
        scb_origin[1] + rotation.apply([0, 0, scb_axis_length])[1], 
        scb_origin[2] + rotation.apply([0, 0, scb_axis_length])[2], 
        'Zscb', color='green')

'''plot a copy of the ECEF reference frame at the SCB origin'''
ax.quiver(scb_origin[0], scb_origin[1], scb_origin[2], ecef_axis_length, 0, 0, arrow_length_ratio=0.3, color='red')
ax.quiver(scb_origin[0], scb_origin[1], scb_origin[2], 0, ecef_axis_length, 0, arrow_length_ratio=0.3, color='blue')
ax.quiver(scb_origin[0], scb_origin[1], scb_origin[2], 0, 0, ecef_axis_length, arrow_length_ratio=0.3, color='green')
ax.text(scb_origin[0] + ecef_axis_length, scb_origin[1], scb_origin[2], 'Xecef', color='red')
ax.text(scb_origin[0], scb_origin[1] + ecef_axis_length, scb_origin[2], 'Yecef', color='blue')
ax.text(scb_origin[0], scb_origin[1], scb_origin[2] + ecef_axis_length, 'Zecef', color='green')

'''Euler angles calculations'''
'''
Rotate the reference frame by Xscb_azimuth about the Zecef axis
Rotate the reference frame by Xscb_elevation about the Yecef axis
Rotate the reference frame by polar_Zscb_for_rotation_about_Xecef about the Xecef axis
'''
azimuth = (ut.cartesian_to_spherical(scb_origin[0], scb_origin[1], scb_origin[2],
                                    scb_coods[0, 0], scb_coods[0, 1], scb_coods[0, 2])[1]) * (180/mt.pi)
polar = (ut.cartesian_to_spherical(scb_origin[0], scb_origin[1], scb_origin[2],
                                    scb_coods[0, 0], scb_coods[0, 1], scb_coods[0, 2])[2]) * (180/mt.pi)

'''DCM logic'''
# Top 4 spaces
if scb_ref_frame_rotated[0, 0] >= 0 and scb_ref_frame_rotated[0, 1] >= 0 and scb_ref_frame_rotated[0, 2] >= 0: # X+ve, Y+ve, Z+ve
   azim_rotation_of_Xscb = -azimuth
   elev_rotation_of_Xscb = (90 - polar)
elif scb_ref_frame_rotated[0, 0] < 0 and scb_ref_frame_rotated[0, 1] >= 0 and scb_ref_frame_rotated[0, 2] >= 0: # X-ve, Y+ve, Z+ve
   azim_rotation_of_Xscb = -(180 - azimuth)
   elev_rotation_of_Xscb = (90 - polar)
elif scb_ref_frame_rotated[0, 0] >= 0 and scb_ref_frame_rotated[0, 1] < 0 and scb_ref_frame_rotated[0, 2] >= 0: # X+ve, Y-ve, Z+ve
   azim_rotation_of_Xscb = azimuth
   elev_rotation_of_Xscb = (90 - polar)
elif scb_ref_frame_rotated[0, 0] < 0 and scb_ref_frame_rotated[0, 1] < 0 and scb_ref_frame_rotated[0, 2] >= 0: # X-ve, Y-ve, Z+ve
   azim_rotation_of_Xscb = (180 - azimuth)
   elev_rotation_of_Xscb = (90 - polar)
# Bottom 4 spaces   
elif scb_ref_frame_rotated[0, 0] >= 0 and scb_ref_frame_rotated[0, 1] >= 0 and scb_ref_frame_rotated[0, 2] < 0: # X+ve, Y+ve, Z-ve
   azim_rotation_of_Xscb = -azimuth
   elev_rotation_of_Xscb = -(90 - polar)
elif scb_ref_frame_rotated[0, 0] < 0 and scb_ref_frame_rotated[0, 1] >= 0 and scb_ref_frame_rotated[0, 2] < 0: # X-ve, Y+ve, Z-ve
   azim_rotation_of_Xscb = -(180 - azimuth)
   elev_rotation_of_Xscb = -(90 - polar)
elif scb_ref_frame_rotated[0, 0] >= 0 and scb_ref_frame_rotated[0, 1] < 0 and scb_ref_frame_rotated[0, 2] < 0: # X+ve, Y-ve, Z-ve
   azim_rotation_of_Xscb = azimuth
   elev_rotation_of_Xscb = -(90 - polar)
elif scb_ref_frame_rotated[0, 0] < 0 and scb_ref_frame_rotated[0, 1] < 0 and scb_ref_frame_rotated[0, 2] < 0: # X-ve, Y-ve, Z-ve
   azim_rotation_of_Xscb = (180 - azimuth)
   elev_rotation_of_Xscb = -(90 - polar)

azim_rotation_of_Xscb = azim_rotation_of_Xscb * (mt.pi/180)
elev_rotation_of_Xscb = elev_rotation_of_Xscb * (mt.pi/180)
dcm_rotation_abt_Zscb = np.array([[mt.cos(azim_rotation_of_Xscb), -mt.sin(azim_rotation_of_Xscb), 0],
                                  [mt.sin(azim_rotation_of_Xscb), mt.cos(azim_rotation_of_Xscb), 0],
                                  [0, 0, 1]])
rotation_about_Zscb = np.array(np.matmul(dcm_rotation_abt_Zscb, scb_ref_frame_rotated.T)).T
# ax.quiver(scb_ref_frame_start[:, 0], scb_ref_frame_start[:, 1], scb_ref_frame_start[:, 2], 
#           rotation_about_Zscb[:, 0], rotation_about_Zscb[:, 1], rotation_about_Zscb[:, 2],
#           arrow_length_ratio=0.4, linestyle=":",
#           color=['red', 'blue', 'green'])

dcm_rotation_abt_Yscb = np.array([[mt.cos(elev_rotation_of_Xscb), 0, mt.sin(elev_rotation_of_Xscb)],
                                  [0, 1, 0],
                                  [-mt.sin(elev_rotation_of_Xscb), 0, mt.cos(elev_rotation_of_Xscb)]])
rotation_about_Yscb = np.array(np.matmul(dcm_rotation_abt_Yscb, rotation_about_Zscb.T)).T
scb_coods_Zscb_rotated_twice = scb_ref_frame_start + rotation_about_Yscb
# ax.quiver(scb_ref_frame_start[:, 0], scb_ref_frame_start[:, 1], scb_ref_frame_start[:, 2], 
#           rotation_about_Yscb[:, 0], rotation_about_Yscb[:, 1], rotation_about_Yscb[:, 2],
#           arrow_length_ratio=0.4, linestyle="-.",
#           color=['red', 'blue', 'green'])

polar_Zscb_for_rotation_about_Xecef = (ut.cartesian_to_spherical(scb_origin[0], scb_origin[1], scb_origin[2], 
                                                                scb_coods_Zscb_rotated_twice[2, 0], scb_coods_Zscb_rotated_twice[2, 1], scb_coods_Zscb_rotated_twice[2, 2])[2]) * (180/mt.pi)

print(f"Polar of Zscb before final rotation: {polar_Zscb_for_rotation_about_Xecef}")
print(f"rotation_about_Yscb: {rotation_about_Yscb}")
# Xecef positive
if rotation_about_Yscb[2, 0] >= 0 and rotation_about_Yscb[2, 1] >= 0 and rotation_about_Yscb[2, 2] >= 0:
   polar_Zscb_for_rotation_about_Xecef = polar_Zscb_for_rotation_about_Xecef
elif rotation_about_Yscb[2, 0] >= 0 and rotation_about_Yscb[2, 1] >= 0 and rotation_about_Yscb[2, 2] < 0:
   polar_Zscb_for_rotation_about_Xecef = (180 - polar_Zscb_for_rotation_about_Xecef)
elif rotation_about_Yscb[2, 0] >= 0 and rotation_about_Yscb[2, 1] < 0 and rotation_about_Yscb[2, 2] < 0:
   polar_Zscb_for_rotation_about_Xecef = -(180 - polar_Zscb_for_rotation_about_Xecef)
elif rotation_about_Yscb[2, 0] >= 0 and rotation_about_Yscb[2, 1] < 0 and rotation_about_Yscb[2, 2] >= 0:
   polar_Zscb_for_rotation_about_Xecef = -(polar_Zscb_for_rotation_about_Xecef)
# Xecef negative
elif rotation_about_Yscb[2, 0] < 0 and rotation_about_Yscb[2, 1] >= 0 and rotation_about_Yscb[2, 2] >= 0:
   polar_Zscb_for_rotation_about_Xecef = polar_Zscb_for_rotation_about_Xecef
elif rotation_about_Yscb[2, 0] < 0 and rotation_about_Yscb[2, 1] >= 0 and rotation_about_Yscb[2, 2] < 0:
   polar_Zscb_for_rotation_about_Xecef = (180 - polar_Zscb_for_rotation_about_Xecef)
elif rotation_about_Yscb[2, 0] < 0 and rotation_about_Yscb[2, 1] < 0 and rotation_about_Yscb[2, 2] < 0:
   polar_Zscb_for_rotation_about_Xecef = -(180 - polar_Zscb_for_rotation_about_Xecef)
elif rotation_about_Yscb[2, 0] < 0 and rotation_about_Yscb[2, 1] < 0 and rotation_about_Yscb[2, 2] >= 0:
   polar_Zscb_for_rotation_about_Xecef = - (polar_Zscb_for_rotation_about_Xecef)

polar_Zscb_for_rotation_about_Xecef = polar_Zscb_for_rotation_about_Xecef * (mt.pi/180)

dcm_rotation_abt_Xscb = np.array([[1, 0, 0],
                                  [0, mt.cos(polar_Zscb_for_rotation_about_Xecef), -mt.sin(polar_Zscb_for_rotation_about_Xecef)],
                                  [0, mt.sin(polar_Zscb_for_rotation_about_Xecef), mt.cos(polar_Zscb_for_rotation_about_Xecef)]])
rotation_about_Xscb = np.array(np.matmul(dcm_rotation_abt_Xscb, rotation_about_Yscb.T)).T
ax.quiver(scb_ref_frame_start[:, 0], scb_ref_frame_start[:, 1], scb_ref_frame_start[:, 2], 
          rotation_about_Xscb[:, 0], rotation_about_Xscb[:, 1], rotation_about_Xscb[:, 2],
          arrow_length_ratio=0.4, linestyle="-",
          color=['red', 'blue', 'green'])
ax.text(scb_origin[0] + rotation_about_Xscb[0, 0], 
        scb_origin[1] + rotation_about_Xscb[0, 1], 
        scb_origin[2] + rotation_about_Xscb[0, 2], 
        'Xscb_R', color='red')
ax.text(scb_origin[0] + rotation_about_Xscb[1, 0], 
        scb_origin[1] + rotation_about_Xscb[1, 1], 
        scb_origin[2] + rotation_about_Xscb[1, 2], 
        'Yscb_R', color='blue')
ax.text(scb_origin[0] + rotation_about_Xscb[2, 0], 
        scb_origin[1] + rotation_about_Xscb[2, 1], 
        scb_origin[2] + rotation_about_Xscb[2, 2], 
        'Zscb_R', color='green')

print(" ")
print(f"azim_rotation_of_Xscb: {azim_rotation_of_Xscb * (180/mt.pi)}")
print(f"elev_rotation_of_Xscb: {elev_rotation_of_Xscb * (180/mt.pi)}")
print(f"polar_Zscb_for_rotation_about_Xecef: {polar_Zscb_for_rotation_about_Xecef * (180/mt.pi)}")

'''display chart'''
plt.show()





