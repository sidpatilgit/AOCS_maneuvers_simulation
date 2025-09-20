""" Create a rotating Earth; create the exact change in tilt, wobble, angular velocity based on how Earth is currently behaving"""

'''Create a rotating, titled Earth'''

import numpy as np #generating data
import matplotlib.pyplot as plt #for plotting
import matplotlib.ticker as mtk
from mpl_toolkits import mplot3d #for 3D plots
import math as mt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# Plot all points on Earth
ax = plt.figure(figsize=(10, 10)).add_subplot(111, projection="3d")
radius_earth = 1
latitude_data = np.linspace(-np.pi/2, np.pi/2, 30) #(pi/2)-polar angle (or elevation angle) in radians
longitude_data = np.linspace(-180*(mt.pi/180), 180*(mt.pi/180), 30) #azimuth in radians

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

ax.plot_wireframe(x_data, y_data, z_data, alpha=0.002, color='gray', rstride=1, cstride=1)

def plot_continents_on_sphere(ax):
    # Use Cartopy's PlateCarree projection for geographic data
    proj = ccrs.PlateCarree()
    
    # Get continent boundaries from Cartopy
    continents = cfeature.NaturalEarthFeature(
        category='physical', name='land', scale='110m', facecolor='none', edgecolor='black'
    )
    
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
                ax.plot(x, y, z, color='black', linewidth=1)
plot_continents_on_sphere(ax)
plt.show()
