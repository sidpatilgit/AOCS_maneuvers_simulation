'''Import the required libraries'''
import numpy as np
import matplotlib.pyplot as mp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mtk
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import random as rd
from scipy.spatial.transform import Rotation as R
import math as m

formatter = mtk.FormatStrFormatter('%.01f')

## Plot the Earth Centered/Earth Fixed (ECEF) frame --> this reference frame rotates as the Earth rotates, Xecef points towards the prime meridian
# convert lat long to spherical, then azimutal and elevation angles in tooltip
# Set up the figure and 3D axes
fig = mp.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Create a sphere for Earth
u = np.linspace(0, 2 * np.pi, 30)  # Azimuthal angles
v = np.linspace(0, np.pi, 30)      # Polar angles
x = np.outer(np.cos(u), np.sin(v))  # x = cos(u) * sin(v)
y = np.outer(np.sin(u), np.sin(v))  # y = sin(u) * sin(v)
z = np.outer(np.ones(np.size(u)), np.cos(v))  # z = cos(v)

# Plot the wireframe sphere
ax.plot_wireframe(x, y, z, color='black', alpha=0.05, rstride=1, cstride=1)

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
mp.show()