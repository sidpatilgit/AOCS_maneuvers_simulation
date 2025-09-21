import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Create a 3D figure
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Step 1: Create a sphere (Earth) with blue color for oceans
u = np.linspace(0, 2 * np.pi, 30)  # Azimuthal angles (longitude)
v = np.linspace(0, np.pi, 30)      # Polar angles (colatitude)
x = np.outer(np.cos(u), np.sin(v))  # x = cos(lon) * sin(lat)
y = np.outer(np.sin(u), np.sin(v))  # y = sin(lon) * sin(lat)
z = np.outer(np.ones(np.size(u)), np.cos(v))  # z = cos(lat)

# Plot the sphere with blue color (oceans)
ax.plot_surface(x, y, z, color='blue', alpha=0.5, rstride=1, cstride=1)

# Step 2: Add continents in green
land = cfeature.LAND.with_scale('110m')  # Low-resolution land geometries
for geom in land.geometries():
    if geom.geom_type == 'Polygon':
        coords = np.array(geom.exterior.coords)
    elif geom.geom_type == 'MultiPolygon':
        for poly in geom:
            coords = np.array(poly.exterior.coords)
            # Convert lon/lat to 3D Cartesian coordinates
            lon = coords[:, 0] * np.pi / 180  # Longitude to radians
            lat = (90 - coords[:, 1]) * np.pi / 180  # Latitude to colatitude
            x_poly = np.cos(lon) * np.sin(lat)
            y_poly = np.sin(lon) * np.sin(lat)
            z_poly = np.cos(lat)
            # Create 3D polygon vertices
            verts = [list(zip(x_poly, y_poly, z_poly))]
            # Add filled polygon (green land)
            ax.add_collection3d(plt.Poly3DCollection(verts, facecolor='green', edgecolor='black', alpha=0.7))
    else:
        continue

# Step 3: Customize axes
ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
ax.set_xlim([-1.5, 1.5])
ax.set_ylim([-1.5, 1.5])
ax.set_zlim([-1.5, 1.5])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Show the plot
plt.title('3D Globe: Green Land, Blue Ocean')
plt.show()