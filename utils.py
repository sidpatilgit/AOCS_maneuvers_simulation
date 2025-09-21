"""Contains all the custom functions, classes, and onjects for this project"""

def spherical_to_cartesian(radius, latitude, longitude): # (radius, elevation, azimuthal)
  import numpy as np
  import math as mt
  latitude_radians = latitude * (mt.pi/180)
  longitude_radians = longitude * (mt.pi/180)
  x_cood = radius * mt.cos(latitude_radians) * mt.cos(longitude_radians)
  y_cood = radius * mt.cos(latitude_radians) * mt.sin(longitude_radians)
  z_cood = radius * mt.sin(latitude_radians)
  #print(f"(The x, y, and z co-ordinates of the point ({radius}, {latitude}, {longitude}) are {x_cood}, {y_cood}, {z_cood})")
  return x_cood, y_cood, z_cood

spherical_to_cartesian(1, 34.0467, -118.5464) # 34.0467°N, 118.5464°W
