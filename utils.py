"""Contains all the custom functions, classes, and onjects for this project"""

'''Spherical to Cartesian co-ordinate conversion'''
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
print(spherical_to_cartesian(1, 34.0467, -118.5464)[0])
print(spherical_to_cartesian(1, 34.0467, -118.5464)[1])
print(spherical_to_cartesian(1, 34.0467, -118.5464)[2])

'''Cartesian to Spherical co-ordinate conversion'''
def cartesian_to_spherical(x, y, z):
  import math as mt
  length = mt.sqrt(x**2 + y**2 + z**2)
  azimuth = mt.atan2(y/(mt.sqrt(x**2 + y**2)), x/(mt.sqrt(x**2 + y**2))) * 180/mt.pi
  elevation = 90 - (mt.atan2(mt.sqrt(x**2 + y**2)/mt.sqrt(x**2 + y**2 + z**2), z/mt.sqrt(x**2 + y**2 + z**2)) * 180/mt.pi)
  return length, elevation, azimuth

cartesian_to_spherical(1, 34.0467, -118.5464) # 34.0467°N, 118.5464°W
print(cartesian_to_spherical(-0.3959544967436572, -0.7278511977670105, 0.5598684402764683)[0])
print(cartesian_to_spherical(-0.3959544967436572, -0.7278511977670105, 0.5598684402764683)[1])
print(cartesian_to_spherical(-0.3959544967436572, -0.7278511977670105, 0.5598684402764683)[2])



