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

'''Cartesian to Spherical co-ordinate conversion'''
def cartesian_to_spherical(x_origin, y_origin, z_origin, x_pt, y_pt, z_pt):
  import math as mt
  x = abs(x_pt - x_origin)
  y = abs(y_pt - y_origin)
  z = abs(z_pt - z_origin)
  length = mt.sqrt(x**2 + y**2 + z**2)
  azimuth = mt.atan2(y/(mt.sqrt(x**2 + y**2)), x/(mt.sqrt(x**2 + y**2)))
  polar = mt.atan2(mt.sqrt(x**2 + y**2)/mt.sqrt(x**2 + y**2 + z**2), z/mt.sqrt(x**2 + y**2 + z**2))
  return length, azimuth, polar

