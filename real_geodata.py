import multilat_lib
import numpy as np
from multilat_lib import *

# Real-world anchor GPS data
anchors_geo = [
    (42.39379, -72.52856, 71.13),
    (42.39381, -72.52861, 71.12),
    (42.39383, -72.52862, 71.09)
]

tag_geo = (42.39382, -72.52859, 71.11)

distances = np.array([4.15, 1.99, 2.71])  # in meters
anchors = geo_to_local_xyz(anchors_geo)

estimated_position = brute_force(anchors, distances)

lat_est, lon_est, alt_est = local_to_geo(estimated_position, anchors_geo[0])

# Project actual tag position to base plane
vertical_error = tag_geo[2] - alt_est

# Output
print(f"\nEstimated Geo Coordinates:\nLatitude: {lat_est}\nLongitude: {lon_est}\nAltitude: {alt_est}")
print("Actual Tag Position:", tag_geo)
print("Vertical Error (Tag Z - Plane Z):", vertical_error)




