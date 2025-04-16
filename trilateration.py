import csv
import numpy as np
from sympy import Point3D, Plane
from multilat_lib import trilateration


''' Simulate Using Blender Data '''
# Read data from Blender
data_file = "tag_data.csv"
points = {}

with open(data_file, "r") as file:
    reader = csv.DictReader(file)
    for row in reader:
        label = row["Point"]
        x, y, z = float(row["X"]), float(row["Y"]), float(row["Z"])
        if label != "Sphere":
            points[label] = {
                "location": np.array([x, y, z]),
                "distance": float(row["Distance to Sphere"])
            }
        else:
            sphere = np.array([x, y, z])

# Extract anchor locations and distances
A = points["CubeA"]["location"]
B = points["CubeB"]["location"]
C = points["CubeC"]["location"]
HA = points["CubeA"]["distance"]
HB = points["CubeB"]["distance"]
HC = points["CubeC"]["distance"]

anchors = np.array([A,B,C])
distances = np.array([HA,HB,HC])

# Run trilateration
estimated_position = trilateration(anchors, distances)

# Project actual tag position to base plane
error = sphere[2] - estimated_position[2]

# Output
print("Estimated Tag Position (from trilateration):", estimated_position)
print("Actual Sphere Position (used in simulation):", sphere)
print("Error (Tag Z - Plane Z):", error)
