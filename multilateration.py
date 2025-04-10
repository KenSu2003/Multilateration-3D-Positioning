from multilat_lib import *
import csv
import numpy as np


# Load anchor positions and distances
anchors = []
distances = []
sphere_actual = None

with open("tag_data.csv", "r") as file:
    reader = csv.DictReader(file)
    for row in reader:
        if row["Point"] != "Sphere":
            anchors.append([float(row["X"]), float(row["Y"]), float(row["Z"])])
            distances.append(float(row["Distance to Sphere"]))
        else:
            sphere_actual = [float(row["X"]), float(row["Y"]), float(row["Z"])]

anchors = np.array(anchors)
distances = np.array(distances)

estimated_position = brute_force(anchors, distances)
print(f"Estimated Position: {estimated_position}\n")
if sphere_actual:
    error = np.linalg.norm(np.array(sphere_actual) - estimated_position)
    print("Actual Coordinates of Sphere:", sphere_actual)
    print("Estimation Error (Euclidean distance):", error)
