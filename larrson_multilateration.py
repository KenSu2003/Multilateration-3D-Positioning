import csv
import numpy as np
import multilat_lib.py

# --- Load data from CSV ---
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





# --- Solve position ---


# estimated_position = trilateration_minimum_squared(anchors, distances)
estimated_position = multilateration_minimum_squared(anchors, distances)
print("Estimated Coordinates of Sphere:", estimated_position)

if sphere_actual:
    error = np.linalg.norm(np.array(sphere_actual) - estimated_position)
    print("Actual Coordinates of Sphere:", sphere_actual)
    print("Estimation Error:", error)