import bpy
import math
import csv

# Function to calculate the Euclidean distance
def calculate_distance(point1, point2):
    return math.sqrt(
        (point2[0] - point1[0]) ** 2 +
        (point2[1] - point1[1]) ** 2 +
        (point2[2] - point1[2]) ** 2
    )

# Get locations
sphere = bpy.data.objects["Sphere"].location
cubeA = bpy.data.objects["Cube"].location
cubeB = bpy.data.objects["Cube.001"].location
cubeC = bpy.data.objects["Cube.002"].location
cubeD = bpy.data.objects["Cube.003"].location

# Calculate distances
HA = calculate_distance(sphere, cubeA)
HB = calculate_distance(sphere, cubeB)
HC = calculate_distance(sphere, cubeC)
HD = calculate_distance(sphere, cubeD)

# Save to file
output_file = "/Users/kensu/Desktop/SDP/Multilateration-3D-Positioning/tag_data.csv"

with open(output_file, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Point", "X", "Y", "Z", "Distance to Sphere"])
    writer.writerow(["CubeA", cubeA[0], cubeA[1], cubeA[2], HA])
    writer.writerow(["CubeB", cubeB[0], cubeB[1], cubeB[2], HB])
    writer.writerow(["CubeC", cubeC[0], cubeC[1], cubeC[2], HC])
    writer.writerow(["CubeD", cubeD[0], cubeD[1], cubeD[2], HD])
    writer.writerow(["Sphere", sphere[0], sphere[1], sphere[2], ""])

print(f"Data written to {output_file}")
