import bpy
import math

# Function to calculate the Euclidean distance
def calculate_distance(point1, point2):
    return math.sqrt(
        (point2[0] - point1[0]) ** 2 +
        (point2[1] - point1[1]) ** 2 +
        (point2[2] - point1[2]) ** 2
    )

sphere = bpy.data.objects["Sphere"].location
cubeA = bpy.data.objects["Cube"].location
cubeB = bpy.data.objects["Cube.001"].location
cubeC = bpy.data.objects["Cube.002"].location

print(f"Sphere Location: {sphere} \n CubeA Location: {cubeA} \n CubeB Location: {cubeB} \n CubeC Location: {cubeC}")

AB = calculate_distance(cubeA,cubeB)
AC = calculate_distance(cubeA,cubeC)
BC = calculate_distance(cubeB,cubeC)
# print(f"AB = {AB}")

HA = calculate_distance(sphere, cubeA)
HB = calculate_distance(sphere, cubeB)
HC = calculate_distance(sphere, cubeC)
# print(f"HA = {HA}, HB = {HB}, HC = {HC}")