import csv
from sympy import Point, Matrix, sqrt, Point3D, Plane
from math import sqrt as numeric_sqrt

O = Point(0, 0)

# ———————————————————————— Get Data from Blender ————————————————————————

# Read data from CSV
data_file = "tag_data.csv"
points = {}

with open(data_file, "r") as file:
    reader = csv.DictReader(file)
    for row in reader:
        if row["Point"] != "Sphere":
            points[row["Point"]] = {
                "location": Point(float(row["X"]), float(row["Y"]), float(row["Z"])),
                "distance": float(row["Distance to Sphere"])
            }
        else:
            sphere = Point(float(row["X"]), float(row["Y"]), float(row["Z"]))

# Set Data
A, B, C = points["CubeA"]["location"], points["CubeB"]["location"], points["CubeC"]["location"]
HA, HB, HC = points["CubeA"]["distance"], points["CubeB"]["distance"], points["CubeC"]["distance"]

# Define the three base station points as 3D points
A3D = Point3D(A.x, A.y, A.z)
B3D = Point3D(B.x, B.y, B.z)
C3D = Point3D(C.x, C.y, C.z)

# Create a plane using the three base stations
base_plane = Plane(A3D, B3D, C3D)

# ———————————————————————— Height Calculation ————————————————————————

AB = A.distance(B)
AC = A.distance(C)
BC = B.distance(C)

p = (AB + AC + BC) / 2
S_ABC = sqrt(p * (p - AB) * (p - AC) * (p - BC))  # Area of the base triangle

a = HA
b = HB
c = HC
l = B.distance(C)
m = A.distance(C)
n = A.distance(B)

V_HABC_matrix = Matrix([
    [a**2, (a**2 + b**2 - n**2) / 2, (a**2 + c**2 - m**2) / 2],
    [(a**2 + b**2 - n**2) / 2, b**2, (b**2 + c**2 - l**2) / 2],
    [(a**2 + c**2 - m**2) / 2, (b**2 + c**2 - l**2) / 2, c**2]
])

V_HABC_squared = V_HABC_matrix.det() / 36


V_HABC = numeric_sqrt(abs(float(V_HABC_squared)))

# Calculate HT
HT = 3 * (V_HABC / float(S_ABC))

# Print the results
print("Semi-perimeter of triangle (p):", p)
print("Area of triangle (S_ABC):", float(S_ABC))
print("Volume of tetrahedron (V_HABC):", V_HABC)
print("Height (HT):", HT)


# ———————————————————————— Check Height ————————————————————————
print(f"Check = {sphere.z - HT}")

# Find signed distance from sphere to plane
signed_distance = base_plane.distance(sphere) 

# If the sphere is above, signed_distance will be positive; if below, it will be negative
if signed_distance > 0:
    print(f"The tag is ABOVE the base station plane by {signed_distance:.3f} units.")
elif signed_distance < 0:
    print(f"The tag is BELOW the base station plane by {abs(signed_distance):.3f} units.")
else:
    print("The tag is ON the same plane as the base stations.")