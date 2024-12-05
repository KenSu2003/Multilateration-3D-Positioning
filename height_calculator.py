from sympy import Point, Matrix, sqrt
from math import sqrt as numeric_sqrt
# import tag_response as tag
# from tag_response import HA, HB, HC, cubeA, cubeB, cubeC, sphere

O = Point(0, 0)

# A, B, C = tag.cubeA, tag.cubeB, tag.cubeC
# HA, HB, HC = tag.HA, tag.HB, tag.HC

A = Point(10.0000, 1.0000, 0)
B = Point(1.0000, 10.0000, 0)
C = Point(1.0000, 1.0000, 0)

HA = 8.433780570691637
HB = 7.589317959056357
HC = 7.633360750059722



# ———————————————————————— DO NOT EDIT CODE BELOW ————————————————————————

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


# Check
