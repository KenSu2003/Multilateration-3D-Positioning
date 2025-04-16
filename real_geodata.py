import numpy as np
from multilat_lib import geo_to_local_xyz, local_to_geo
from multilat_lib import *

# Anchor class + 2D trilateration as before
class Anchor:
    def __init__(self, id, x, y):
        self.id = id
        self.x_coord = x
        self.y_coord = y

def trilaterate(anchors, distances):
    x1, x2, x3 = anchors[0].x_coord, anchors[1].x_coord, anchors[2].x_coord
    y1, y2, y3 = anchors[0].y_coord, anchors[1].y_coord, anchors[2].y_coord
    r1, r2, r3 = distances

    A = -2*x1 + 2*x2
    B = -2*y1 + 2*y2
    C = r1**2 - r2**2 - x1**2 + x2**2 - y1**2 + y2**2
    D = -2*x2 + 2*x3
    E = -2*y2 + 2*y3
    F = r2**2 - r3**2 - x2**2 + x3**2 - y2**2 + y3**2

    X = (C*E - F*B) / (E*A - B*D)
    Y = (C*D - A*F) / (B*D - A*E)
    return (X, Y)


def error_analysis(estimated,actual):
    return (estimated-actual)/actual

# --- Input Coordinates ---
anchors_geo = [
    (42.39366207403457, -72.5296460245343, 71.13234),
    (42.39388705549964, -72.52889805962968, 71.12234),
    (42.39356330142929, -72.52869249311617, 71.01239)
]

tag_geo_true = (42.39363871500632, -72.52916719231453, 71.143205)
origin_geo = min(anchors_geo, key=lambda x: (x[0], x[1]))

# Convert to local coordinates
anchor_coords = geo_to_local_xyz(anchors_geo, origin_geo)
tag_local = geo_to_local_xyz([tag_geo_true], origin_geo)[0]
anchor_list = [Anchor(i, x, y) for i, (x, y, z) in enumerate(anchor_coords)]

# Simulate distances
def simulate_distances_3d(anchor_coords, tag, noise_std=0.05):
    dists = np.linalg.norm(anchor_coords - tag, axis=1)
    return dists + np.random.normal(0, noise_std, size=dists.shape)

# distances = simulate_distances_3d(anchor_coords, tag_local)
distances = simulate_distances_3d(anchor_coords, tag_local, noise_std=0.0)

# 1. Estimate XY using 2D trilateration
estimated_x, estimated_y = trilaterate(anchor_list, distances)

# 2. Estimate Z using your hybrid brute_force function
estimated_z = brute_force(anchor_coords=anchor_coords, x=estimated_x, y=estimated_y, distances=distances)
# estimated_z = brute_force(anchor_coords=anchor_coords, distances=distances)[2]

# Combine into full 3D estimate
estimated_local = np.array([estimated_x, estimated_y, estimated_z])
lat_est, lon_est, alt_est = local_to_geo(estimated_local, origin_geo)

# Output
print("\n=== Hybrid Trilateration (XY: 2D + Z: Brute Force) ===")
print(f"Estimated Local: x={estimated_x:} m, y={estimated_y:} m, z={estimated_z:} m")

print("\nEstimated Tag GPS:")
print(f"Latitude:  {lat_est}")
print(f"Longitude: {lon_est}")
print(f"Altitude:  {alt_est} m")


lat_actual, long_actual, alt_actual = tag_geo_true
print("\nActual Tag GPS:")
print(f"Latitude:  {tag_geo_true[0]}")
print(f"Longitude: {tag_geo_true[1]}")
print(f"Altitude:  {tag_geo_true[2]} m")


print("\n=== Error Analysis ===")

print("\nLocal XYZ Comparison")
print(f"Estimated Local: x={estimated_local[0]:.3f}, y={estimated_local[1]:.3f}, z={estimated_local[2]:.3f}")
print(f"Actual Local:    x={tag_local[0]:.3f}, y={tag_local[1]:.3f}, z={tag_local[2]:.3f}")

print("\nError in Local Coordinates (Estimated - Actual):")
print(f"Δx = {estimated_local[0] - tag_local[0]:.3f} m")
print(f"Δy = {estimated_local[1] - tag_local[1]:.3f} m")
print(f"Δz = {estimated_local[2] - tag_local[2]:.3f} m")

print("\nError in Global Coordinates:")
print("Latitude Error:", f"{abs(error_analysis(lat_est,lat_actual))} %")
print("Longitude Error:", f"{abs(error_analysis(lon_est,long_actual))} %")
print("Altitude Error:", f"{abs(error_analysis(alt_est,alt_actual))} m")
