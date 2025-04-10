import csv
import numpy as np

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


# def trilateration_closed_form(anchors, distances):
#     '''
#     Larsson 2022 — 4.3 + 4.7: Linear Trilateration with Double Compaction Matrix (M)
#     '''

#     if len(anchors) != 4:
#         raise ValueError("This closed-form solver requires exactly 4 anchors!")

#     assert anchors.shape == (4, 3), # This solver needs exactly 4 anchors in 3D
#     assert distances.shape == (4,), # Need 4 distances

#     # Shift an achor as the origin to be the reference point
#     origin = anchors[0]                 # pick 1st anchor as origin (reference point)
#     distances = distances**2            # work with squared distances

#     # Shift the anchors to the origin, and calculate the distance matrix
#     Anchor_offset = anchors - origin   

#     # Find the Matrix of anchor offsets (local coordinate differences)
#     A = Anchor_offset[1:]                           # anchors 2-4

#     # 
#     b = 0.5 * (
#         distances[0]                            # anchor 1 (origin) squared distance
#         - distances[1:]                         # anchors 2-4 squared distances
#         + np.sum(D[1:]**2, axis=1) # anchors 2-4 squared distances
#     )

#     # Solve and shift back
#     estimated_shifted = np.linalg.solve(A, b)   # Solve the Linear System
#     estimated = estimated_shifted + origin      # Shift back to global coordinates (x,y,z)
#     return estimated




def trilateration_minimum_squared(anchors, distances):
    '''
    Larsson 2022 — Linear Trilateration with 4+ Anchors (Least Squares Form)

    Solves A @ t = b for the tag position t, using anchor positions (a_i) and distances (d_i).
    Works with 4 or more anchors in 3D space.

    - t      = unknown tag position
    - a_i    = anchor positions (shape: [n, 3])
    - d_i    = measured distances (shape: [n])
    '''

    assert anchors.shape[1] == 3, "Anchors must be 3D"
    assert anchors.shape[0] >= 4, "Need at least 4 anchors"
    assert anchors.shape[0] == distances.shape[0], "Anchor and distance count mismatch"

    # Step 1: Shift all anchors relative to the first (reference)
    origin = anchors[0]                 # set the first anchor as the origin, reference point
    D = anchors - origin                # shift the anchors relative to the origin
    distances_squared = distances**2    # find the euclidean distance squared

    # Step 2: Build matrix A (anchor offset vectors)
    '''
    A = ai - a1     # Matrix of anchor offsets (local coordinate differences)
        a_i = (x_i, y_i, z_i) , a1 = (x1, y1, z1)
    A = [
        [x2 - x1, y2 - y1, z2 - z1],
        [x3 - x1, y3 - y1, z3 - z1],
        [x4 - x1, y4 - y1, z4 - z1]
    ]
    '''
    A = D[1:]  # shape: (n-1, 3)
    

    # Step 3: Build b vector using linearized distance equation
    b = 0.5 * (
        distances_squared[0] - distances_squared[1:] + np.sum(D[1:]**2, axis=1)
    )  # shape: (n-1,)

    # Step 4: Solve A @ t = b in least squares sense
    t_shifted, *_ = np.linalg.lstsq(A, b, rcond=None)
    t = t_shifted + origin

    return t
    


# --- Solve position ---


# estimated_position = trilateration_minimum_squared(anchors, distances)
estimated_position = trilateration_minimum_squared(anchors, distances)
print("Estimated Coordinates of Sphere:", estimated_position)

if sphere_actual:
    error = np.linalg.norm(np.array(sphere_actual) - estimated_position)
    print("Actual Coordinates of Sphere:", sphere_actual)
    print("Estimation Error:", error)