from scipy.optimize import minimize
from sympy import Point3D, Plane, Matrix, sqrt as symbolic_sqrt
import numpy as np
import math



def brute_force(anchor_coords=None, x=None, y=None, distances=None):
    '''
    Hybrid brute force:
    - Full 3D: estimate [x, y, z] when x and y are None
    - Z-only: estimate z when x and y are provided
    '''

    if x is None and y is None:
        # Full 3D optimization
        def objective(q):
            return np.sum((np.linalg.norm(anchor_coords - q, axis=1) - distances) ** 2)

        initial_guess = np.mean(anchor_coords, axis=0)
        res = minimize(objective, initial_guess)
        return res.x  # full [x, y, z]

    else:
        # Estimate Z only, given fixed x and y
        def objective_z(z, anchor_coords, distances):
            q = np.array([x, y, z[0]])
            return np.sum((np.linalg.norm(anchor_coords - q, axis=1) - distances) ** 2)

        res = minimize(objective_z, [0.0], args=(anchor_coords, distances))
        return res.x[0]
        

def multilateration_minimum_squared(anchors, distances):
    '''
    Larsson 2022 — Linear Trilateration with 4+ Anchors (Least Squares Form)
    Solves A @ t = b for the tag position t.
    '''

    if anchors.shape[1] != 3:
        raise ValueError("Anchors must be 3D (shape: [n, 3])")
    if anchors.shape[0] < 4:
        raise ValueError("At least 4 anchors are required for multilateration.")
    if anchors.shape[0] != distances.shape[0]:
        raise ValueError(f"Anchor and distance count mismatch: {anchors.shape[0]} anchors, {distances.shape[0]} distances")

    # Step 1: Shift all anchors relative to the first (reference)
    origin = anchors[0]
    D = anchors - origin
    distances_squared = distances**2

    # Step 2: Build matrix A
    A = D[1:]

    # Step 3: Build vector b
    b = 0.5 * (
        distances_squared[0] - distances_squared[1:] + np.sum(D[1:]**2, axis=1)
    )

    # Step 4: Solve the system
    t_shifted, *_ = np.linalg.lstsq(A, b, rcond=None)
    t = t_shifted + origin

    return t


def multilateration_closed_form(anchors, distances):
    '''
    Requires exactly 4 anchors
    Larsson 2022 — Closed-form Trilateration
    '''
    if anchors.shape != (4, 3):
        raise ValueError("This solver needs exactly 4 anchors in 3D (shape must be [4, 3])")
    if distances.shape != (4,):
        raise ValueError("This solver needs exactly 4 distances (shape must be [4,])")

    origin = anchors[0]
    distances = distances**2
    Anchor_offset = anchors - origin
    A = Anchor_offset[1:]

    b = 0.5 * (
        distances[0] - distances[1:] + np.sum(Anchor_offset[1:]**2, axis=1)
    )

    estimated_shifted = np.linalg.solve(A, b)
    estimated = estimated_shifted + origin
    return estimated

def trilateration(anchors, distances, ignore=False):
    '''
    Estimate the 3D position of a tag using geometric trilateration with exactly 3 UWB anchors.

    This implementation is based on the algorithm proposed in the paper
    "A Precise 3D Positioning Approach Based on UWB with Reduced Base Stations" (Xu et al., 2021),
    which reduces the number of required base stations from four to three by projecting the tag height
    from a known triangle plane formed by the anchors.
    '''

    if anchors.shape != (3, 3):
        raise ValueError(f"This solver needs exactly 3 anchors in 3D (shape must be [3, 3]).\n There are exactly {len(anchors)} anchors")
    if distances.shape != (3,):
        raise ValueError(f"This solver needs exactly 3 distances (shape must be [3,]).\n There are exactly {len(distances)} distances.")

    # Enforce anchor height uniformity
    z_coords = anchors[:, 2]
    if not ignore and not np.allclose(z_coords, z_coords[0], atol=1e-3):
        raise ValueError("Anchors must be on the same horizontal plane (same Z coordinate).")

    # Unpack anchors and distances
    A, B, C = anchors
    HA, HB, HC = distances

    # Convert anchor points to sympy Point3D
    A3D, B3D, C3D = Point3D(*A), Point3D(*B), Point3D(*C)

    # Define the triangle plane
    base_plane = Plane(A3D, B3D, C3D)

    # Triangle side lengths
    AB = A3D.distance(B3D)
    AC = A3D.distance(C3D)
    BC = B3D.distance(C3D)

    # Semi-perimeter and area
    p = (AB + AC + BC) / 2
    S_ABC = symbolic_sqrt(p * (p - AB) * (p - AC) * (p - BC))

    # Convert sympy distances to float
    n = float(AB)
    m = float(AC)
    l = float(BC)

    # Cayley-Menger determinant-based volume
    V_HABC_matrix = Matrix([
        [HA**2, (HA**2 + HB**2 - n**2) / 2, (HA**2 + HC**2 - m**2) / 2],
        [(HA**2 + HB**2 - n**2) / 2, HB**2, (HB**2 + HC**2 - l**2) / 2],
        [(HA**2 + HC**2 - m**2) / 2, (HB**2 + HC**2 - l**2) / 2, HC**2]
    ])
    
    V_HABC_squared = V_HABC_matrix.det() / 36
    V_HABC = np.sqrt(abs(float(V_HABC_squared)))

    # Height from tag to triangle
    HT = 3 * V_HABC / float(S_ABC)

    # Estimate tag position: centroid + height * normal_vector
    base_centroid = np.mean(anchors, axis=0)
    normal_vec = np.array(base_plane.normal_vector).astype(np.float64)
    normal_vec /= np.linalg.norm(normal_vec)

    estimated_position = base_centroid + HT * normal_vec

    return estimated_position



def local_to_geo(local_coords, origin_geo):
    """
    Converts local XYZ (meters) back to geodetic coordinates (lat, lon, alt).
    
    Parameters:
    -----------
    local_coords : np.ndarray of shape (3,)
        The local ENU coordinates in meters.
    origin_geo : tuple
        The geodetic origin (lat0, lon0, alt0)

    Returns:
    --------
    lat, lon, alt : float
        Estimated geographic coordinates.
    """
    lat0, lon0, alt0 = origin_geo

    # Latitude and longitude scales
    lat_scale = 111_000
    lon_scale = (40_075_000 / 360) * math.cos(math.radians(lat0))

    dx, dy, dz = local_coords
    lat = lat0 + (dy / lat_scale)
    lon = lon0 + (dx / lon_scale)
    alt = alt0 + dz

    return lat, lon, alt


# Convert lat/lon/alt to local XYZ
def geo_to_local_xyz(anchor_coords):
    """
    Converts geodetic coordinates (lat, lon, alt) to local Cartesian (x, y, z) in meters.
    Uses the first anchor as the local origin (0, 0, 0).
    
    Parameters:
    -----------
    anchor_coords : np.ndarray or list of tuples
        Each row is (latitude, longitude, altitude)

    Returns:
    --------
    local_coords : np.ndarray of shape (n, 3)
        Local ENU (East, North, Up) coordinates in meters
    """
    lat0, lon0, alt0 = anchor_coords[0]

    # Latitude and Longitude scale in meters per degree
    lat_scale = 111_000  # approx constant
    lon_scale = (40_075_000 / 360) * math.cos(math.radians(lat0))  # varies with latitude

    local_coords = []
    for lat, lon, alt in anchor_coords:
        dx = (lon - lon0) * lon_scale
        dy = (lat - lat0) * lat_scale
        dz = alt - alt0
        local_coords.append([dx, dy, dz])

    return np.array(local_coords)