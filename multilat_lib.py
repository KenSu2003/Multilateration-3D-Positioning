from scipy.optimize import minimize
from sympy import Point3D, Plane, Matrix, sqrt as symbolic_sqrt
import numpy as np



def brute_force(anchors, distances):
    '''
    Literally just finding the best fit. 
    '''

    # Objective function: minimize squared difference between measured and calculated distances
    def objective(q):
        return np.sum((np.linalg.norm(anchors - q, axis=1) - distances) ** 2)

    # Initial guess: average of anchor positions
    initial_guess = np.mean(anchors, axis=0)

    # Optimization
    res = minimize(objective, initial_guess)

    # Output
    estimated_position = res.x
    print("Estimated Coordinates of Sphere:", estimated_position)

    return estimated_position


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

    Parameters:
    -----------
    anchors : np.ndarray of shape (3, 3)
        The coordinates of the three base stations (anchors), where each anchor is a 3D point [x, y, z].
        All three anchors must lie in the same horizontal plane (i.e., have approximately the same Z value).
    
    distances : np.ndarray of shape (3,)
        The distances from each anchor to the target tag (HA, HB, HC), typically derived from TOF measurements.
    
    ignore : boolean
        Ignores whether the anchors are on the same height.

    Returns:
    --------
    estimated_position : np.ndarray of shape (3,)
        The estimated 3D coordinates of the tag [x, y, z], calculated by projecting from the base triangle
        using geometric volume and plane normal.

    Raises:
    -------
    ValueError:
        - If the shape of `anchors` or `distances` is incorrect.
        - If the anchors are not coplanar in the XY plane (i.e., Z-coordinates are not approximately equal).

    Notes:
    ------
    - Assumes a triangle base defined by three anchors lying flat on the same Z-plane.
    - The method calculates the volume of the tetrahedron formed by the three anchors and the tag,
      then uses that volume and the triangle area to determine the height of the tag above the plane.
    - Accuracy significantly degrades if anchors are not level.

    Reference:
    ----------
    Zhiqiang Xu et al., "A Precise 3D Positioning Approach Based on UWB with Reduced Base Stations", 2021.
    '''

    if anchors.shape != (3, 3):
        raise ValueError("This solver needs exactly 3 anchors in 3D (shape must be [3, 3])")
    if distances.shape != (3,):
        raise ValueError("This solver needs exactly 3 distances (shape must be [3,])")

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
