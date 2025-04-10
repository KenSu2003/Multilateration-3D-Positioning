from scipy.optimize import minimize
import numpy as np

def brute_force(anchors, distances):

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
