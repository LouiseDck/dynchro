import numpy as np
from numpy import ndarray


def apply_gaussian_kernel(x: ndarray, counts: ndarray, bw = 0.05) -> ndarray:
    print(bw)
    d = get_gaussian_transition_matrix(x, bw)
    return counts.T @ d


def get_gaussian_transition_matrix(x: ndarray, bw = 0.05) -> ndarray:
    d = np.exp(-1 * (x**2) / (bw**2))
    return d / sum(d)


def apply_epanechnikov_kernel(x: ndarray, counts: ndarray) -> ndarray:
    d = np.where(x <= 1, 0.75 * (1 - x**2), 0)
    d /= sum(d)
    return counts.T @ d


def apply_triangular_kernel(x: ndarray, counts: ndarray) -> ndarray:
    x = np.where(x > 1, 0, x)
    d = 1 - x
    d /= sum(d)
    return counts.T @ d


def apply_boxcar_kernel(x: ndarray, counts: ndarray) -> ndarray:
    d = x
    d = np.where(d <= 0.1, 0.05, d)
    d = np.where(d > 0.1, 0, d)
    d /= sum(d)
    return counts.T @ d
