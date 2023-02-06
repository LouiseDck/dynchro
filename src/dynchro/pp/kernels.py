import numpy as np
from numpy import ndarray


def apply_gaussian_kernel(x: ndarray, counts: ndarray) -> ndarray:
    d = np.exp(-1 * (x**2) / (0.05**2))

    # FOR VISUALISATION PURPOSES
    # import matplotlib.pyplot as plt
    # plt.style.use('seaborn-whitegrid')

    # plt.plot(x, color = "red")
    # plt.plot(d, color = "blue")
    # plt.figure()
    # plt.show()

    d = d / sum(d)
    # plt.plot(d, color = "green")

    # plt.show()

    # @ = matrix multiplication (is transposition necessary or could the arguments be switched?)
    return counts.T @ d


def get_gaussian_transition_matrix(x: ndarray, counts: ndarray) -> ndarray:
    d = np.exp(-1 * (x**2) / (0.05**2))
    return d / sum(d)


def apply_epanechnikov_kernel(x: ndarray, counts: ndarray) -> ndarray:
    # Can this be done faster? Vectorized? -> first make everything 0 that is >= 1?
    d = [0.75 * (1 - y**2) if y <= 1 else 0 for y in x]
    # x = np.where(x > 0.005, math.sqrt(0.005), x)
    # d = 0.75 * (0.005 - x ** 2)
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
