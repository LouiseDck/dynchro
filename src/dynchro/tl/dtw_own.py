import numpy as np
import scipy.spatial


def skip_dtw(x, y, distance="correlation"):
    r, c = len(x), len(y)
    if distance == "euclidean":
        cost = scipy.spatial.distance.cdist(x, y, "euclidean")
    else:
        cost = scipy.spatial.distance.cdist(x, y, "correlation")
    distances = np.zeros((r + 1, c + 1))
    # distances[1:, 1:] = cost
    distances[0, 1:], distances[1:, 0] = np.inf, np.inf
    distances[0, 0] = 0

    for i in range(1, r + 1):
        for j in range(1, c + 1):
            options = [distances[i - 1, j - 1], distances[i - 1, j], distances[i, j - 1]]
            distances[i, j] = min(options) + cost[i - 1, j - 1]

    return distances[-1, -1], cost, distances  # , path


def skip_traceback(D):
    i, j = np.array(D.shape) - 1

    # j = np.argmin(D[:,-1])

    p, q = [i], [j]
    while (i > 0) or (j > 0):
        tb = np.argmin((2 * D[i - 1, j - 1], D[i, j - 1], D[i - 1, j]))
        if tb == 0:
            i -= 1
            j -= 1
        elif tb == 1:
            j -= 1
        elif tb == 2:  # (tb == 2):
            i -= 1
        p.insert(0, i)
        q.insert(0, j)
    return np.array(p), np.array(q)


def skip_traceback2(D):
    i, j = np.array(D.shape) - 1

    # j = np.argmin(D[:,-1])

    p, q = [i], [j]
    while (i > 0) or (j > 0):
        tb = np.argmin((2 * D[i - 1, j - 1], D[i, j - 1], D[i - 1, j], 2 * D[i - 2, j - 2], D[i, j - 2], D[i - 2, j]))
        if tb == 0:
            i -= 1
            j -= 1
        elif tb == 1:
            j -= 1
        elif tb == 2:  # (tb == 2):
            i -= 1
        elif tb == 3:
            i -= 2
            j -= 2
        elif tb == 4:
            j -= 2
        elif tb == 5:
            i -= 2
        p.insert(0, i)
        q.insert(0, j)
    return np.array(p), np.array(q)
