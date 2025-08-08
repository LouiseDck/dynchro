import numpy as np
import scipy.spatial


def dtw(x, y, distance="euclidean"):
    cost = scipy.spatial.distance.cdist(x, y, distance)
    tol = 1e-15
    cost[abs(cost) < tol] = 0

    r, c = len(x), len(y)

    distances = np.zeros((r + 1, c + 1))
    distances[0, 1:], distances[1:, 0] = np.inf, np.inf
    distances[0, 0] = 0

    for i in range(1, r + 1):
        for j in range(1, c + 1):
            options = [distances[i - 1, j - 1], distances[i - 1, j], distances[i, j - 1]]
            distances[i, j] = min(options) + cost[i - 1, j - 1]

    return distances[-1, -1], cost, distances[1:, 1:]


def traceback_start(D, i, j):
    p, q = [i], [j]
    while (i > 0) or (j > 0):
        tb = np.argmin((D[i - 1, j - 1], D[i, j - 1], D[i - 1, j]))
        if tb == 0:
            i -= 1
            j -= 1
        elif tb == 1:
            j -= 1
        elif tb == 2:
            i -= 1
        p.insert(0, i)
        q.insert(0, j)
    return np.array(p), np.array(q)


def traceback(D):
    i, j = np.array(D.shape) - 1

    return traceback_start(D, i, j)


def traceback_yedge(D):
    i = np.argmin(D[:, -1])
    j = D.shape[1] - 1

    return traceback_start(D, i, j)


def traceback_xedge(D):
    j = np.argmin(D[-1, :])
    i = D.shape[0] - 1

    return traceback_start(D, i, j)



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


def amerced_dtw(x, y, distance="correlation", penalty=0):
    r, c = len(x), len(y)
    if distance == "euclidean":
        cost = scipy.spatial.distance.cdist(x, y, "euclidean")
    else:
        cost = scipy.spatial.distance.cdist(x, y, "correlation")

    tol = 1e-15
    cost[abs(cost) < tol] = 0

    distances = np.zeros((r + 1, c + 1))
    distances[0, 1:], distances[1:, 0] = np.inf, np.inf
    distances[0, 0] = 0

    for i in range(1, r + 1):
        for j in range(1, c + 1):
            opt1 = distances[i - 1, j - 1] + cost[i - 1, j - 1]
            opt2 = distances[i - 1, j] + cost[i - 1, j - 1] + penalty
            opt3 = distances[i, j - 1] + cost[i - 1, j - 1] + penalty
            distances[i, j] = min(opt1, opt2, opt3)
