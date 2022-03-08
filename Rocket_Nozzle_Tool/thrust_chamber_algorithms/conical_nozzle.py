import numpy as np


def conical_nozzle_with_ratio(omega, s, R, epsilon, n):
    '''calculate the coordinates of points, which build a conical nozzle.

       Args:
           omega: angle of the nozzle.
           s: radius of throat.
           R: radius of rounding arc at throat.
           epsilon: nozzle ratio (exit area / throat area).
           n: number of points on the nozzle curve.

       Returns:
           **x**; x-coordinates.
           **y**; y-coordinates.
           '''
    # number of points on rounding arc at throat.
    n_omega = 10
    # the central angle between 2 adjacent points on rounding arc.
    d_omega = omega / n_omega

    # coordinate of the point at throat.
    x_throat = 0
    y_throat = s

    # set list for coordinates of the points,
    # which build the rounding arc.
    x_rounding = list()
    y_rounding = list()

    # calculate the coordinates for every point on the rounding arc.
    for i in range(n_omega):
        x_rounding.append(x_throat + R * np.sin(i * d_omega))
        y_rounding.append(y_throat + R * (1 - np.cos(i * d_omega)))

    L = (np.sqrt(s * epsilon) - s) / np.tan(omega)

    # calculate the length of the conical part of the nozzle.
    L_curve = L - R * np.sin(omega)

    # get the distance between 2 adjacent points on conical part.
    d_L = L_curve / n

    # set list for coordinates of the points,
    # which build the conical part of the nozzle.
    x_curve = list()
    y_curve = list()

    # the coordinate of the transitional point between rounding arc
    # and conical part.
    x_transition = R * np.sin(omega)
    y_transition = s + R * (1 - np.cos(omega))

    # calculate the coordinate of every point on the conical part.
    for i in range(n + 1):
        x_curve.append(x_transition + i * d_L)
        y_curve.append(y_transition + i * d_L * np.tan(omega))

    # assemble the coordinates of all points,
    # which build the conical nozzle.
    x = x_rounding + x_curve
    y = y_rounding + y_curve

    return x, y