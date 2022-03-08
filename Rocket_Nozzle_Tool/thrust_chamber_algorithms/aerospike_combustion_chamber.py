import numpy as np
import math


def inner_chamber_contour(x0, y0, alpha, r, x_end, d):
    '''generate the inner contour of combustion chamber.

       Args:
           x0, y0: coordinate of the first point on inner contour.
           alpha: the total turning angle of aerospike.[rad]
           r: rounding radius at throat.
           x_end: x-coordinate of the last point on inner contour.
           d: distance between two points.

       Returns:
           **x_chamber**; list for x-coordinates of the points on
                      inner contour.
           **y_chamber**; list for y-coordinates of the points on
                      inner contour.
                    '''
    # calculate the angle between horizontal axis and the line, which
    # between circular center and the first point on circular arc.
    b = 0.5 * np.pi - abs(alpha)

    y_r = y0 - r * np.sin(b)  # y coordinate of center of circle.
    x_r = x0 - r * np.cos(b)  # x coordinate of center of circle.

    # list for coordinate of point on rounding.
    r_x = list()
    r_y = list()

    # coordinate for the last point of the rounding arc /
    # the first point of the parallel part.
    x_r_end = x_r
    y_r_end = y_r + r

    r_l = r * alpha  # length of the rounding circular arc.
    n = math.ceil(r_l / d)  # the circular arc be partitioned in n parts.
    a_step = alpha / n  # the angle for every part of circular arc.

    # calculate the coordinate for every point on the circular arc.
    for i in range(n - 1):
        r_x.append(x_r + r * np.cos(b + (i + 1) * a_step))
        r_y.append(y_r + r * np.sin(b + (i + 1) * a_step))

    l_x = list()  # list for x-coordinates of the points on parallel part.
    l_y = list()  # list for y-coordinates of the points on parallel part.
    l_x.append(x_r_end)  # set the first point on the parallel part.
    l_y.append(y_r_end)  # set the first point on the parallel part.

    # set the coordinate for every point on the parallel part.
    i2 = 0
    while l_x[i2] > x_end:
        l_x.append(l_x[i2] - d)
        l_y.append(l_y[i2])
        i2 = i2 + 1

    # assemble all points on inner contour of combustion chamber.
    x_chamber = r_x + l_x
    y_chamber = r_y + l_y

    # reverse the order of all points. Because the points should be
    # outputted with increasing x-coordinate.
    x_chamber.reverse()
    y_chamber.reverse()

    return x_chamber, y_chamber


def outer_chamber_contour(x_exp, y_exp, alpha, r, l1, l2, d):
    '''generate the outer contour of combustion chamber.

       Args:
           x_exp, y_exp: coordinate of expansion point.
           alpha: total turning angle of aerospike./convergence angle[rad]
           r: rounding radius between convergence part and
              parallel part.
           l1: length of parallel part.
           l2: length of convergence part.
           d: distance between two points.

       Returns:
           **x_chamber**; list for x-coordinates of the points on
                      outer contour.
           **y_chamber**; list for y-coordinates of the points on
                      outer contour.
                           '''

    # list for coordinates of points on convergence part.
    l2_x = list()
    l2_y = list()

    # get quantity of the points on convergence part.
    n_l2 = math.floor(l2 / d)

    # set coordinate for every points on convergence part.
    for i in range(n_l2):
        l2_x.append(x_exp - i * d * np.cos(alpha))
        l2_y.append(y_exp + i * d * np.sin(alpha))

    # set the last point on convergence part. Hold the exact location of
    # the last point.
    l2_x.append(x_exp - l2 * np.cos(alpha))
    l2_y.append(y_exp + l2 * np.sin(alpha))

    # coordinate of the first point on the rounding circular arc between
    # convergence and parallel part.
    x0 = l2_x[n_l2]
    y0 = l2_y[n_l2]

    # calculate the angle between horizontal axis and the line, which
    # between circular center and the first point on circular arc.
    b = 0.5 * np.pi - alpha

    # coordinate of center of the circular arc.
    y_r = y0 - r * np.sin(b)
    x_r = x0 - r * np.cos(b)

    # lists for coordinates of points on the circular arc.
    r_x = list()
    r_y = list()

    # the last point of the circular arc/
    # the first point of the parallel part.
    x_r_end = x_r
    y_r_end = y_r + r

    r_l = r * alpha  # length of the circular arc.
    n = math.ceil(r_l / d)  # circular arc can be partitioned in n parts.
    a_step = alpha / n  # the angle for every part of circular arc.

    # set coordinate for every point on the circular arc.
    for i in range(n - 1):
        r_x.append(x_r + r * np.cos(b + (i + 1) * a_step))
        r_y.append(y_r + r * np.sin(b + (i + 1) * a_step))

    # set lists for coordinates of points on parallel part.
    l1_x = list()
    l1_y = list()

    # set the first point on the parallel part.
    l1_x.append(x_r_end)
    l1_y.append(y_r_end)

    # get the x-coordinate of the last point on the parallel part.
    x_end = x_r_end - l1

    # set the coordinate for every point on the parallel part.
    i2 = 0
    while l1_x[i2] > x_end:
        l1_x.append(l1_x[i2] - d)
        l1_y.append(l1_y[i2])
        i2 = i2 + 1

    # assemble all points on outer contour of combustion chamber.
    x_chamber = l2_x + r_x + l1_x
    y_chamber = l2_y + r_y + l1_y

    # reverse the order of all points. Because the points should be
    # outputted with increasing x-coordinate.
    x_chamber.reverse()
    y_chamber.reverse()

    return x_chamber, y_chamber