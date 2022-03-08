import numpy as np
import ideal_nozzle
from scipy.optimize import fsolve


def tictop_nozzle_with_length(gamma, Ma_tic, L_tic, R, L, s, epsilon, n):
    '''calculate a tictop-nozzle with existing rounding arc at throat.

       Args:
           gamma: heat capacity of the gas.
           Ma_tic: designed Mach number of a perfect nozzle at exit.
                   (curve of TIC-part)
           L_tic: length of TIC-part.
           R: radius of rounding arc at throat.
           L: total length of the nozzle.
           s: throat radius.
           epsilon: nozzle ratio.(exit area / throat area)
           n: number of points on the nozzle curve.

       Returns:
           **x;** x-coordinates.
           **y;** y-coordinates.
           **theta_E;** angle between flow direction and symmetry axis at exit.
           or
           **1,0,0 and similar**; signal for corresponding faulty situation.
           '''

    # use the function to generate a perfect nozzle-curve.
    tic = ideal_nozzle.IdealNozzle()
    x_perfect, y_perfect = \
        tic.ideal_nozzle_with_R(gamma, Ma_tic, s, n, R)

    x_tic = list()
    y_tic = list()
    i = 0
    while (x_perfect[i] < L_tic):
        x_tic.append(x_perfect[i])
        y_tic.append(y_perfect[i])
        i = i + 1
    n_tic = len(x_tic) - 1

    # get the coordinate of point at exit.
    x_exit = L
    y_exit = np.sqrt(s * epsilon)

    # if the cross-sectional area at transition point is already
    # larger than exit area, the entered designed Mach number of TIC is
    # too large or the length of TIC-part is too long.
    # return (0,1,0) as signal for the above situation.
    if y_tic[n_tic] > y_exit:
        return 0, 1, 0

    else:
        # get the slope at transition point.
        slope_transition = (y_tic[n_tic] - y_tic[n_tic - 1]) / \
                           (x_tic[n_tic] - x_tic[n_tic - 1])

        # the parabolic part is a second-degree polynomial:
        # y(x) = a * x**2 + b * x + c
        # the coordinates of transition point and exit point,
        # and the slope at transition point are known,
        # so we can get 3 equations.
        # equation 1:
        # y_transition = a * x_transition ** 2 + b * x_transition + c
        # equation 2: y_exit = a * x_exit ** 2 + b * x_exit + c
        # equation 3: slope_transition_point = 2 * a * x_transition + b
        # solve the 3 equations through matrix.
        x_matrix = np.array([[x_tic[n_tic] ** 2, x_tic[n_tic], 1],
                             [x_exit ** 2, x_exit, 1],
                             [2 * x_tic[n_tic], 1, 0]])
        y_matrix = np.array([[y_tic[n_tic]],
                             [y_exit],
                             [slope_transition]],
                            dtype='object')
        x_matrix_inverse = np.linalg.inv(x_matrix)
        coefficients = np.dot(x_matrix_inverse, y_matrix)
        a = float(coefficients[0])
        b = float(coefficients[1])
        c = float(coefficients[2])

        # if a > 0, the designed Mach number of TIC-nozzle is too small.
        # return (1,0,0) as signal for the above situation.
        if a > 0:
            return 1, 0, 0

        else:
            # if the angle between flow direction and symmetry axis is
            # smaller than 0 at exit, the entered nozzle length is too big,
            # or the length of TIC-part is too long,
            # calculate the available maximum value of nozzle length with
            # this TIC-length, and return (0,0,maximum length) as signal
            # for the above situation.
            if (2 * a * x_exit + b) < 0:
                L_max = -(b / (2 * a))
                return 0, 0, L_max
            else:
                # get the angle between flow direction
                # and symmetry axis at exit.
                theta_E = np.arctan(2 * a * x_exit + b)

                # get the step length for the points on parabolic curve.
                L_parabola = L - x_tic[n_tic]
                n_parabola = n - n_tic
                d_L_parabola = L_parabola / n_parabola

                # set list for points on parabolic curve.
                x_parabola = list()
                y_parabola = list()

                # calculate every point on parabolic curve.
                for i in range(n_parabola):
                    x_parabola.append(
                        x_tic[n_tic] + (i + 1) * d_L_parabola)
                    y_parabola.append(
                        a * (x_parabola[i] ** 2) + b * x_parabola[i] + c)

                # assemble all points.
                x = x_tic + x_parabola
                y = y_tic + y_parabola

                return x, y, theta_E


def tictop_nozzle_with_exit_angle(gamma, Ma_tic, L_tic, R, theta_e,
                                  s, epsilon, n):
    '''calculate a tictop-nozzle with existing rounding arc at throat.

       Args:
           gamma: heat capacity of the gas.
           Ma_tic: designed Mach number of a perfect nozzle at exit.
                   (curve of TIC-part)
           L_tic: length of TIC-part.
           R: radius of rounding arc at throat.
           L: total length of the nozzle.
           s: throat radius.
           epsilon: nozzle ratio.(exit area / throat area)
           n: number of points on the nozzle curve.

       Returns:
           **x**; x-coordinates.
           **y**; y-coordinates.
           **theta_E**; angle between flow direction and symmetry axis at exit.
           or
           **1,0,0 and similar**; signal for corresponding faulty situation.
           '''

    # use the function to generate a perfect nozzle-curve.
    tic = ideal_nozzle.IdealNozzle()
    x_perfect, y_perfect = \
        tic.ideal_nozzle_with_R(gamma, Ma_tic, s, n, R)

    x_tic = list()
    y_tic = list()
    i = 0
    while (x_perfect[i] < L_tic):
        x_tic.append(x_perfect[i])
        y_tic.append(y_perfect[i])
        i = i + 1
    n_tic = len(x_tic) - 1

    # get the coordinate of point at exit.
    y_exit = np.sqrt(s * epsilon)

    # if the cross-sectional area at transition point is already
    # larger than exit area, the entered designed Mach number of TIC is
    # too large or the length of TIC-part is too long.
    # return (0,1,0) as signal for the above situation.
    if y_tic[n_tic] > y_exit:
        return 0, 1, 0

    else:
        # get the slope at transition point.
        slope_transition = (y_tic[n_tic] - y_tic[n_tic - 1]) / \
                           (x_tic[n_tic] - x_tic[n_tic - 1])

        # the parabolic part is a second-degree polynomial:
        # y(x) = a * x**2 + b * x + c
        # the coordinates of transition point and exit point,
        # and the slope at transition point are known,
        # so we can get 3 equations.
        # equation 1:
        # y_transition = a * x_transition ** 2 + b * x_transition + c
        # equation 2: y_exit = a * x_exit ** 2 + b * x_exit + c
        # equation 3: slope_transition_point = 2 * a * x_transition + b
        # solve the 3 equations through matrix.

        def func(x):
            a, b, c, x_e = x[0], x[1], x[2], x[3]
            return [a * x_tic[n_tic] ** 2 + b * x_tic[n_tic] + c -
                    y_tic[n_tic],
                    2 * a * x_tic[n_tic] + b - slope_transition,
                    2 * a * x_e + b - np.tan(theta_e),
                    a * x_e ** 2 + b * x_e + c - y_exit]

        coefficients = fsolve(func, [1, 1, 1, x_tic[n_tic]])
        a, b, c, x_exit = coefficients[0], coefficients[1], \
                          coefficients[2], coefficients[3]

        L = x_exit

        # if a > 0, the designed Mach number of TIC-nozzle is too small.
        # return (1,0,0) as signal for the above situation.
        if a > 0:
            return 1, 0, 0

        else:
            # if the angle between flow direction and symmetry axis is
            # smaller than 0 at exit, the entered nozzle length is too big,
            # or the length of TIC-part is too long,
            # calculate the available maximum value of nozzle length with
            # this TIC-length, and return (0,0,maximum length) as signal
            # for the above situation.
            if (2 * a * x_exit + b) < 0:
                L_max = -(b / (2 * a))
                return 0, 0, L_max
            else:
                # get the angle between flow direction
                # and symmetry axis at exit.
                theta_E = np.arctan(2 * a * x_exit + b)

                # get the step length for the points on parabolic curve.
                L_parabola = L - x_tic[n_tic]
                n_parabola = n - n_tic
                d_L_parabola = L_parabola / n_parabola

                # set list for points on parabolic curve.
                x_parabola = list()
                y_parabola = list()

                # calculate every point on parabolic curve.
                for i in range(n_parabola):
                    x_parabola.append(
                        x_tic[n_tic] + (i + 1) * d_L_parabola)
                    y_parabola.append(
                        a * (x_parabola[i] ** 2) + b * x_parabola[i] + c)

                # assemble all points.
                x = x_tic + x_parabola
                y = y_tic + y_parabola

                return x, y, L