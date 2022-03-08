import numpy as np
from scipy.optimize import fsolve


def parabola_nozzle_with_length(omega, s, L, epsilon, n):
    '''calculate the coordinates of points,
       which build a parabolic nozzle.

       Args:
           omega: central angle of rounding arc at throat.
           s: throat radius.
           L: length of the nozzle.
           epsilon: nozzle ratio, (area exit / area throat).
           n: number of points on parabolic curve.

       Returns:
           **x**; x-coordinates.
           **y**; y-coordinates.
           **theta_E**; angle between flow direction and symmetry axis at exit.
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

    # radius of rounding arc after throat.
    R = 0.382 * s

    # calculate the coordinates for every point on the rounding arc.
    for i in range(n_omega):
        x_rounding.append(x_throat + R * np.sin(i * d_omega))
        y_rounding.append(y_throat + R * (1 - np.cos(i * d_omega)))

    # get the coordinate of transit point, which is between
    # rounding arc and parabolic curve.
    x_transition = R * np.sin(omega)
    y_transition = y_throat + R * (1 - np.cos(omega))

    # get the coordinate of point at exit.
    x_exit = L
    y_exit = np.sqrt(s * epsilon)

    # the nozzle curve is a second-degree polynomial(parabola):
    # y(x) = a * x**2 + b * x + c
    # the coordinates of transition point and exit point,
    # and the slop(tan(omega)) at transition point are known,
    # so we can get 3 equations.
    # equation 1:
    # y_transition = a * x_transition ** 2 + b * x_transition + c
    # equation 2: y_exit = a * x_exit ** 2 + b * x_exit + c
    # equation 3: tan(omega) = 2 * a * x_transition + b
    # solve the 3 equations through matrix.
    x_matrix = np.array([[x_transition ** 2, x_transition, 1],
                         [x_exit ** 2, x_exit, 1],
                         [2 * x_transition, 1, 0]])
    y_matrix = np.array([[y_transition],
                         [y_exit],
                         [np.tan(omega)]],
                        dtype='object')
    x_matrix_inverse = np.linalg.inv(x_matrix)
    coefficients = np.dot(x_matrix_inverse, y_matrix)
    a = float(coefficients[0])
    b = float(coefficients[1])
    c = float(coefficients[2])

    # if the angle between flow direction and symmetry axis is smaller
    # than 0 at exit, the entered nozzle length is too big,
    # or the entered central angle of the rounding arc (omega) is too big.
    # calculate the available maximum value of nozzle length with
    # this omega, and return 0,0,maximum length as signal
    # for the above situation.
    if (2 * a * x_exit + b) < 0:
        L_max = -(b / (2 * a))
        return 0, 0, L_max

    else:
        # get the angle between flow direction
        # and symmetry axis at exit.
        theta_E = np.arctan(2 * a * x_exit + b)

        # get the step length for the points on parabolic curve.
        L_parabola = L - x_transition
        d_L_parabola = L_parabola / n

        # set list for points on parabolic curve.
        x_parabola = list()
        y_parabola = list()

        # calculate every point on parabolic curve.
        for i in range(n + 1):
            x_parabola.append(x_transition + i * d_L_parabola)
            y_parabola.append(
                a * (x_parabola[i] ** 2) + b * x_parabola[i] + c)

        # assemble all points.
        x = x_rounding + x_parabola
        y = y_rounding + y_parabola

        return x, y, theta_E


def parabola_nozzle_with_exit_angle(omega, s, theta_e, epsilon, n):
    '''calculate the coordinates of points,
       which build a parabolic nozzle.

       Args:
           omega: central angle of rounding arc at throat.
           s: throat radius.
           theta_e: angle at nozzle exit.
           epsilon: nozzle ratio, (area exit / area throat).
           n: number of points on parabolic curve.

       Returns:
           **x**; x-coordinates.
           **y**; y-coordinates.
           **theta_E**; angle between flow direction and symmetry axis at exit.
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

    # radius of rounding arc after throat.
    R = 0.382 * s

    # calculate the coordinates for every point on the rounding arc.
    for i in range(n_omega):
        x_rounding.append(x_throat + R * np.sin(i * d_omega))
        y_rounding.append(y_throat + R * (1 - np.cos(i * d_omega)))

    # get the coordinate of transit point, which is between
    # rounding arc and parabolic curve.
    x_transition = R * np.sin(omega)
    y_transition = y_throat + R * (1 - np.cos(omega))

    # get the coordinate of point at exit.
    y_exit = np.sqrt(s * epsilon)

    # the nozzle curve is a second-degree polynomial(parabola):
    # y(x) = a * x**2 + b * x + c
    # the coordinates of transition point, the slop(tan(omega))
    # at transition point,
    # the y coordinate at exit and the slop(tan(theta_e)) at exit
    # are known, so we can get 4 equations.
    # equation 1:
    # y_transition = a * x_transition ** 2 + b * x_transition + c
    # equation 2: y_exit = a * x_exit ** 2 + b * x_exit + c
    # equation 3: tan(omega) = 2 * a * x_transition + b
    # equation 4: tan(theta_e) = 2 * a * x_exit + b
    # solve the 4 equations to get the coefficients a, b, c for parabola
    # and the x-coordinate at exit point.
    def func(x):
        a, b, c, x_e = x[0], x[1], x[2], x[3]
        return [a * x_transition ** 2 + b * x_transition + c -
                y_transition,
                2 * a * x_transition + b - np.tan(omega),
                2 * a * x_e + b - np.tan(theta_e),
                a * x_e ** 2 + b * x_e + c - y_exit]

    coefficients = fsolve(func, [1, 1, 1, x_transition])
    a, b, c, x_exit = coefficients[0], coefficients[1], coefficients[2], \
                      coefficients[3]

    L = x_exit

    # if the angle between flow direction and symmetry axis is smaller
    # than 0 at exit, then return 0,0,1 as signal for this error.
    if theta_e < 0:

        return 0, 0, 1

    else:

        # get the step length for the points on parabolic curve.
        L_parabola = L - x_transition
        d_L_parabola = L_parabola / n

        # set list for points on parabolic curve.
        x_parabola = list()
        y_parabola = list()

        # calculate every point on parabolic curve.
        for i in range(n + 1):
            x_parabola.append(x_transition + i * d_L_parabola)
            y_parabola.append(
                a * (x_parabola[i] ** 2) + b * x_parabola[i] + c)

        # assemble all points.
        x = x_rounding + x_parabola
        y = y_rounding + y_parabola

        return x, y, L