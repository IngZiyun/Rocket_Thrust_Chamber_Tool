import numpy as np


def cylinder_combustion_chamber(s, r1, r2, R_cylinder, L_cylinder, alpha):
    '''calculate the coordinates of points, which build a cylinder
       combustion chamber.

       Args:
           s: throat radius.
           r1: radius of rounding arc between throat and convergent part.
           r2: radius of rounding arc between cylinder part and
               convergent part.
           R_cylinder: radius of the cylinder combustion chamber.
           L_cylinder: length of the cylinder combustion chamber.
           alpha: convergent angle of the convergent part.

       Returns:
           **x**; x-coordinates.
           **y**; y-coordinates.
           '''

    # calculate the length of straight line of convergent part.
    L_convergence = ((R_cylinder - s) -
                     r1 * (1 - np.cos(alpha)) -
                     r2 * (1 - np.cos(alpha))) / np.sin(alpha)

    # calculate the total length of combustion chamber in the direction
    # of symmetry axis. Length of cylinder part plus
    # length of convergent part in direction of symmetry axis.
    L_total = L_cylinder + (r1 + r2) * np.sin(alpha) + \
              L_convergence * np.cos(alpha)

    # set the number of points, which build the cylinder part,
    # and get the step length between 2 points.
    n_cylinder = 50
    d_L_cylinder = L_cylinder / n_cylinder

    # list for coordinates of those points.
    x_cylinder = list()
    y_cylinder = list()

    # get the coordinates of every point of cylinder part.
    for i in range(n_cylinder + 1):
        x_cylinder.append(-L_total + i * d_L_cylinder)
        y_cylinder.append(R_cylinder)

    # set the number of points, which build the rounding arc between
    # cylinder part and convergent part,
    # and get the step length between 2 points.
    n_r2 = 10
    d_alpha_2 = alpha / n_r2

    # list for coordinates of those points.
    x_r2 = list()
    y_r2 = list()

    # get the coordinates of every point on straight part
    # of convergent part.
    for i in range(n_r2):
        x_r2.append(x_cylinder[n_cylinder] +
                    r2 * np.sin((i + 1) * d_alpha_2))
        y_r2.append(y_cylinder[n_cylinder] -
                    r2 * (1 - np.cos((i + 1) * d_alpha_2)))

    # get number of points, which build the straight part of
    # convergent part, and get the step length between 2 points.
    n_convergence = 10
    d_L_convergence = L_convergence / n_convergence

    # list for coordinates of those points.
    x_convergence = list()
    y_convergence = list()

    # get the coordinates of every point on the straight part of
    # convergent part.
    for i in range(n_convergence):
        x_convergence.append(x_r2[n_r2 - 1] +
                             (i + 1) * d_L_convergence * np.cos(alpha))
        y_convergence.append(y_r2[n_r2 - 1] -
                             (i + 1) * d_L_convergence * np.sin(alpha))

    # set the number of points, which build the rounding arc between
    # convergent part and throat, and get the step length between 2 points.
    n_r1 = 10
    d_alpha_1 = alpha / n_r1

    # list for coordinates of those points.
    x_r1 = list()
    y_r1 = list()

    # get the coordinates of every point on the rounding arc between
    # convergent part and throat.
    for i in range(n_r1 - 1):
        x_r1.append(x_convergence[n_convergence - 1] +
                    r1 * (np.sin(alpha) -
                          np.sin(alpha - (i + 1) * d_alpha_1)))
        y_r1.append(y_convergence[n_convergence - 1] -
                    r1 * (1 - np.cos(alpha)) +
                    r1 * (1 - np.cos(alpha - (i + 1) * d_alpha_1)))

    # assemble all points.
    x = x_cylinder + x_r2 + x_convergence + x_r1
    y = y_cylinder + y_r2 + y_convergence + y_r1

    # get the volume of the generated combustion chamber.
    chamber_volume = L_cylinder * (np.pi * R_cylinder ** 2) + \
                     (((R_cylinder - s) / np.tan(alpha)) * np.pi / 3) * \
                     (R_cylinder ** 2 + s ** 2 + R_cylinder * s)

    return x, y, chamber_volume