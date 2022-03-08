import numpy as np


def ideal_aerospike_contour(gamma, N, Ma_e, Short, form):
    '''calculate the contour of a aerospike.

       Args:
           gamma: heat capacity of the gas.
           N: quantity of points on aerospike.
           Ma_e: Mach number at exit.
           Short: Truncate of aerospike.
           form: 0 for linear, other value for annular aerospike.

       Returns:
            **x_Ye or x_out**; list for relationship of x/y_ex
                           x is the x-coordinates of points on aerospike
                           y_ex is the y-coordinate of expansion point.
            **y_Ye or y_out**; list for relationship of x/y_ex
                           y is the y-coordinates of points on aerospike
                           y_ex is the y-coordinate of expansion point.
            **nu_e**; total turning angle.
                    '''
    L = (100 - Short) / 100  # find the completeness of the spike contour.
    Ma = list()  # list for Mach number.
    Ma.append(1)  # first value of Mach number is 1.

    # list for the angle between flow velocity and mach wave [rad].
    mu = list()

    nu = list()  # list for turning angle.
    phi = list()  # list for phi.
    y_Ye = list()  # list for value y/Ye.
    x_Ye = list()  # list for x/Ye.
    d_Ma = (Ma_e - 1) / N  # get delta Mach number.

    # get the total flow turning angle.
    nu_e = ((gamma + 1) / (gamma - 1)) ** 0.5 * np.arctan(((gamma - 1) /
           (gamma + 1) * (Ma_e ** 2 - 1)) ** 0.5) -                           \
           np.arctan((Ma_e ** 2 - 1) ** 0.5)

    # get the expansion ratio.
    exp_r = 1 / (Ma_e / ((1 + ((gamma - 1) / (gamma + 1)) *
                          ((Ma_e ** 2) - 1)) ** ((gamma + 1) /
                                                 (2 * (gamma - 1)))))

    for i in range(int(N)):
        mu.append(np.arcsin(1 / Ma[i]))  # mu value for i-th Mach number.

        # nu value for i-th Mach number.
        nu.append(((gamma + 1) / (gamma - 1)) ** 0.5 *
                  np.arctan(((gamma - 1) /
                             (gamma + 1) * (Ma[i] ** 2 - 1)) ** 0.5) -
                  np.arctan((Ma[i] ** 2 - 1) ** 0.5))

        # find the phi value.
        phi.append(nu_e - nu[i] + mu[i])

        # get the value for Rx/Re or (Rx/Re)^2
        value = ((1 - (2 / (gamma + 1) *
                       (1 + (gamma - 1) / 2 * Ma[i] ** 2)) **
                  ((gamma + 1) / (2 * gamma - 2)) * np.sin(phi[i]) / exp_r))

        # if linear Aerospike.
        if form == 0:
            # set the value in list y/Ye.
            y_Ye.append(value)

        # if annular Aerospike
        else:
            # set the value in list y/Ye.
            y_Ye.append(value ** 0.5)

        # find the x/Ye, X / Ye = (1 - y / Ye) / tan(phi)
        x_Ye.append((1 - y_Ye[i]) / np.tan(phi[i]))

        # set the next Mach number in the list for Mach number.
        Ma.append(Ma[i] + d_Ma)

    # if without truncate.
    if L == 1:
        return x_Ye, y_Ye, nu_e

    # if with truncate.
    else:
        x_out = list()  # set list for x/Ye.
        y_out = list()  # set list for y/Ye.
        x_max = max(x_Ye)  # get the maximal value of x/Ye.

        # get the maximal value of x/Ye after truncate
        x_max_out = x_max * L

        i0 = 0  # set index to 0.

        # set the value in new lists for x/Ye and y/Ye until maximal
        # x/Ye value.
        while x_Ye[i0] < x_max_out:
            x_out.append(x_Ye[i0])
            y_out.append(y_Ye[i0])
            i0 = i0 + 1

        return x_out, y_out, nu_e