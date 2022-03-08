import numpy as np
import gasdynamic
from scipy.optimize import fsolve

def aerospike_approximate_wang(gamma, theta_p, Ma_e, n, form):

    # form = 0 means linear aerospike
    R_E = 1
    x_E = 0

    nu_e = ((gamma + 1) / (gamma - 1)) ** 0.5 * np.arctan(((gamma - 1) /
           (gamma + 1) * (Ma_e ** 2 - 1)) ** 0.5) -                           \
           np.arctan((Ma_e ** 2 - 1) ** 0.5)

    alpha = np.arcsin(1 / Ma_e)

    nu_p = nu_e - alpha - theta_p

    def fun_get_mach_from_nu(x, angle):

        nu = angle[0]
        return ((gamma + 1) / (gamma - 1)) ** 0.5 * np.arctan(((gamma - 1) /
           (gamma + 1) * (x ** 2 - 1)) ** 0.5) -                           \
           np.arctan((x ** 2 - 1) ** 0.5) - nu

    Ma_p = fsolve(fun_get_mach_from_nu, Ma_e, args=[nu_p])[0]

    x_F = R_E / np.tan(alpha)

    R_F = 0

    gas_dynamic = gasdynamic.Gasdynamic(gamma)

    area_ratio = gas_dynamic.area_ratio_from_ma(Ma_e)
    area_ratio_primary = gas_dynamic.area_ratio_from_ma(Ma_p)
    area_ratio_end_to_primary = area_ratio / area_ratio_primary

    # In the approximate methode of wang, the expansion angle at the
    # transition point D is a half of total expansion angle at aerospike.
    nu_D = nu_e - 0.5 * alpha

    Ma_D = fsolve(fun_get_mach_from_nu, Ma_p, args=[nu_D])[0]

    # angle of expansion wave is at transition area a half of total
    # expansion wave, and the angle between expansion wave and x-axis is
    # minus [total expansion angle(nu_e) -
    # expansion angle of current area(0.5*nu_e)].
    phi_D = nu_D - nu_p + np.arcsin(1 / Ma_D)


    if not (theta_p + alpha) > (nu_D - nu_p):
        return 0, 0, 0, 0


    area_ratio_D = gas_dynamic.area_ratio_from_ma(Ma_D)
    
    area_ratio_D_to_primary = area_ratio_D / area_ratio_primary


    def fun_parabola(x, coordinates):
        a, b, c = x[0], x[1], x[2]
        x_C, x_D, R_C, R_D = coordinates[0],      \
                                 coordinates[1], coordinates[2],              \
                                 coordinates[3]

        return [a * x_C ** 2 + b * x_C + c - R_C,
                2 * a * x_D + b - np.tan(-(nu_D-nu_p)),
                a * x_D ** 2 + b * x_D + c -
                R_D]

    def fun_cubic(x, coordinates):
        a, b, c, d = x[0], x[1], x[2], x[3]
        x_D, x_F, R_D, R_F = coordinates[0],      \
                                 coordinates[1], coordinates[2],              \
                                 coordinates[3]

        return [a * x_D ** 3 + b * x_D ** 2 + c *
                x_D + d - R_D,
                3 * a * x_D ** 2 + 2 * b * x_D + c -
                np.tan(-(nu_D-nu_p)),
                a * x_F ** 3 + b * x_F ** 2 + c * x_F + d - R_F,
                3 * a * x_F ** 2 + 2 * b * x_F + c - 0]

    if form == 0:

        area_EC = R_E * area_ratio_end_to_primary
        def fun_C_point_for_linear(x):
            x, R = x[0], x[1]
            return [np.sqrt(x ** 2 + (R - 1) ** 2) - area_EC,
                    np.tan(0.5 * np.pi - alpha) * x + 1 - R]
        C_point = fsolve(fun_C_point_for_linear, [0, 0])
        x_C = C_point[0]
        R_C = C_point[1]

        R_D = R_E - (area_EC / area_ratio_D_to_primary)
        x_D = (R_E - R_D) / np.tan(phi_D)

    else:

        area_EC = (np.pi * R_E ** 2) * area_ratio_end_to_primary

        def fun_C_point_for_annular(x):

            x, R = x[0], x[1]
            return [np.tan(0.5 * np.pi - alpha) * x + 1 - R,
                    np.pi * (R + R_E) * np.sqrt((x - x_E) ** 2 +
                                    (R - R_E) ** 2) - area_EC]

        C_point = fsolve(fun_C_point_for_annular, [0, 0])
        x_C = C_point[0]
        R_C = C_point[1]

        R_D =                                                        \
        np.sqrt(R_E ** 2 - (area_EC / area_ratio_D_to_primary) / np.pi)
        x_D = (R_E - R_D) / np.tan(phi_D)

    coeff_pa = \
        fsolve(fun_parabola, [1, 1, 1],
               args=[x_C, x_D, R_C, R_D])

    a_pa, b_pa, c_pa = coeff_pa[0], coeff_pa[1], \
                       coeff_pa[2]

    coeff_cubic = fsolve(fun_cubic, [1, 1, 1, 1],
                         args=[x_D, x_F, R_D, R_F])

    a_cubic, b_cubic, c_cubic, d_cubic = \
        coeff_cubic[0], coeff_cubic[1], \
        coeff_cubic[2], coeff_cubic[3]

    # if second order derivatives of the cubic polynomial is less than 0,
    # the theta_p is too small, return 0,0,0,0 as signal for this case.
    if 6 * a_cubic * x_D + 2 * b_cubic < 0:
        return 0, 0, 0, 0

    d_x = (x_F - x_C) / (n - 1)
    x_coordinates = list()
    R_coordinates = list()

    for i in range(n):
        x = i * d_x + x_C
        x_coordinates.append(x)
        if x < x_D:
            R = a_pa * x ** 2 + b_pa * x + c_pa

        else:
            R = a_cubic * x ** 3 + b_cubic * x ** 2 + c_cubic * x + d_cubic

        R_coordinates.append(R)

    return x_coordinates, R_coordinates, Ma_p, area_ratio_primary

