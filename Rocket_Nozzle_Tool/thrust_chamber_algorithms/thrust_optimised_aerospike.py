import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def thrust_optimised_aerospike_contour(epsilon, L, Ma_E, theta_E, gamma,
                                      n, diagram):


    # set the coordinate for the lip point.
    R_E = 1
    x_E = 0

    theta_star = theta_E - np.sqrt((gamma + 1) / (gamma - 1)) * np.arctan(
        np.sqrt(((gamma - 1) * (Ma_E ** 2 - 1)) / (gamma + 1))) + np.arctan(
        np.sqrt(Ma_E ** 2 - 1))

    phi_star = theta_star - (np.pi / 2)

    alpha_E = np.arcsin(1 / Ma_E)
    phi_E = theta_E - alpha_E

    Ma_E_star = critical_mach_from_mach(Ma_E, gamma)

    constant_A = Ma_E_star * (
                np.cos(theta_E) + np.tan(alpha_E) * np.sin(-theta_E))
    constant_B = R_E * ((Ma_E_star * np.sin(theta_E)) ** 2) * np.tan(
        alpha_E) * ((1 + ((gamma - 1) / 2) * Ma_E ** 2) ** (- 1 / (gamma - 1)))


    def func_MaStar_R(x, R_control):

        R = R_control[0]

        return R * ((1 + ((gamma - 1) / 2) * (((2 * x ** 2) / (gamma + 1)) /
               (1 - (((gamma - 1) * x ** 2) / (gamma + 1))))) **
               (-1 / (gamma - 1))) * ((x * ((((-2 * constant_A) / x) *
               (np.sqrt((1 - ((gamma - 1) * (x ** 2) / (gamma + 1))) /
               (x ** 2 - 1))) + (np.sqrt((4 * constant_A ** 2 / (x ** 2)) *
               ((1 - (((gamma - 1) * (x ** 2)) / (gamma + 1))) / (x ** 2 - 1))
               - (4 * ((2 * (x ** 2)) / (gamma + 1)) / (x ** 2 - 1)) *
               (((constant_A ** 2) / x ** 2) - 1)))) / (2 * ((2 * x ** 2 /
               (gamma + 1)) / (x ** 2 - 1))))) ** 2) * (np.sqrt((1 -
               (((gamma - 1) * (x ** 2)) / (gamma + 1))) / (x ** 2 - 1)))     \
               - constant_B



    def sin_theta_at_control_surface(x):
        '''x: critical Mach number'''

        sin_theta = (((-2 * constant_A) / x) * (np.sqrt((1 - ((gamma - 1) *
                    (x ** 2) / (gamma + 1))) / (x ** 2 - 1))) +
                    (np.sqrt((4 * constant_A ** 2 / (x ** 2)) * ((1 -
                    (((gamma - 1) * (x ** 2)) / (gamma + 1))) / (x ** 2 - 1)) -
                    (4 * ((2 * (x ** 2)) / (gamma + 1)) / (x ** 2 - 1)) *
                    (((constant_A ** 2) / x ** 2) - 1)))) / (2 * ((2 * x ** 2 /
                    (gamma + 1)) / (x ** 2 - 1)))

        # original 4 * constant_A, without **2 !!!
        return sin_theta

    control_x = list()
    control_R = list()

    # set a array for all data of the points on the control surface.
    data_control_surface = np.zeros(shape=(n, 10))

    # data format: x(0), R(1), Ma_star(2), theta(3), alpha(4), lambda_L(5),
    # lambda_R(6), eta(7), beta_L(8), beta_R(9)

    d_x = L / (n - 1)

    data_control_surface[0, 0] = 0

    control_x.append(data_control_surface[0, 0])

    data_control_surface[0, 1] = 1

    control_R.append(data_control_surface[0, 1])

    data_control_surface[0, 2] = critical_mach_from_mach(Ma_E, gamma)
    data_control_surface[0, 3] = theta_E
    data_control_surface[0, 4] = np.abs(
        np.arctan(tan_alpha_from_critical_mach(data_control_surface[0, 2],
                                               gamma)))
    data_control_surface[0, 5] = np.tan(
        data_control_surface[0, 3] + data_control_surface[0, 4])
    data_control_surface[0, 6] = np.tan(
        data_control_surface[0, 3] - data_control_surface[0, 4])
    data_control_surface[0, 7] = 1 / data_control_surface[0, 2] * np.tan(
        data_control_surface[0, 4])
    data_control_surface[0, 8] = \
        (np.sin(data_control_surface[0, 3]) * np.sin(
            data_control_surface[0, 4])) / (
                    data_control_surface[0, 1] * np.cos(
                data_control_surface[0, 3] + data_control_surface[0, 4]))
    data_control_surface[0, 9] = \
        (np.sin(data_control_surface[0, 3]) * np.sin(
            data_control_surface[0, 4])) / (
                    data_control_surface[0, 1] * np.sin(
                data_control_surface[0, 3] - data_control_surface[0, 4]))
    for i in range(1, n):
        data_control_surface[i, 0] = d_x * i
        control_x.append(d_x * i)
        data_control_surface[i, 1] = data_control_surface[i - 1, 1] + d_x * (
            np.tan(data_control_surface[i - 1, 3] - data_control_surface[
                i - 1, 4]))
        control_R.append(data_control_surface[i, 1])
        data_control_surface[i, 2] = fsolve(func_MaStar_R,
                                            data_control_surface[i - 1, 2],
                                            args=[data_control_surface[
                                                      i, 1]])[0]
        data_control_surface[i, 3] = np.arcsin(
            sin_theta_at_control_surface(data_control_surface[i, 2]))
        data_control_surface[i, 4] = np.abs(np.arctan(
            tan_alpha_from_critical_mach(data_control_surface[i, 2], gamma)))
        data_control_surface[i, 5] = np.tan(
            data_control_surface[i, 3] + data_control_surface[i, 4])
        data_control_surface[i, 6] = np.tan(
            data_control_surface[i, 3] - data_control_surface[i, 4])
        data_control_surface[i, 7] = 1 / data_control_surface[i, 2] * np.tan(
            data_control_surface[i, 4])
        data_control_surface[i, 8] = \
            (np.sin(data_control_surface[i, 3]) * np.sin(
                data_control_surface[i, 4])) / (
                        data_control_surface[i, 1] * np.cos(
                    data_control_surface[i, 3] + data_control_surface[i, 4]))
        data_control_surface[i, 9] = \
            (np.sin(data_control_surface[i, 3]) * np.sin(
                data_control_surface[i, 4])) / (
                        data_control_surface[i, 1] * np.sin(
                    data_control_surface[i, 3] - data_control_surface[i, 4]))

    # get the throat-point.
    R_T = np.sqrt(R_E ** 2 - ((np.cos(theta_star) * R_E ** 2) / epsilon))
    A_E = np.pi * R_E ** 2
    A_star = A_E / epsilon

    def func_x_T(x):
        return np.pi * (R_E + R_T) * np.sqrt(
            (x_E - x) ** 2 + (R_E - R_T) ** 2) - A_star

    x_T = fsolve(func_x_T, -1.0)[0]

    # delta Mach number around the lip point E.
    d_Ma = (Ma_E - 1) / (n - 1)

    # data form: x(0), R(1), Ma_star(2), theta(3), alpha(4), lambda_R(5),
    # eta(6), beta_R(7)
    data_at_lip = np.zeros(shape=(n, 8))
    for i in range(n):
        Ma = Ma_E - i * d_Ma
        theta = theta_E - np.sqrt((gamma + 1) / (gamma - 1)) * np.arctan(
            np.sqrt(((gamma - 1) * ((Ma_E - (n - 1 - i) * d_Ma) ** 2 - 1)) / (
                        gamma + 1))) + np.arctan(
            np.sqrt((Ma_E - (n - 1 - i) * d_Ma) ** 2 - 1))
        alpha = np.arcsin(1 / Ma)
        Ma_star = critical_mach_from_mach(Ma, gamma)
        data_at_lip[i, 0] = 0
        data_at_lip[i, 1] = 1
        data_at_lip[i, 2] = Ma_star
        data_at_lip[i, 3] = theta
        data_at_lip[i, 4] = alpha
        data_at_lip[i, 5] = np.tan(theta - alpha)
        data_at_lip[i, 6] = 1 / (np.tan(alpha) * Ma_star)
        data_at_lip[i, 7] = (np.sin(theta) * np.sin(alpha)) / (
                    1 * np.sin(theta - alpha))

    def characteristic(x_R, R_R, Ma_star_R, theta_R, alpha_R, lambda_R, eta_R,
                       beta_R, x_L, R_L, Ma_star_L, theta_L, alpha_L, lambda_L,
                       eta_L, beta_L):

        def iteration(x_N, R_N, Ma_star_N, theta_N, lambda_N, eta_N, beta_N):

            x_ite = ((lambda_R * x_R - lambda_N * x_N) + (R_N - R_R)) / (
                    lambda_R - lambda_N)
            R_ite = R_N - lambda_N * (x_N - x_ite)
            Ma_star_ite = (theta_R - theta_N + eta_N * Ma_star_N + eta_R *
                          Ma_star_R - beta_R * (R_R - R_ite) - beta_N *
                          (x_N - x_ite)) / (eta_N + eta_R)
            theta_ite = theta_N - eta_N * (Ma_star_ite - Ma_star_N) + beta_N *\
            (x_N - x_ite)
            alpha_ite =                                                       \
            np.abs(np.arctan(tan_alpha_from_critical_mach(Ma_star_ite, gamma)))
            lambda_L_ite = np.tan(theta_ite + alpha_ite)
            lambda_R_ite = np.tan(theta_ite - alpha_ite)
            eta_ite = 1 / (Ma_star_ite * np.tan(alpha_ite))
            beta_L_ite = np.sin(theta_ite) * np.sin(alpha_ite) / (
                        R * np.cos(theta_ite + alpha_ite))
            beta_R_ite = np.sin(theta_ite) * np.sin(alpha_ite) / (
                        R * np.sin(theta_ite - alpha_ite))

            return x_ite, R_ite, Ma_star_ite, theta_ite, alpha_ite,           \
                   lambda_L_ite, lambda_R_ite, eta_ite, beta_L_ite, beta_R_ite

        x = ((lambda_R * x_R - lambda_L * x_L) + (R_L - R_R)) / (
                lambda_R - lambda_L)
        R = R_L - lambda_L * (x_L - x)
        Ma_star = (theta_R - theta_L + eta_L * Ma_star_L + eta_R *
                   Ma_star_R - beta_R * (R_R - R) - beta_L * (x_L - x)) / (
                   eta_L + eta_R)
        if Ma_star < 1:
            theta = 0
            alpha = 0
            lambda_L_N = 0
            lambda_R_N = 0
            eta = 0
            beta_L_N = 0
            beta_R_N = 0
        else:
            theta = theta_L - eta_L * (Ma_star_L - Ma_star) + beta_L * (
                        x_L - x)
            alpha = np.abs(np.arctan(tan_alpha_from_critical_mach(Ma_star,
                                                                  gamma)))
            lambda_L_N = np.tan(theta + alpha)
            lambda_R_N = np.tan(theta - alpha)
            eta = 1 / (Ma_star * np.tan(alpha))
            beta_L_N = np.sin(theta) * np.sin(alpha) / (
                        R * np.cos(theta + alpha))
            beta_R_N = np.sin(theta) * np.sin(alpha) / (
                        R * np.sin(theta - alpha))
            tolerance_theta = 1
            times = 0

            while tolerance_theta > 0.000001:
                theta_0 = theta
                x, R, Ma_star, theta, alpha, lambda_L_N, lambda_R_N, eta,     \
                beta_L_N, beta_R_N =                                          \
                iteration(x, R, Ma_star, theta, lambda_L_N, eta, beta_L_N)
                tolerance_theta = np.abs(theta_0 - theta)
                times = times + 1
                if times > 5:
                    tolerance_theta = 0.0000001

        return x, R, Ma_star, theta, alpha, lambda_L_N, lambda_R_N, eta,      \
               beta_L_N, beta_R_N

    # set a list for the data of all characteristic points.
    data_cha = list()

    # data format: x(0), R(1), Ma_star(2), theta(3), alpha(4), lambda_L(5),
    # lambda_R(6), eta(7), beta_L(8), beta_R(9)

    # set the first element of data_cha as a array for the points on the
    # first left-running characteristic line.
    data_cha.append(np.zeros(shape=(n, 10)))

    for i in range(10):
        data_cha[0][0, i] = data_control_surface[1, i]

    # set list for coordinates of points on the
    # first left-running characteristic line.

    chara_R1_x = list()
    chara_R1_R = list()
    chara_R1_x.append(data_cha[0][0, 0])
    chara_R1_R.append(data_cha[0][0, 1])

    for i in range(n - 1):

        if data_cha[0][i, 2] < 1:
            break

        data_cha[0][i + 1, 0], data_cha[0][i + 1, 1], data_cha[0][i + 1, 2], \
        data_cha[0][i + 1, 3], data_cha[0][i + 1, 4], data_cha[0][i + 1, 5], \
        data_cha[0][i + 1, 6], data_cha[0][i + 1, 7], data_cha[0][i + 1, 8], \
        data_cha[0][i + 1, 9] = \
            characteristic(data_at_lip[i + 1, 0], data_at_lip[i + 1, 1],
                           data_at_lip[i + 1, 2], data_at_lip[i + 1, 3],
                           data_at_lip[i + 1, 4], data_at_lip[i + 1, 5],
                           data_at_lip[i + 1, 6], data_at_lip[i + 1, 7],
                           data_cha[0][i, 0], data_cha[0][i, 1],
                           data_cha[0][i, 2],
                           data_cha[0][i, 3], data_cha[0][i, 4],
                           data_cha[0][i, 5],
                           data_cha[0][i, 7], data_cha[0][i, 8])

        n_R_total = i

        # set coordinates for the points on the first left-running
        # characteristic line.
        chara_R1_x.append(data_cha[0][i + 1, 0])
        chara_R1_R.append(data_cha[0][i + 1, 1])

    n_R = n_R_total

    chara_R_x = list()
    chara_R_R = list()

    # data format: x(0), R(1), Ma_star(2), theta(3), alpha(4), lambda_L(5),
    # lambda_R(6), eta(7), beta_L(8), beta_R(9)
    for i in range(n - 2):

        data_cha.append(np.zeros(shape=(n_R, 10)))
        chara_R_x.clear()
        chara_R_R.clear()

        chara_R_x.append(data_control_surface[i + 2, 0])
        chara_R_R.append(data_control_surface[i + 2, 1])

        for j in range(n_R - 1):

            # set data of first point for every left-running
            # characteristic line.
            for i2 in range(10):
                data_cha[i + 1][0, i2] = data_control_surface[i + 2, i2]
            # characteristic(Right-Parameters(0,1,2,3,4,6,7,9),
            # Left-Parameters(0,1,2,3,4,5,7,8))
            data_cha[i + 1][j + 1, 0], data_cha[i + 1][j + 1, 1], \
            data_cha[i + 1][j + 1, 2], \
            data_cha[i + 1][j + 1, 3], data_cha[i + 1][j + 1, 4], \
            data_cha[i + 1][j + 1, 5], \
            data_cha[i + 1][j + 1, 6], data_cha[i + 1][j + 1, 7], \
            data_cha[i + 1][j + 1, 8], \
            data_cha[i + 1][j + 1, 9] = \
                characteristic(data_cha[i][j + 1, 0], data_cha[i][j + 1, 1],
                               data_cha[i][j + 1, 2],
                               data_cha[i][j + 1, 3], data_cha[i][j + 1, 4],
                               data_cha[i][j + 1, 6],
                               data_cha[i][j + 1, 7], data_cha[i][j + 1, 9],
                               data_cha[i + 1][j, 0], data_cha[i + 1][j, 1],
                               data_cha[i + 1][j, 2],
                               data_cha[i + 1][j, 3], data_cha[i + 1][j, 4],
                               data_cha[i + 1][j, 5],
                               data_cha[i + 1][j, 7], data_cha[i + 1][j, 8])

            chara_R_x.append(data_cha[i + 1][j + 1, 0])
            chara_R_R.append(data_cha[i + 1][j + 1, 1])

        # if the characteristic net is wanted, then plot the every
        # left-running characteristic line.
        if diagram == 1:

            plt.plot(chara_R_x, chara_R_R, 'k', linewidth=0.3)
            plt.plot([0, data_cha[i + 1][n_R - 1, 0]],
                     [1, data_cha[i + 1][n_R - 1, 1]], 'k', linewidth=0.3)

        # the points on the left-running characteristic line are less than
        # the points on the last left-running characteristic line.
        if i > (n - 1 - n_R_total):
            n_R = n_R - 1

    # list for coordinates of the points on aerospike contour.
    contour_x = list()
    contour_R = list()

    def get_point_on_contour(x_L, R_L, x_1, R_1, x_2, R_2, theta_1, theta_2):
        '''get the coordinate of the point on aerospike contour (streamline)
        through the characteristic line.'''

        def func_get_contour_x(x):
            return ((R_L - R_1 + (x_1 - x) * ((R_1 - R_2) / (x_1 - x_2))) / (
                    x_L - x)) - (np.tan((((x_1 - x) * theta_2) / (x_1 - x_2))
                    + (((x - x_2) * theta_1) / (x_1 - x_2))))

        # get the coordinate x and R(y) of the point on aerospike contour.
        P_x = fsolve(func_get_contour_x, x_2)[0]
        P_R = R_1 - (x_1 - P_x) * ((R_1 - R_2) / (x_1 - x_2))

        return P_x, P_R

    # set the coordinate of the first point(point D, exit point)
    # on the contour in the coordinate-list.
    contour_x.append(data_control_surface[n - 1, 0])
    contour_R.append(data_control_surface[n - 1, 1])

    # set the coordinate of the first point(point D, exit point)
    # on the contour as the first value of calculate loop.
    point_x = data_control_surface[n - 1, 0]
    point_R = data_control_surface[n - 1, 1]

    # set beginning number, which determines the characteristic points for
    # every contour point, for the loop.
    number = 0
    number2 = 1
    number1 = 0

    while point_x > x_T:
        xL, RL = point_x, point_R
        x2, R2, theta2 = data_cha[n - 2 - number][number2, 0], \
                         data_cha[n - 2 - number][number2, 1], \
                         data_cha[n - 2 - number][number2, 3]
        x1, R1, theta1 = data_cha[n - 3 - number][number1, 0], \
                         data_cha[n - 3 - number][number1, 1], \
                         data_cha[n - 3 - number][number1, 3]
        point_x, point_R = get_point_on_contour(xL, RL, x1, R1, x2, R2, theta1,
                                                theta2)
        if point_R < 1:
            contour_x.append(point_x)
            contour_R.append(point_R)

        if point_R < data_cha[n - 3 - number][number1 + 1, 1]:
            number1 = number1 + 1
            number2 = number2 + 1
            if number2 > (n_R_total - 2):
                point_x, point_R = x_T, R_T
                contour_x.append(point_x)
                contour_R.append(point_R)
        else:
            number = number + 1
            if (n - 3 - number) < 0:
                point_x, point_R = x_T, R_T
                contour_x.append(point_x)
                contour_R.append(point_R)


    if diagram == 1:

        plt.plot(control_x, control_R, 'b', linewidth=2,
                 label='control surface')
        plt.plot(contour_x, contour_R, 'r', linewidth=3,
                 label='contour of aerospike')
        plt.plot(0, 1, 'k', label='characteristic net')
        plt.title('Thrust Optimised Aerospike Nozzle', fontsize='large')
        plt.legend()
        plt.show()

    return contour_x, contour_R


def mach_from_critical_mach(ma_star, gamma):
    ma = np.sqrt(
        (ma_star ** 2) / (1 - (((gamma - 1) * (ma_star ** 2 - 1)) / 2)))
    return ma


def critical_mach_from_mach(ma, gamma):
    ma_star = np.sqrt(
        (ma ** 2) / (1 + (((gamma - 1) / (gamma + 1)) * (ma ** 2 - 1))))
    return ma_star


def tan_alpha_from_critical_mach(critical_Ma, gamma):
    tan_alpha = np.sqrt(
        (1 - (((gamma - 1) * critical_Ma ** 2) / (gamma + 1))) / (
                    critical_Ma ** 2 - 1))

    return tan_alpha

def rho_from_ma(Ma, gamma):

    roh = ((1+((Ma**2)*(gamma-1)/2))**(1/(gamma-1)))
        # calculate density for every Mach number and set the value
        # in the list for density.

    return roh



def get_parameter_through_mach_exit_angle(Ma_E, theta_E, gamma):


    n = 1000
    R_E = 1
    alpha_E = np.arcsin(1 / Ma_E)
    phi_E = theta_E - alpha_E

    Ma_E_star = critical_mach_from_mach(Ma_E, gamma)

    constant_A = Ma_E_star * (
            np.cos(theta_E) + np.tan(alpha_E) * np.sin(-theta_E))
    constant_B = R_E * ((Ma_E_star * np.sin(theta_E)) ** 2) * np.tan(
        alpha_E) * ((1 + ((gamma - 1) / 2) * Ma_E ** 2) ** (- 1 / (gamma - 1)))

    def func_MaStar_R(x, R_control):

        R = R_control[0]

        return R * ((1 + ((gamma - 1) / 2) * (((2 * x ** 2) / (gamma + 1)) /
               (1 - (((gamma - 1) * x ** 2) / (gamma + 1))))) **
               (-1 / (gamma - 1))) * ((x * ((((-2 * constant_A) / x) *
               (np.sqrt((1 - ((gamma - 1) * (x ** 2) / (gamma + 1))) /
               (x ** 2 - 1))) + (np.sqrt((4 * constant_A ** 2 / (x ** 2)) *
               ((1 - (((gamma - 1) * (x ** 2)) / (gamma + 1))) / (x ** 2 - 1))
               - (4 * ((2 * (x ** 2)) / (gamma + 1)) / (x ** 2 - 1)) *
               (((constant_A ** 2) / x ** 2) - 1)))) / (2 * ((2 * x ** 2 /
               (gamma + 1)) / (x ** 2 - 1))))) ** 2) * (np.sqrt((1 -
               (((gamma - 1) * (x ** 2)) / (gamma + 1))) / (x ** 2 - 1)))     \
               - constant_B



    def sin_theta_at_control_surface(x):
        '''x: critical Mach number'''

        sin_theta = (((-2 * constant_A) / x) * (np.sqrt((1 - ((gamma - 1) *
                    (x ** 2) / (gamma + 1))) / (x ** 2 - 1))) +
                    (np.sqrt((4 * constant_A ** 2 / (x ** 2)) * ((1 -
                    (((gamma - 1) * (x ** 2)) / (gamma + 1))) / (x ** 2 - 1)) -
                    (4 * ((2 * (x ** 2)) / (gamma + 1)) / (x ** 2 - 1)) *
                    (((constant_A ** 2) / x ** 2) - 1)))) / (2 * ((2 * x ** 2 /
                    (gamma + 1)) / (x ** 2 - 1)))

        # original 4 * constant_A, without **2 !!!
        return sin_theta

    data_control_surface = np.zeros(shape=(n, 10))

    # data format: x(0), R(1), Ma_star(2), theta(3), alpha(4), lambda_L(5),
    # lambda_R(6), eta(7), beta_L(8), beta_R(9)

    d_R = R_E / n

    data_control_surface[0, 0] = 0
    data_control_surface[0, 1] = 1
    data_control_surface[0, 2] = critical_mach_from_mach(Ma_E, gamma)
    data_control_surface[0, 3] = theta_E
    data_control_surface[0, 4] = np.abs(
        np.arctan(tan_alpha_from_critical_mach(data_control_surface[0, 2],
                                               gamma)))

    for i in range(1, n):
        data_control_surface[i, 1] = R_E - i * d_R
        data_control_surface[i, 0] = (data_control_surface[i, 1] -
                                      data_control_surface[i - 1, 1]) / (
            np.tan(data_control_surface[i - 1, 3] - data_control_surface[
                i - 1, 4])) + data_control_surface[i - 1, 0]
        data_control_surface[i, 2] = fsolve(func_MaStar_R,
                                            data_control_surface[i - 1, 2],
                                            args=[data_control_surface[
                                                      i, 1]])[0]
        data_control_surface[i, 3] = np.arcsin(
            sin_theta_at_control_surface(data_control_surface[i, 2]))
        data_control_surface[i, 4] = np.abs(np.arctan(
            tan_alpha_from_critical_mach(data_control_surface[i, 2], gamma)))

    different = list()
    for i in range(n):
        ma = mach_from_critical_mach(data_control_surface[i, 2], gamma)

        value = ((2 * np.sqrt(ma ** 2 - 1))/(gamma * ma ** 2)) -              \
        np.sin(-2 * data_control_surface[i, 3])

        different.append(np.abs(value))

    index_D = different.index(min(different))

    x_D = data_control_surface[index_D, 0]
    R_D = data_control_surface[index_D, 1]

    length = x_D

    roh_throat = 1 / rho_from_ma(1, gamma)

    roh_roh_star = list()
    for i in range(n):
        mach = mach_from_critical_mach(data_control_surface[i, 2], gamma)
        roh =1 / rho_from_ma(mach, gamma)
        roh_roh_star.append(roh / roh_throat)

    value_0 = 0
    for i in range(index_D + 1):

        d_value = roh_roh_star[i] * data_control_surface[i, 2] *              \
        (np.sin(data_control_surface[i, 4]) /
         np.sin(data_control_surface[i, 4] - data_control_surface[i, 3])) *   \
        2 * data_control_surface[i, 1] * d_R

        value_0 = value_0 + d_value

    epsilon = 1 / value_0


    return length, epsilon


def thrust_optimised_aerospike_nozzle(Ma_E, theta_E, gamma, n, diagram):

    length, epsilon =                                                         \
    get_parameter_through_mach_exit_angle(Ma_E, theta_E, gamma)
    x_coordinate, y_coordinate =                                              \
    thrust_optimised_aerospike_contour(epsilon, length, Ma_E, theta_E, gamma,
                                       n, diagram)

    return x_coordinate, y_coordinate