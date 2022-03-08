import numpy as np


class IdealNozzle():

    def expansion_angle_from_Ma(self, gamma, Ma):
        '''calculate the prandtl-mayer expansion angle through the Mach number.
           formule (25).

           Args:
               gamma: heat capacity ratio of the gas.
               Ma: Mach number.

           Returns:
               **psi**; prandtl-mayer expansion angle.
               '''

        # calculate the expansion angle with the formule.
        psi = 1/2 * (np.sqrt((gamma + 1)/(gamma - 1)) *
                     np.arctan(np.sqrt((Ma**2 - 1)*(gamma - 1) / (gamma + 1)))
                     -np.arctan(np.sqrt(Ma**2 - 1)))

        return psi


    def Ma_from_expansion_angle(self, gamma, psi, delta):
        '''calculate the Mach number through prandtl-mayer expansion angle.
           formule (25).

           Args:
               gamma: heat capacity ratio of the gas.
               psi: prandtl-mayer expansion.
               delta: toleranz of the calculate result.

           Returns:
               **Ma**; Mach number.
               '''
        # set the Mach number = 2 for first value of psi.
        Ma_1 = 5
        psi_cal = self.expansion_angle_from_Ma(gamma, Ma_1)

        # find, in which interval is the Mach number to be calculated.
        if psi_cal > psi:
            Ma_L = 1
            Ma_R = Ma_1
        else:
            while(psi_cal < psi):
                Ma_1 = Ma_1 + 5
                psi_cal = self.expansion_angle_from_Ma(gamma, Ma_1)
            Ma_L = Ma_1
            Ma_R = Ma_1 + 5

        # set the initial value of Mach number and psi_cal
        Ma = 0.5 * (Ma_L + Ma_R)
        psi_cal = self.expansion_angle_from_Ma(gamma, Ma)
        while(np.absolute(psi_cal - psi) > delta):
            if psi_cal - psi < 0:
                Ma_L = Ma
                Ma = 0.5 * (Ma_L + Ma_R)
                psi_cal = self.expansion_angle_from_Ma(gamma, Ma)
            else:
                Ma_R = Ma
                Ma = 0.5 * (Ma_L + Ma_R)
                psi_cal = self.expansion_angle_from_Ma(gamma, Ma)

        return Ma


    def tau_from_mach(self, gamma, Ma):
        '''calculate tau through Mach number, formule (13a).

           Args:
               gamma: heat capacity ratio.
               Ma: Mach number.

           Returns:
               **tau**; tau'''
        # calculate tau through Ma number, formula (13a).
        tau = np.sqrt((((2 / (gamma + 1)) +
                          ((gamma - 1) * Ma ** 2) / (gamma + 1)) **
                         ((gamma + 1) / (2 * (gamma - 1)))) / Ma)

        return tau



    def ideal_nozzle_curve(self, gamma, Ma_E, s, n):
        '''
           calculate the coordinates of the points on the curve of a perfect
           shock-free supersonic nozzle. The contour of a perfect shock-free
           nozzle can be determined through integration of the characteristic
           equations of axially symmetric flow. Foelsch had developed a
           mathematical methode to approximate the calculation of the Laval
           nozzle contour, which produces parallel and uniform supersonic flow.

           Args:
               gamma(float): heat capacity ratio. :math:`\gamma`
               Ma_E(float): Mach number at exit. :math:`Ma_E`
               s(float): radius of throat. :math:`s`
               n(float): number of the points. :math:`n`

           Returns:
               **x(list)**; x-coordinate of the points on nozzle curve.
               **y(list)**; y-coordinate of the points on nozzle curve.

           | For more details:
           | [1] Foelsch - The Analytical Design of an Axially Symmetric
                           Laval Nozzle for a Parallel and Uniform Jet
           | [2] Shapiro - The Dynamics and Thermodynamics of
                           Compressible Fluid Flow - Volume 2 -
                           Axially Symmetric Supersonic Flow

           4 parameters are required for this calculation,
           the heat capacity ratio of the gas :math:`\gamma`
           the Mach number at exit of the nozzle :math:`Ma_E`
           the throat radius :math:`s`
           the number of points on the nozzle curve  :math:`n`

           tau is determined for every Mach number as:
           formula (13a) in [1]

           :math:`\\tau\ =\\sqrt{\\frac{\\left(\\frac{2}{\\gamma+1}+
           \\frac{\\gamma-1}{\\gamma+1}\\bullet{\\rm Ma}^2\\right)^
           \\frac{\\left(\\gamma+1\\right)}{2\\left(\\gamma-1\\right)}}{Ma}}`



          :math:`F_{\\theta_P}` is determined for every Mach number as:
           formula (16a) in [1]


          :math:`F\\left(\\theta_P\\right)\\ =\\sqrt{\\sin^2\\left(\\theta_P\\right)+2\\bullet\\left[\\cos{\\left(\\theta_P\\right)-\\cos{\\left(\\omega\\right)}}\\right]\\bullet\\left[\\sqrt{{{\\rm Ma}_P}^2-1}\\bullet\\sin{\\left(\\theta_P\\right)+\\cos{\\left(\\theta_P\\right)}}\\right]}`


           y coordinate of every point on the curve is determined as:
           formula (20) in [1]

           :math:`y\\ =\\ \\left[\\frac{D}{4\\bullet\\sin{\\left(\\frac{\\omega}{2}\\right)}}\\right]\\bullet\\left(\\frac{\\tau_P}{\\tau_E}\\right)F\\left(\\theta_P\\right)`

           x coordinate of every point on the curve is determined as:
           formula (21) in [1]

           :math:`x\\ =\\ \\left[\\frac{s\\bullet\\tau_P}{2\\bullet\\sin{\\left(\\frac{\\omega}{2}\\right)}}\\right]\\bullet\\frac{1+\\left[\\cos{\\left(\\theta_P\\right)\\bullet\\sqrt{{\\mathrm{Ma}}^2-1}}-\\sin{\\left(\\theta_P\\right)}\\right]F\\left(\\theta_P\\right)}{\\sin{\\left(\\theta_P\\right)}\\bullet\\sqrt{{\\mathrm{Ma}}^2-1}+\\cos{\\left(\\theta_P\\right)}}`

                  '''

        #get the final expansion angle(psi_E) through final Mach number.
        psi_E = self.expansion_angle_from_Ma(gamma, Ma_E)

        # get the tau at exit(self.tau_E)
        self.tau_E = self.tau_from_mach(gamma, Ma_E)

        # get the angle between direction of flow and symmetric axis(theta_A)
        # at point A. formula (27)
        theta_A = 0.5 * psi_E

        #get the expansion angle at point A(psi_A). formule (27)
        psi_A = theta_A

        #get the omega. formule (27)
        self.omega = psi_A

        #get the Mach number at point A(Ma_A)
        Ma_A = self.Ma_from_expansion_angle(gamma, psi_A, 0.00000001)

        #get the tau at pointA(self.tau_A)
        self.tau_A = self.tau_from_mach(gamma, Ma_A)

        # set initial value.
        Ma_p = Ma_A

        # get the step length of Mach number.
        d_Ma = (Ma_E - Ma_A) / n

        # set lists for x- and y-coordinates.
        x = list()
        y = list()

        # calculate every coordinates on the nozzle.
        for p in range(n + 1):

            psi_p = self.expansion_angle_from_Ma(gamma, Ma_p)

            # formula (26)
            theta_p = psi_E - psi_p

            tau_p = self.tau_from_mach(gamma, Ma_p)

            # formula (16a)
            F_p = np.sqrt(np.sin(theta_p) ** 2 +
                  2 * (np.cos(theta_p) - np.cos(self.omega)) *
                           (np.sqrt(Ma_p ** 2 - 1) *
                            np.sin(theta_p) + np.cos(theta_p)))

            # formula (28) into formula (20)
            # (replace D in formula 34 with formula (28) D = 2s*tau_E )

            y.append((s / (2 * np.sin(self.omega / 2))) * tau_p * F_p)

            # formula (28) into formula (21)
            # (replace D in formula 34 with formula (28) D = 2s*tau_E )
            x.append((s / (2 * np.sin(self.omega / 2))) * tau_p
                      * ((1 + (np.cos(theta_p) * np.sqrt(Ma_p ** 2 - 1) -
                      np.sin(theta_p)) * F_p) / ((np.sin(theta_p) *
                      np.sqrt(Ma_p ** 2 - 1)) + np.cos(theta_p))))

            # new Mach number for next step.
            Ma_p = Ma_p + d_Ma

        return x, y


    def ideal_nozzle_throat(self, R, lambda_throat, s):
        '''calculate the coordinates of points on the throat area
           of a perfect nozzle.

           Args:
               R: radius of the arc at throat.
               lambda_throat: length of lambda, straight line
                              between arc at throat and nozzle curve.
               s: radius of the throat.

           Returns:
               **x_throat**; relationship between x-coordinate of points
                    on throat area and throat radius.
               **y_throat**; relationship between y-coordinate of points
                    on throat area and throat radius.
                    '''

        # set list for x- and y-coordinates.
        x_throat = list()
        y_throat = list()

        # coordinate of the point at throat.
        x_T = 0
        y_T = s

        # set the number of points on rounding arc at throat.
        n_omega = 10
        # calculate the central angle between 2 adjacent points
        # on rounding arc.
        d_omega = self.omega / n_omega

        # get the
        for i in range(n_omega):
            x_throat.append(x_T + R * np.sin(i * d_omega))
            y_throat.append(y_T + R * (1 - np.cos(i * d_omega)))

        # get the coordinate of point D, formula (30)
        x_D = R * np.sin(self.omega)
        y_D = s + R * (1 - np.cos(self.omega))

        n_lambda = 10
        d_lambda = lambda_throat / n_lambda

        for i in range(n_lambda):
            x_throat.append(x_D + i * d_lambda * np.cos(self.omega))
            y_throat.append(y_D + i * d_lambda * np.sin(self.omega))

        return x_throat, y_throat


    def ideal_nozzle_with_R(self, gamma, Ma_E, s, n, R):
        '''calculate the coordinates of the points which build a perfect
           nozzle with existing rounding throat-radius.

           Args:
               gamma: heat capacity of the gas.
               Ma_E: Mach number at exit.
               s: radius of the throat.
               n: number of points on the nozzle curve.
               R: rounding radius at throat.

           Returns:
               **x**; x-coordinates of the points, which build
                  a perfect nozzle.
               **y**; y-coordinates of the points, which build
                  a perfect nozzle.
               or
               **R_max**; available maximum value of the
                      rounding radius at throat.
               **0**; 0
           '''

        # get the coordinates of points on nozzle curve.
        x_curve, y_curve = self.ideal_nozzle_curve(gamma, Ma_E, s, n)

        # get the length of lambda between arc of throat and nozzle curve.
        # formula (28) into formula (32)
        lambda_throat = (s / (2 * np.sin(self.omega / 2))) *                  \
                        (self.tau_A - (1 / np.cos(self.omega / 2))) -         \
                        R * np.tan(self.omega / 2)

        # if lambda_throat < 0, the entered R is larger than the
        # available maximum value, so set the Lambda = 0 into formula (32),
        # to calculate the maximum value of R.
        if lambda_throat < 0:
            R_max = (s / (2 * np.sin(self.omega / 2))) *                      \
                        (self.tau_A - (1 / np.cos(self.omega / 2))) /         \
                     np.tan(self.omega / 2)
            return R_max, 0

        # else, calculate the coordinates of points on the throat area, and
        # assemble the coordinates of all points,
        # which build the perfect nozzle.
        else:
            # calculate x_0 with formula (31)
            x_0 = s / np.tan(self.omega) - R * np.tan(self.omega / 2)

            # get the coordinates in throat area.
            x_throat, y_throat =                                              \
            self.ideal_nozzle_throat(R, lambda_throat, s)

            # translate the coordinates of points on nozzle curve
            # in new coordinate-system.
            for i in range(len(x_curve)):
                x_curve[i] = x_curve[i] - x_0

            # assemble coordinates of all points.
            x = x_throat + x_curve
            y = y_throat + y_curve

            return x, y


    def ideal_nozzle_with_lambda(self, gamma, Ma_E, s, n, lambda_throat):
        '''calculate the coordinates of the points, which build a perfect
           nozzle, with existing straight line between throat rounding arc
           and nozzle curve.

           Args:
               gamma: heat capacity of the gas.
               Ma_E: Mach number at exit.
               s: radius of the throat.
               n: number of points on the nozzle curve.
               lambda_throat: length of the straight line between
                             rounding arc at throat and nozzle curve.

           Returns:
               **x**; x-coordinates of the points, which build
                  a perfect nozzle.
               **y**; y-coordinates of the points, which build
                  a perfect nozzle.
               or
               **0**; 0
               **lambda_max**; available maximum value of the length of the
                           straight line between rounding arc and nozzle curve.
                   '''

        # get the coordinates of points on nozzle curve.
        x_curve, y_curve = self.ideal_nozzle_curve(gamma, Ma_E, s, n)

        # get the radius of arc at throat. formula (28) into formula (32).
        R = ((s / (2 * np.sin(self.omega / 2))) *
            (self.tau_A - (1 / np.cos(self.omega / 2))) - lambda_throat ) /   \
            np.tan(self.omega / 2)

        # if R < 0, the entered lambda_throat is larger than the
        # available maximum value, so set the R = 0 into formula (32),
        # to calculate the maximum value of lambda_throat.
        if R < 0:
            lambda_max = (s / (2 * np.sin(self.omega / 2))) *                 \
                        (self.tau_A - (1 / np.cos(self.omega / 2)))
            return 0, lambda_max

        # else, calculate the coordinates of points on the throat area, and
        # assemble the coordinates of all points,
        # which build the perfect nozzle.
        else:
            # calculate x_0 through formula (31).
            x_0 = s / np.tan(self.omega) - R * np.tan(self.omega / 2)

            # get the coordinates in throat area.
            x_throat, y_throat =                                              \
            self.ideal_nozzle_throat(R, lambda_throat, s)

            # translate the coordinates of points on nozzle curve
            # in new coordinate-system.
            for i in range(len(x_curve)):
                x_curve[i] = x_curve[i] - x_0

            # assemble the coordinates of all points.
            x = x_throat + x_curve
            y = y_throat + y_curve

            return x, y