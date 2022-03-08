import numpy as np
from scipy.optimize import fsolve
import math


class Aerospike_Geometry():

    def __init__(self, form):

        '''to confirm annular or linear.
               Args:
                   form: 0 for linear, others for annular.
                   '''
        self.form = form # get the parameters for aerospike-form.



    def supfic_area_cone(self, x1, x2, y1, y2):

        '''Function for calculation of the superficial area of the
           truncated cone.
                Args:
                    x1, y1: x1, y1: x- and y-coordinate of point 1.
                    x2, y2: x2, y2: x- and y-coordinate of point 2.
                Returns:
                    A: the superficial.
                    '''
        A = np.pi * (y1 + y2) * np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
        # equation for calculate the superficical area of a truncated cone.

        return A



    def circle_detemination(self, x1, y1, x2, y2, x3, y3):

        '''calculate the center and radius of a circle through 3 points
                Args:
                    x1, y1: x- and y-coordinate of point 1.
                    x2, y2: x- and y-coordinate of point 2.
                    x3, y3: x- and y-coordinate of point 3.
                Returns:
                    x, y: x- and y-coordinate of the center of the circle
                    r: radius of the circle.
                '''
        x = ((y2 - y1) * (y3 * y3 - y1 * y1 + x3 * x3 - x1 * x1) - (y3 - y1) *
             (y2 * y2 - y1 * y1 + x2 * x2 - x1 * x1)) / (2 * ((x3 - x1) *
             (y2 - y1) - (x2 - x1) * (y3 - y1))) # equation for x-coordinate.

        y = ((x2 - x1) * (x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1) - (x3 - x1) *
             (x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1)) / (2 * ((y3 - y1) *
            (x2 - x1) - (y2 - y1) * (x3 - x1))) # equation for y-coordinate.

        r = np.sqrt((x - x1) ** 2 + (y - y1) ** 2) # equation for radius.

        return x, y, r



    def find_minimal_dis(self, x, y, x_tofind, y_tofind):

        '''find the point form points-group 'B' for every point
           of points-group 'A', which build the minimal distance/
           cross-sectional area with it.
                Args:
                    x:        list of the x-coordinates of the points in A.
                    y:        list of the y-coordinates of the points in A.
                    x_tofind: list of the x-coordinates of the points in B.
                    y_tofind: list of the y-coordinates of the points in B.
                Returns:
                    x_found: list for the x-coordinates
                             of the found points in B.
                    y_found: list for the y-coordinates
                             of the found points in B.
                    '''
        n = len(x) # get the quantity of elements in group A.
        m = len(x_tofind) # get the quantity of elements in group B.
        x_found = list() # set the list for x-coordinates of found points.
        y_found = list() # set the list for y-coordinates of found points.
        dis = list() # set a list for the saving of distance/area.
        if self.form == 0:
            # if linear aerospike, just to find the minimal distance.
            for i in range(n): # for i-th point in group A.
                dis.clear() # clear the list for distance.
                for im in range(m): # for im-th point in group B.
                                    # calculate the distance between i-th point
                                    # in A and every point in B.
                    dis.append(np.sqrt((x[i] - x_tofind[im]) ** 2 +
                                       (y[i] - y_tofind[im]) ** 2))
                    # calculate the distance between i-th point in A and im-th
                    # point in B, then save the distance in distance list.
                dismin = min(dis)# find the minimal value in the distance list.
                n_dismin = dis.index(dismin)
                # find the index of the minimal value in distance list, the
                # index of minimal distance is also the index of the point in
                # B, which build the minimal distance with the i-th point in A.
                x_found.append(x_tofind[n_dismin])# set the x-coordinate.
                y_found.append(y_tofind[n_dismin])# set the y-coordinate.

        else: # if annular aerospike, the minimal area must be found.
            for i in range(n): # for i-th point in group A.
                dis.clear() # clear the list for area.
                for im in range(m): # for im-th point in group B.
                                    # calculate the builded area between i-th
                                    # point in A and every point in B.
                    dis.append(self.supfic_area_cone(x[i], x_tofind[im], y[i],
                                                     y_tofind[im]))
                    # calculate the area between i-th point in A and im-th
                    # point in B, then save the area in area list.
                dismin = min(dis) # find the minimal value in the area list.
                n_dismin = dis.index(dismin)
                # find the index of the minimal value in area list, the
                # index of minimal area is also the index of the point in
                # B, which build the minimal area with the i-th point in A.
                x_found.append(x_tofind[n_dismin])# set the x-coordinate.
                y_found.append(y_tofind[n_dismin])# set the y-coordinate.

        return x_found, y_found



    def separate_inner(self, xin, yin, x_ex, y_ex):
        '''separete the points in combustion chamber, throat-point and the
           points on aerospike from all points on the inner contour.
                Args:
                     xin, yin:   lists for x- and y-coordinates of the points
                                 on inner contour.
                     x_ex, y_ex: coordinate of the expansion-point.
                Returns:
                    x_throat[0], y_throat[0]: coordinate of throat-point.
                    x_cham, y_cham: lists for x- and y-coordiantes of the
                                    points in combustion chamber.
                    x_spike, y_spike: list for x- and y-coordinates of the
                                      points on aerospike.
                         '''
        n = len(xin)  # get the quantity of points on inner contour.
        x_cham = list() # set a list for x-coordinates of points in chamber.
        y_cham = list() # set a list for y-coordinates of points in chamber.
        x_spike = list() # set a list for x-coordinates of points on spike.
        y_spike = list() # set a list for y-coordinates of points on spike.
        xex = list()
        yex = list()
        xex.append(x_ex)
        yex.append(y_ex)
        #translate the x- and y-coordinate of the expansion-point to two lists.
        x_throat, y_throat = self.find_minimal_dis(xex, yex, xin, yin)
        # find the throat point from the all points on inner contour. Because
        # the throat-point and expansion-point build the smallest
        # cross-sectional area.
        m = xin.index(x_throat[0]) # get the index on the throat point.
        for i1 in range(m - 1):
            # set the x- and y-coordinates for points in chamber. The points
            # before throat are in chamber.
            x_cham.append(xin[i1])
            y_cham.append(yin[i1])
        for i2 in range(m, n, 1):
            # set the x- and y-coordinates for points on aerospike. The points
            # after throat are on aerospike.
            x_spike.append(xin[i2])
            y_spike.append(yin[i2])
        return x_throat[0], y_throat[0], x_cham, y_cham, x_spike, y_spike



    def spike_area_ratio(self, x_exp, y_exp, xp, yp):
        '''calculate the area ratio on aerospike.
                Args:
                    x_exp, y_exp: coordinate of the expansion-point.
                    xp, yp: lists for x- and y-coordinates of the points on
                            aerospike
                Returns:
                    area_ratio: list for area ratio on aerospike.
                    '''
        area_ratio = list() # set a list for saving values of area ratio.
        n = len(xp) # get the quantity of the points on aerospike.
        x_ex = list()
        x_ex.append(x_exp)
        y_ex = list()
        y_ex.append(y_exp)
        #translate the x- and y-coordinate of the expansion-point to two lists.
        x_throat, y_throat = self.find_minimal_dis(x_ex, y_ex, xp, yp)
        # get the throat-point.
        x_throat = x_throat[0] # change the list to float for calculation.
        y_throat = y_throat[0] # change the list to float for calculation.

        if self.form == 0: # if linear aerospike.
            A_throat = (np.sqrt((x_exp - x_throat) ** 2 +
                                 (y_exp - y_throat) ** 2))
            # get the throat cross-sectional area.
            for iar in range(n):
                # calculate all area ratios, and set the values in the list
                # for area ratios.
                value = np.sqrt((x_exp - xp[iar]) ** 2 +
                                    (y_exp - yp[iar]) ** 2)
                area_ratio.append(A_throat / value)

            return area_ratio

        else:# if annular aerospike.
            A_throat = self.supfic_area_cone(x_exp, x_throat, y_exp, y_throat)
            # get the throat cross-sectional area.
            for iara in range(n):
                # calculate all area ratios, and set the values in the list
                # for area ratios.
                #value = self.supfic_area_cone(x_exp, xp[iara], y_exp, yp[iara])
                value = self.supfic_area_cone(x_exp, xp[iara], y_exp, yp[iara])
                area_ratio.append(A_throat / value)

            return area_ratio



    def full_area_ratio(self, xin, yin, xout, yout):
        '''calculate the all area ratios along the contour.
                Args:
                    xin, yin: lists for x- and y-coordinates of the points
                              on inner contour(inner points).
                    xout, yout: lists for x- and y-coordinates of the points
                                on the outer contour(outer points).
                Returns:
                    A_ratio_lin or A_ratio_anu: list for all area ratios.
                    '''
        x_ex_in, y_ex_in = self.find_minimal_dis(xout, yout, xin, yin)
        # for every outer point find the inner point, which can build the
        # minimal cross-sectional area with the outer point.
        n_out = len(xout) # get the quantity of the outer points
        l = list() # set a list for saving the values of
                   # distance or cross-sectional areas.
        if self.form == 0: # if linear aerospike, just find the distance.
            for i in range(n_out):
                # calculate the distance of every point-couple,
                # set it in the distance list.
                l.append((xout[i] - x_ex_in[i]) ** 2 +
                         (yout[i] - y_ex_in[i]) ** 2)
        else: # if annular aerospike, find the cross-sectional area.
            for i in range(n_out):
                # calculate the cross-sectional area of every point couple,
                # set it in the list for cross-sectional area.
                l.append(self.supfic_area_cone(xout[i], x_ex_in[i],
                                               yout[i], y_ex_in[i]))
        n_exp = l.index(min(l)) # get the index of minimal value from
                                # the list l.
        self.x_ex = xout[n_exp] # get the coordinate of expansion-point.
                                # the index of the point, which can build the
                                # smallst cross-sectional area with a inner
        self.y_ex = yout[n_exp] # point, is the expansion-point.
        x_out_cham = list() # set lists for x- and y-coordinates of the points
                            # on outer contour of combustion chamber.
        y_out_cham = list()
        for icham in range(n_exp):
            # set the values in lists for x- and y-coordinates of the points
            # on the outer contour of combustion chamber.
            x_out_cham.append(xout[icham])
            y_out_cham.append(yout[icham])
        self.x_throat, self.y_throat, x_in_cham, y_in_cham, self.x_spike,     \
        self.y_spike=self.separate_inner(xin,yin,self.x_ex,self.y_ex)
        # separete the points in combustion chamber, throat-point and the
        # points on aerospike from all points on the inner contour.
        n_p = len(x_in_cham) # get the quantity of points on inner contour of
                             # combustion chamber.
        self.x_r_throat, self.y_r_throat, self.r_throat =                     \
        self.circle_detemination(self.x_throat, self.y_throat,
                                 x_in_cham[n_p - 1],
                                 y_in_cham[n_p - 1],
                                 x_in_cham[n_p - 3],
                                 y_in_cham[n_p - 3])
        # get the rounding radius at throat through the throat point and
        # other 2 points on the rounding.
        if n_exp<len(x_in_cham):# if the points on outer contour of combustion
                                # chamber are less than the points on inner
                                # contour of combustion chamber.
            self.x_in_cham_end, self.y_in_cham_end =                          \
            self.find_minimal_dis(x_out_cham,y_out_cham,x_in_cham,y_in_cham)
            # find the point from all points on inner contour of combustion
            # chamber for every point on the outer contour of combustion
            # chamber, which build the minimal cross-sectional area with it,
            # such points are the final points on inner contour of combustion
            # chamber.
            self.x_out_cham_end = x_out_cham
            self.y_out_cham_end = y_out_cham
            # set lists for coordinates of final points on outer contour of
            # combustion chamber
        else: # if the points on inner contour of combustion chamber are less
              # than the points on outer contour of combustion chamber.
            self.x_out_cham_end, self.y_out_cham_end =                        \
            self.find_minimal_dis(x_in_cham,y_in_cham,x_out_cham,y_out_cham)
            # find the point from all points on outer contour of combustion
            # chamber for every point on the inner contour of combustion
            # chamber, which build the minimal cross-sectional area with it,
            # such points are the final points on outer contour of combustion
            # chamber.
            self.x_in_cham_end = x_in_cham
            self.y_in_cham_end = y_in_cham
              # set lists for coordinates of final points on inner contour of
              # combustion chamber
        n = len(self.x_in_cham_end) # get the quantity of points on
                                    # inner contour of combustion chamber.
        if self.form == 0: # if linear Aerospike
            A_throat_lin =                                                    \
            np.sqrt((self.x_ex-self.x_throat)**2+(self.y_ex-self.y_throat)**2)
            # get the throat area
            self.A_ratio_lin_cham = list()#set list for area ratios in chamber.
            for i in range(n):
                # calculate the area ratios in chamber.
                self.A_ratio_lin_cham.append(A_throat_lin/(np.sqrt(
                    (self.x_in_cham_end[i]-self.x_out_cham_end[i])**2+
                    (self.y_in_cham_end[i]-self.y_out_cham_end[i])**2)))

            self.A_ratio_lin_spike = self.spike_area_ratio(
                           self.x_ex, self.y_ex, self.x_spike, self.y_spike)
            # calculate the area ratios on aerospike.

            A_ratio_lin = self.A_ratio_lin_cham + self.A_ratio_lin_spike
            # assemble the area ratios together.

            return A_ratio_lin
        else: # if annular Aerospike
            A_throat_anu = self.supfic_area_cone(
                self.x_throat, self.x_ex, self.y_throat, self.y_ex)
            # calculate the throat area.
            self.A_ratio_anu_cham = list()#set list for area ratios in chamber.
            for ia in range(n):
                # calculate the area ratios in chamber.
                value = self.supfic_area_cone(self.x_in_cham_end[ia],
                self.x_out_cham_end[ia], self.y_in_cham_end[ia],
                                              self.y_out_cham_end[ia])
                self.A_ratio_anu_cham.append(A_throat_anu/value)

            self.A_ratio_anu_spike = self.spike_area_ratio(self.x_ex,
                                        self.y_ex, self.x_spike, self.y_spike)
            # calculate the area ratios on aerospike.

            A_ratio_anu = self.A_ratio_anu_cham + self.A_ratio_anu_spike
            # assemble the area ratios together.

            return A_ratio_anu



    def points_arrange(self, xin, yin, xout, yout):
        '''arrange the contour-points.
                Args:
                    xin, yin: lists for x- and y-coordinates of the points on
                              inner contour.
                    xout, yout: lists for x- and y-coordinates of the points on
                                outer contour.
                Returns:
                    self.x_in_cham_end: list for x-coordinates of the points on
                                        inner contour of combustion chamber.
                    self.y_in_cham_end: list for y-coordinates of the points on
                                        inner contour of combustion chamber.
                    self.x_throat: x-coordinate of throat point.
                    self.y_throat: y-coordinate of throat point.
                    self.x_spike: list for x-coordinates of the points on
                                  aerospik.
                    self.y_spike: list for y-coordinates of the points on
                                  aerospik.
                    self.x_ex: x-coordinate of expansion point.
                    self.y_ex: y-coordinate of expansion point.
                    self.x_out_cham_end: list for x-coordinates of points on
                                         outer contour of combustion chamber.
                    self.y_out_cham_end: list for y-coordinates of points on
                                         outer contour of combustion chamber.
                                         '''

        self.full_area_ratio(xin, yin, xout, yout)
        # run the function 'full_area-ratio' then the wanted parameters can
        # be got.

        return self.x_in_cham_end, self.y_in_cham_end, self.x_throat,         \
        self.y_throat, self.x_spike, self.y_spike, self.x_ex, self.y_ex,      \
        self.x_out_cham_end, self.y_out_cham_end, self.x_r_throat,            \
        self.y_r_throat, self.r_throat



    def chamber_area_ratio(self, xin, yin, xout, yout, x_throat,
                           y_throat, x_exp, y_exp):
        '''calculate the area ratio in combustion chamber.
                Args:
                    xin, yin: lists for x- and y-coordinates of points on
                              inner contour of combustion chamber.
                    xout, yout: lists for x- and y-coordinates of points on
                                outer contour of combustion chamber.
                    x_throat, y_throat: coordinate of throat point.
                    x_exp, y_exp: coordinate of expansion point.
                Returns:
                    area_ratio: list for area ratio on aerospike.
                    '''
        n = len(xin) # get the quantity of points on inner contour of chamber.
        a_ratio =list()# set a list for area ratio.
        if self.form == 0:# if linear aerospike.
            a_throat = np.sqrt((x_exp-x_throat)**2+(y_exp-y_throat)**2)
            # calculate the throat area.
            for i in range(n):
                # calculate all area ratios.
                a_ratio.append(a_throat/np.sqrt((xin[i]-xout[i])**2+
                                                (yin[i]-yout[i])**2))
        else:# if annular aerospike.
            a_throat = self.supfic_area_cone(x_exp, x_throat, y_exp, y_throat)
            #calculate the throat area.
            for i in range(n):
                # calculate all area ratios.
                a = self.supfic_area_cone(xin[i], xout[i], yin[i], yout[i])
                a_ratio.append(a_throat/a)

        return a_ratio



    def divided_area_ratio(self, xin, yin, xout, yout):
        '''get the area ratios of different parts of a aerospike
           thrust-chamber.
                Args:
                    xin, yin: lists for x- and y-coordinates of points on
                              inner contour.
                    xout, yout: lists for x- and y-coordinates of points on
                                outer contour.
                Returns:
                    self.A_ratio_lin_cham or_anu_cham: list for area ratios in
                                                       chamber.
                    self.A_ratio_lin_spike or_anu_spike: list for area ratios
                                                         on aerospike.
                                                       '''
        self.full_area_ratio(xin, yin, xout, yout)
        # run the function 'full_area_ratio', then the wanted parameters can
        # be got.

        if self.form == 0:# if linear aerospike.

            return self.A_ratio_lin_cham,                                     \
            self.A_ratio_lin_spike, self.x_ex, self.y_ex

        else: # if annular aerospike.

            return self.A_ratio_anu_cham,                                     \
            self.A_ratio_anu_spike, self.x_ex, self.y_ex



    def anu_hydraulic_diameter(self, x_in, y_in, x_out, y_out):
        '''calculate the local hydraulic diameter of
           annular aerospike thrust-chamber.
                Args:
                    x_in, y_in: coordinate of the point on inner contour.
                    x_out, y_out: coordinate of the point on outer contour.
                Returns:
                    dt: local hydraulic diameter.
                    '''
        A = self.supfic_area_cone(x_in, x_out, y_in, y_out)
        #calculate the area.
        U = 2 * np.pi * y_in + 2 * np.pi * y_out# calculate the scope.
        dt = 4 * (A / U) # calculate the hydraulic diameter.

        return dt



    def lin_hydraulic_diameter(self, x_in, y_in, x_out, y_out, depth):
        '''calculate the local hydraulic diameter of
           linear aerospike thrust-chamber.
                Args:
                    x_in, y_in: coordinate of the point on inner contour.
                    x_out, y_out: coordinate of the point on outer contour.
                    depth: depth of the linear aerospike.
                Returns:
                    dt: local hydraulic diameter.
                            '''
        U = 2 * depth + 2 * np.sqrt((x_in - x_out) ** 2 + (y_in - y_out) ** 2)
        # calculate the scope.
        A = depth * np.sqrt((x_in - x_out) ** 2 + (y_in - y_out) ** 2)
        # calculate the area.
        dt = 4 * (A / U) # calculate the hydraulic diameter.

        return dt



    def spike_contour(self, k, N, Ma_e, Short):
        '''calculate the contour of a aerospike.
                Args:
                    k: isentropic exponente.
                    N: quantity of points on aerospike.
                    Ma_e: Mach number at exit.
                    Short: Truncate of aerospike.
                Reterns:
                    x_Ye or x_out: list for relationship of x/y_ex
                                x is the x-coordinates of points on aerospike
                                y_ex is the y-coordinate of expansion point.
                    y_Ye or y_out: list for relationship of x/y_ex
                                y is the y-coordinates of points on aerospike
                                y_ex is the y-coordinate of expansion point.
                    nu_e: total turning angle.
                                '''
        L = (100 - Short) / 100 # find the completeness of the spike contour.
        Ma = list() # list for Mach number.
        Ma.append(1) # firts value of Mach number is 1.
        mu = list()
        # list for the angle between flow velocity and mach wave [rad].
        nu = list() # list for turning angle.
        phi = list() # list for phi.
        y_Ye = list() # list for value y/Ye.
        x_Ye = list() # list for x/Ye.
        d_Ma = (Ma_e - 1) / N # get delta Mach number.
        nu_e = ((k + 1) / (k - 1)) ** 0.5 * np.arctan(((k - 1) / (k + 1) *
               (Ma_e ** 2 - 1)) ** 0.5) - np.arctan((Ma_e ** 2 - 1) ** 0.5)
        # find the total flow turning angle.
        exp_r = 1 / (Ma_e /((1+ ((k-1)/(k+1)) * ((Ma_e ** 2) - 1)) ** ((k+1) /
                (2 * (k-1)))))
        # find the expansion ratio.
        for i in range(int(N)):
            mu.append(np.arcsin(1 / Ma[i]))  # mu value for i-th Mach number.
            nu.append(((k + 1) / (k - 1)) ** 0.5 * np.arctan(((k - 1) /
                      (k + 1) * (Ma[i] ** 2 - 1)) ** 0.5) -
                      np.arctan((Ma[i] ** 2 - 1) ** 0.5))
            # nu value for i-th Mach number.
            phi.append(nu_e - nu[i] + mu[i]) # find the phi value.
            value = ((1 - (2 / (k + 1) * (1 + (k - 1) / 2 * Ma[i] ** 2)) **
                      ((k + 1) / (2 * k - 2)) * np.sin(phi[i]) / exp_r))
            # find the value for Rx/Re or (Rx/Re)^2
            if self.form == 0: # if linear Aerospike.
                y_Ye.append(value) # set the value in list y/Ye.

            else:  # if annular Aerospike
                y_Ye.append(value ** 0.5) # set the value in list y/Ye.

            x_Ye.append((1 - y_Ye[i]) / np.tan(phi[i]))
            # find the x/Ye, X / Ye = (1 - y / Ye) / tan(phi)
            Ma.append(Ma[i] + d_Ma)
            # set the next Mach number in the list for Mach number.
        if L == 1:# if without truncate.
            return x_Ye, y_Ye, nu_e

        else: # if with truncate.
            x_out = list() # set list for x/Ye.
            y_out = list() # set list for y/Ye.
            x_max = max(x_Ye) # get the maximal value of x/Ye.
            x_max_out = x_max * L
            # get the maximal value of x/Ye after truncate
            i0 = 0 # set index to 0.
            while x_Ye[i0] < x_max_out:
                # set the value in new lists for x/Ye and y/Ye until maximal
                # x/Ye value.
                x_out.append(x_Ye[i0])
                y_out.append(y_Ye[i0])
                i0 = i0 + 1
            return x_out, y_out, nu_e



    def inner_chamber_contour(self, x0, y0, angle_degree, r, x_end, d):
        '''generate the inner contour of combustion chamber.
                Args:
                    x0, y0: coordinate of the first point on inner contour.
                    angle_degree: the total turning angle of aerospike.
                    r: rounding radius at throat.
                    x_end: x-coordinate of the last point on inner contour.
                    d: distance between two points.
                Returns:
                    x_chamber: list for x-coordinates of the points on
                               inner contour.
                    y_chamber: list for y-coordinates of the points on
                               inner contour.
                               '''
        a = np.pi * (angle_degree / 180) # translate the angle to radian.
        b = 0.5 * np.pi - abs(a)
        # calculate the angle between horizontal axis and the line, which
        # between circular center and the first point on circular arc.
        y_r = y0 - r * np.sin(b) # y coordinate of center of circle.
        x_r = x0 - r * np.cos(b) # x coordinate of center of circle.
        r_x = list() # list for x-coordinate of point on rounding.
        r_y = list() # list for y-coordinate of point on rounding.
        x_r_end = x_r # x-coordinate for the last point of the rounding/
                      # the first point of the parallel part.
        y_r_end = y_r + r # x-coordinate for the last point of the rounding.
        r_l = r * a # length of the rounding circular arc.
        n = math.ceil(r_l / d) # the circular arc be partitioned in n parts.
        a_step = a / n # the angle for every part of circular arc.

        for i in range(n-1):
            # calculate the coordinate for every point on the circular arc.
            r_x.append(x_r + r * np.cos(b + (i + 1) * a_step))
            r_y.append(y_r + r * np.sin(b + (i + 1) * a_step))

        l_x = list() # list for x-coordinates of the points on parallel part.
        l_y = list() # list for y-coordinates of the points on parallel part.
        l_x.append(x_r_end) # set the first point on the parallel part.
        l_y.append(y_r_end) # set the first point on the parallel part.

        i2 = 0
        while l_x[i2] > x_end:
            # set the coordinate for every point on the parallel part.
            l_x.append(l_x[i2]-d)
            l_y.append(l_y[i2])
            i2 = i2 + 1

        x_chamber = r_x + l_x
        y_chamber = r_y + l_y
        # assemble all points on inner contour of combustion chamber.
        x_chamber.reverse()
        y_chamber.reverse()
        # reverse the order of all points. Because the points should be
        # outputted with increasing x-coordinate.
        return x_chamber, y_chamber



    def outer_chamber_contour(self, x_exp, y_exp, angle_degree, r, l1, l2, d):
        '''generate the outer contour of combustion chamber.
                Args:
                    x_exp, y_exp: coordinate of expansion point.
                    angle_degree: total turning angle of aerospike./
                                  convergence angle.
                    r: rounding radius between convergence part and
                       parallel part.
                    l1: length of parallel part.
                    l2: length of convergence part.
                    d: distance between two points.
                Returns:
                    x_chamber: list for x-coordinates of the points on
                               outer contour.
                    y_chamber: list for y-coordinates of the points on
                               outer contour.
                               '''
        a = np.pi*(angle_degree/180) # translate the angle to radian.
        l2_x = list() # list for x-coordinates of points on convergence part.
        l2_y = list() # list for y-coordinates of points on convergence part.
        n_l2 = math.floor(l2/d)#get quantity of the points on convergence part.
        for i in range(n_l2):
            #set coordinate for every points on convergence part.
            l2_x.append(x_exp - i * d * np.cos(a))
            l2_y.append(y_exp + i * d * np.sin(a))
        #set the last point on convergence part. Hold the exact location of the
        # last point.
        l2_x.append(x_exp - l2 * np.cos(a))
        l2_y.append(y_exp + l2 * np.sin(a))

        x0 = l2_x[n_l2]
        y0 = l2_y[n_l2]
        # coordinate of the first point on the rounding circular arc between
        # convergence and parallel part.
        b = 0.5 * np.pi - a
        # calculate the angle between horizontal axis and the line, which
        # between circular center and the first point on circular arc.
        y_r = y0 - r * np.sin(b)  # y coordinate of center of the circular arc.
        x_r = x0 - r * np.cos(b)  # x coordinate of center of the circular arc.
        r_x = list()  # lists for coordinates of points on the circular arc.
        r_y = list()
        x_r_end = x_r
        y_r_end = y_r + r
        # the last point of the circular arc/
        # the first point of the parallel part.
        r_l = r * a # length of the circular arc.
        n = math.ceil(r_l / d)  # circular arc can be partitioned in n parts.
        a_step = a / n  # the angle for every part of circular arc.

        for i in range(n-1):
            # set coordinate for every point on the circular arc.
            r_x.append(x_r + r * np.cos(b + (i + 1) * a_step))
            r_y.append(y_r + r * np.sin(b + (i + 1) * a_step))

        l1_x = list()
        l1_y = list()
        #set lists for coordinates of points on parallel part.
        l1_x.append(x_r_end)  # set the first point on the parallel part.
        l1_y.append(y_r_end)

        i2 = 0
        x_end = x_r_end - l1
        # get the x-coordinate of the last point on the parallel part.
        while l1_x[i2] > x_end:
            # set the coorinate for every point on the parallel part.
            l1_x.append(l1_x[i2] - d)
            l1_y.append(l1_y[i2])
            i2 = i2 + 1

        x_chamber = l2_x + r_x + l1_x
        y_chamber = l2_y + r_y + l1_y
        # assemble all points on outer contour of combustion chamber.
        x_chamber.reverse()
        y_chamber.reverse()
        # reverse the order of all points. Because the points should be
        # outputted with increasing x-coordinate.
        return x_chamber, y_chamber



    def lin_chamber_volume(self, x_in_c, y_in_c, x_out_c, y_out_c, depth):
        '''calculte the volume of a linear aerospike combustion chamber.
                Args:
                    x_in_c, y_in_c: lists for x- and y-coordinates of
                                    points on inner contour of
                                   combustion chamber.
                    x_out_c, y_out_c: lists for x- and y-coordinates of
                                      points on outer contour of
                                      combustion chamber.
                    depth: depth of combustion chamber.
                Returns:
                    volume: volume of the combustion chamber.
               '''
        n = len(x_in_c)
        #get quantity of the points on the contour of combustion chamber.
        volume = 0 #set the total volume to 0.
        for i in range(n-1):
            area_aver = 0.5 * (np.sqrt((x_in_c[i] - x_out_c[i]) ** 2 +
                                       (y_in_c[i] -y_out_c[i]) ** 2) +
                               np.sqrt((x_in_c[i+1] - x_out_c[i+1]) ** 2 +
                               (y_in_c[i+1] -y_out_c[i+1]) ** 2)) * depth
            # calculate average area of every slice of the combustion chamber.

            l_aver = 0.5 * (np.sqrt((x_in_c[i] - x_in_c[i+1]) ** 2 +
                                    (y_in_c[i] -y_in_c[i+1]) ** 2) +
                            np.sqrt((x_out_c[i] - x_out_c[i+1]) ** 2 +
                                    (y_out_c[i] -y_out_c[i+1]) ** 2))
            # calculate average width of every slice of the combustion chamber.

            d_v = area_aver * l_aver * 2
            # get the volume of every slice.

            volume = volume + d_v
            # add the volume of slice in total volume.

        return volume



    def anu_chamber_volume(self, x_in_c, y_in_c, x_out_c, y_out_c):
        '''calculte the volume of a annular aerospike combustion chamber.
                Args:
                    x_in_c, y_in_c: lists for x- and y-coordinates of
                                    points on inner contour of
                                    combustion chamber.
                    x_out_c, y_out_c: lists for x- and y-coordinates of
                                      points on outer contour of
                                      combustion chamber.
                Returns:
                    volume: volume of the combustion chamber.
                       '''
        n = len(x_in_c)
        # get quantity of the points on the contour of combustion chamber.
        volume = 0 #set the total volume to 0.
        for i in range(n - 1):
            area_aver = 0.5 * (self.supfic_area_cone(x_in_c[i],
                                                     x_out_c[i], y_in_c[i],
                                                     y_out_c[i]) +
                               self.supfic_area_cone(x_in_c[i+1], x_out_c[i+1],
                                                     y_in_c[i+1], y_out_c[i+1])
                               )
            # calculate average area of every slice of the combustion chamber.

            l_aver = 0.5 * (np.sqrt((x_in_c[i] - x_in_c[i + 1]) ** 2 +
                                    (y_in_c[i] - y_in_c[i + 1]) ** 2) +
                            np.sqrt(
                (x_out_c[i] - x_out_c[i + 1]) ** 2 +
                (y_out_c[i] - y_out_c[i + 1]) ** 2))
            # calculate average width of every slice of the combustion chamber.

            d_v = area_aver * l_aver
            # get the volume of every slice.

            volume = volume + d_v
            # add the volume of slice in total volume.

        return volume




class Gasdynamic():

    def __init__(self, k):
        '''get the isentropic exponent.
                Args:
                    k: isentropic exponent.
                    '''
        self.k = k # set the isentropic exponent.



    def machnr(self, MachL, MachR, Aver):
        '''calculate Mach number and critical Mach number through
           cross-sectional area.
                Args:
                    MachL: left limit of the range of critical Mach number,
                           in which the desired Mach number lies.
                    MachR: right limit of the range of critical Mach number,
                           in which the desired Mach number lies.
                    Aver: list for relationships of cross-sectional area(A*/A).
                Returns:
                    Ma: Mach number.
                    '''
        k = self.k # get the isentropic exponent.
        Ma = list() # set list for Mach number.
        CMa = list() # set list for critical Mach number.
        m = len(Aver)
        # get quantity of the elements in list for area-relationship(A*/A).

        for iMa in range(m):
            # calculate the corresponding Mach number for every A*/A.
            MaL = MachL # set left limit for calculation.
            MaR = MachR # set right limit for calculation.
            Mak = 0.5 * (MaR + MaL)
            # calculate the first value of critical Mach number.
            Av_cal = Mak * (((1 - ((k - 1) / (k + 1)) * Mak ** 2) /
                             (1 - ((k - 1) / (k + 1)))) ** (1 / (k - 1)))
            # calculate the area-relationship with the value Mak.
            if Mak > 1:# if in supersonic area.

                while (np.absolute(Av_cal - Aver[iMa]) > 0.0000001):
                    #calculate the area-relationship, until the difference
                    #between the actual value and calculated value is smaller
                    #than 0.0000001.
                    if Av_cal - Aver[iMa] < 0:
                        # if the calculated value is smaller than the
                        # actual value.
                        MaR = Mak# actual Mak instead of right limit.
                    else:
                        # if the calculated value is bigger than the
                        # actual value.
                        MaL = Mak#actual Mak instead of left limit.
                    Mak = 0.5 * (MaR + MaL)
                    #calculate the new Mak with new limits
                    Av_cal = Mak * (((1 - ((k - 1) / (k + 1)) * Mak ** 2) /
                                (1 - ((k - 1) / (k + 1)))) ** (1 / (k - 1)))
                    #calculate new A*/A with new Mak.
                CMa.append(Mak)
                #set the calculated Mak value for actual A*/A in the list
                # for critical Mach number.
            else: # if in subsonic area.

                while (np.absolute(Av_cal - Aver[iMa]) > 0.0000001):
                    # calculate the area-relationship, until the difference
                    # between the actual value and calculated value is smaller
                    # than 0.0000001.
                    if Av_cal - Aver[iMa] < 0:
                        # if the calculated value is smaller than the
                        # actual value.
                        MaL = Mak# actual Mak instead of left limit.
                    else:
                        # if the calculated value is bigger than the
                        # actual value.
                        MaR = Mak# actual Mak instead of right limit.
                    Mak = 0.5 * (MaR + MaL)
                    # calculate the new Mak with new limits
                    Av_cal = Mak * (((1 - ((k - 1) / (k + 1)) * Mak ** 2) /
                                (1 - ((k - 1) / (k + 1)))) ** (1 / (k - 1)))
                    # calculate new A*/A with new Mak.
                CMa.append(Mak)
                # set the calculated Mak value for actual A*/A in the list
                # for critical Mach number.

        for iMa in range(m):
            # translate all critical Mach numbers to Mach number.
            Ma.append(np.sqrt(CMa[iMa] ** 2 / (1 - 0.5 * (k - 1) *
                                               (CMa[iMa] ** 2 - 1))))

        return Ma



    def ma_a(self, A_ratio):
        '''calculate Mach number through area-relationship(A*/A)
                Args:
                    A_ratio: list for area-relationship.
                Returns:
                    Mach: Mach number for every value of A*/A.
                    '''
        k = self.k #get the isentropic exponent.
        CMamax = np.sqrt((k + 1) / (k - 1))  # maximal critical Mach-number.
        A_ratio_max = max(A_ratio) # find throat point.
        n = A_ratio.index(A_ratio_max) # index throat point.
        m = len(A_ratio) # get quantity of area-relationship.
        A_ratio_cham=list() # set list for values A*/A before throat.
        A_ratio_nozzle=list() # set list for values A*/A after throat.
        if n == 0: # if there is no A*/A before throat.
            MachCham = list()
            # set the list for Mach number in subsonic area as empty.
        else:# if there are values of A*/A before throat.
            for i in range(n):
                #set a list for area-relationship in subsonic.
                A_ratio_cham.append(A_ratio[i])
            MachCham = self.machnr(0, 1, A_ratio_cham)
            # calculate Mach number in subsonic area.
        if m == n:# if there is no A*/A after throat.
            MachNozzle = list()
            #set the list for Mach number in supersonic area as empty.
        else:# if there are values of A*/A after throat.
            for i in range(n + 1, m):
                # set a list for area-relationship in subsonic.
                A_ratio_nozzle.append(A_ratio[i])
            MachNozzle = self.machnr(1, CMamax, A_ratio_nozzle)
            # calculate Mach number in supersonic area.
        Mach = MachCham + [1] + MachNozzle # assemble all Mach numbers.

        return Mach



    def t_ma(self, Tc, Ma):
        '''calculate the Temperatur through Mach number.
        Args:
            Tc: Temperature in Chamber.
            Ma: list of Mach number along the Thrust-Chamber.
        Returns:
            Tem: list of Temperature along the Thrust-Chamber.
        '''
        k = self.k # get the isentropic exponent.
        Tem = list() # set list for Temperature.
        n = len(Ma) # get quantity of Mach numbers.

        for i in range(n):
            # calculate temperature for every Mach number and set the value
            # in the list for temperature.
            Tem.append(Tc / (1 + (0.5 * (k - 1) * Ma[i] ** 2)))

        return Tem



    def p_ma(self, Pc, Ma):
        '''calculate the Pressure through Mach number.
                Args:
                    Pc: Pressure in Chamber.
                    Ma: list of Mach number along the Thrust-Chamber.
                Returns:
                    P: list of Pressure along the Thrust-Chamber.
                '''
        k = self.k # get the isentropic exponent.
        P = list() # set list for Pressure.
        n = len(Ma) # get quantity of Mach numbers.

        for i in range(n):
            # calculate pressure for every Mach number and set the value
            # in the list for pressure.
            P.append(Pc / ((1 + ((k - 1) / 2) * Ma[i] ** 2) ** (k / (k - 1))))

        return P



    def rho_ma(self, roh_c, Ma):
        '''calculate the density through Mach number.
                Args:
                   roh_c: density in Chamber.
                   Ma: list of Mach number along the Thrust-Chamber.
                Returns:
                   roh: list of density along the Thrust-Chamber.
                        '''
        k = self.k # get the isntropic exponent.
        roh = list() # set list for density.
        n = len(Ma) # get quantity of Mach numbers.
        for i in range(n):
            roh.append(roh_c/((1+((Ma[i]**2)*(k-1)/2))**(1/(k-1))))
            # calculate density for every Mach number and set the value
            # in the list for density.

        return roh



    def ma_t(self, Tc, T):
        '''calculate the Mach number through Temperatures.
                Args:
                    Tc: Tempereature in Chamber.
                    T: list of Temperature along the Thrust-Chamber.
                Returns:
                    Ma: list of Mach number along the Thrust-Chamber.
                                '''
        k = self.k# get the isntropic exponent.
        Ma = list()# set list for Mach number.
        n = len(T)# get quantity of Temperature values.
        for i in range(n):
            Ma.append(np.sqrt((2/(k-1))*((Tc/T[i])-1)))
            # calculate Mach number for every Temperature value and set the
            # value in the list for Mach number.

        return Ma



    def ma_p(self, Pc, P):
        '''calculate the Mach number through Pressure.
                Args:
                    Pc: Pressure in Chamber.
                    P: list of Pressure along the Thrust-Chamber.
                Returns:
                    Ma: list of Mach number along the Thrust-Chamber.
                                        '''
        k = self.k# get the isntropic exponent.
        Ma = list()# set list for Mach number.
        n = len(P)# get quantity of Pressure values.
        for i in range(n):
            Ma.append(np.sqrt((2/(k-1))*(((Pc/P[i])**((k-1)/k))-1)))
            # calculate Mach number for every pressure value and set the value
            # in the list for Mach number.

        return Ma



    def ma_roh(self, roh_c, roh):
        '''calculate the Mach number through density.
                Args:
                  roh_c: density in Chamber.
                  roh: list of Pressure along the Thrust-Chamber.
                Returns:
                  Ma: list of Mach number along the Thrust-Chamber.
                                                '''
        k = self.k# get the isntropic exponent.
        Ma = list()# set list for Mach number.
        n = len(roh)# get quantity of density value.
        for i in range(n):
            Ma.append(np.sqrt((2/(k-1))*(((roh_c/roh[i])**(k-1))-1)))
            # calculate Mach number for every density value and set the value
            # in the list for Mach number.

        return Ma



    def charact_velocity(self, R, Tc):
        '''calculate the characteristic velocity.
                Args:
                    R: specific gas constant.
                    Tc: Temperature in Chamber.
                Returns:
                    c_char: characteristic velocity.
                    '''
        k = self.k# get the isentopic exponent.
        c_char =                                                              \
            (np.sqrt(R * Tc)) /                                               \
        (np.sqrt(k * ((2 / (k + 1)) ** ((k + 1) / (k - 1)))))
        #calculate the characteristic velocity.

        return c_char


class Heat_transfer():

    def radiation(self, h2o, co2, p, d, T, Tw):
        '''calculate heat transfer through Radiation.
                Args:
                    h2o: mole fraction of h2o.
                    co2: mole fraction of co2.
                    p: pressure.
                    d: hydraulic diameter.
                    T: gas temperature.
                    Tw: wall temperature.
                Returns:
                    qr: heat flux density through radiation.
                    '''
        qh2o = 4.07 * ((1.02 * (h2o * p / 100000)) ** 0.8) *                  \
        ((0.9 * d) ** 0.6) * ((T / 100) ** 3 - (Tw / 100) ** 3)
        #calculate the radiation through h2o.
        qco2 = 4.07 * ((1.02 * (0.9 * d * co2 * p / 100000)) ** (1 / 3)) *    \
        ((T / 100) ** 3.5 - (Tw / 100) ** 3.5)
        #calculate the radiation through co2.

        qr = qco2 + qh2o
        # assemble the radiation.

        return qr



    def convective(self,c, dt, r_c, viscosity, heatcp, prandtl, champressure,
                   charaC, areaRela, Tc, T_re, Twall, MachNr, kappa):
        '''calculate the convective heat transfer.
                Args:
                    c: constant in Bartz equation.
                    dt: hydraulic diameter of cross-sectional area at throat.
                    r_c: rounding radius at throat.
                    viscosity: gas viscosity.
                    heatcp: heat capacity by constant pressure.
                    prandtl: prandtl number of gas.
                    champressure: pressure in chamber.
                    charaC: characteristic velocity.
                    areaRela: area-relationship(A*/A).
                    Tc: temperature in chamber.
                    T_re: recovery temperature of gas.
                    Twall: wall temperature.
                    MachNr: Mach number.
                    kappa: isentropic exponent.
                Returns:
                    heattrans_coefficient: convective heat transfer
                                           coefficient.
                    heat_flux_density: convective heat flux density.
                    '''
        sig1 = (Twall / Tc)
        sig2 = 1 + 0.5 * (kappa - 1) * MachNr ** 2

        sigma =1 /                                                            \
        (((0.5 * sig1 * sig2 + 0.5) ** (0.68)) * (sig2 ** 0.12))

        heat1 = c / (dt ** 0.2)
        heat2 = (heatcp * (viscosity ** 0.2)) / (prandtl ** 0.6)
        heat3 = (champressure / charaC) ** 0.8
        heat4 = (dt / r_c) ** 0.1

        heattrans_coeffcient = heat1 * heat2 * heat3 * heat4                  \
         * (areaRela ** 0.9) * sigma

        # calculate the convective heat transfer coefficient through
        # Bartz-equation.
        heat_flux_density = heattrans_coeffcient * (T_re - Twall)
        # calculate the convective heat flux density.

        return heattrans_coeffcient, heat_flux_density



class Thermo_load():

    def __init__(self, Pc, Tc, Tw):
        '''get three gas dynamic values.
                Args:
                    Pc: Pressure in chamber.
                    Tc: Temperature in chamber.
                    Tw: wall Temperature.
                    '''
        self.Pc = Pc
        self.Tc = Tc
        self.Tw = Tw
        # get the three values.



    def lin_spike(self, x_spike, y_spike, x_exp, y_exp, depth, r_throat, k,
                  viscosity, cp2, cp1, cp0, R, prandtl, h2o, co2):
        '''calculate the heat transfer on linear aerospike.
                Args:
                    x_spike, y_spike: x- and y-coordinates of the points
                                      on aerospike.
                    x_exp, y_exp: x- and y-coordinates of the expansion point.
                    depth: depth of the linear aerospike.
                    r_throat: rounding radius at throat.
                    k: isentropic exponent.
                    viscosity: gas viscosity.
                    cp0, cp1, cp2: parameters of the 2 degreee polynomial for
                                   heat capacity of gas.
                    R: specific gas constant.
                    prandtl: prandtl number of gas.
                    h2o: mole fraction of h2o.
                    co2: mole fraction of co2.
                Returns:
                    flux_conv: convective heat flux density.
                    flux_radi: radiational heat flux denstiy.
                    flox_together: total heat flux density.
                    heat_spike_together: total heat flux power around
                                         the aerospike.
                    heat_trans_coeff: covective heat transfer coefficient.
                                         '''
        Geo = Aerospike_Geometry(0)#use class 'Aerospike_Geometry' for linear.
        area_ratio = Geo.spike_area_ratio(x_exp, y_exp, x_spike, y_spike)
        # get the area-relationship along the aerospike.
        gasdy = Gasdynamic(k)# use class 'Gasdynamic' with isentropic exponent.
        mach = gasdy.ma_a(area_ratio)#get the Mach number along the aerospike.
        dt = list()# set list for hydraulic diameter.
        n = len(x_spike)#get quantity of points on aerospike.
        for i in range(n):
            #calculate hydraulic diameter for every position on the aerospike.
            dt.append(Geo.lin_hydraulic_diameter(x_spike[i],
                                                 y_spike[i],
                                                 x_exp, y_exp, depth))
        if self.dt_throat is None:
            # if the hydraulic diameter at throat is not given.
            self.dt_throat = dt[0]
            # set the firs value of hydraulic diameter as hydraulic diameter
            # at throat.
        v_charact = gasdy.charact_velocity(R, self.Tc)
        #get characteristic velocity.
        T = gasdy.t_ma(self.Tc, mach)#get temperature at every positions.
        P = gasdy.p_ma(self.Pc, mach)#get pressure at every position.
        r_recovery = prandtl ** 0.33 # recovery factor for turbulent flow
        T_re = list() # set list for recovery temperature of gas.
        heatcp = list()# set list for cp.
        for i in range(n):
            #calculate recovery temperature and cp value at every position.
            heatcp.append(cp2 * T[i] ** 2 + cp1 * T[i] + cp0)
            T_re.append(r_recovery * (self.Tc - T[i]) + T[i])
        heat_trans = Heat_transfer()#use class 'Heat_transfer'.
        heat_trans_coeff = list()
        #set list for convective heat transfer coeffizient.
        flux_conv = list()#set list for convective heat flux density.
        flux_radi = list()#set list for raditional heat flux density.
        flux_together = list()#set list for total heat flux density.
        heat_spike_together = 0#set the total heat flux power as 0.
        for i in range(n):
            #get the convective heat transfer coeffizient,
            # convective heat flux density, total heat flux density
            # at every position and set they in the corresponding lists.
            # for linear aerospike set constant in Bartz equation c=0.035
            coeffcient, conv_flux_density = heat_trans.convective(0.032,
                self.dt_throat, r_throat, viscosity, heatcp[i], prandtl,
                self.Pc, v_charact, area_ratio[i], self.Tc, T_re[i], self.Tw,
                                                                  mach[i], k)
            heat_trans_coeff.append(coeffcient)
            flux_conv.append(conv_flux_density)

            flux_radi.append(heat_trans.radiation(h2o, co2, P[i],
                                                  dt[i], T[i], self.Tw))
            #get raditional heat flux density at every position and set it in
            # the list for raditional heat flux density.
            flux_together.append(flux_conv[i] + flux_radi[i])
            #calculate the total heat flux density at every postion, and set
            #it in corresponding list.

        for i in range(n-1):
            #calculate the total heat flux power.
            spike_area = depth * np.sqrt((x_spike[i]-x_spike[i+1])**2+
                                         (y_spike[i]-y_spike[i+1])**2)
            #calculate the superficial area of a slice.

            flux_average = (flux_together[i]+flux_together[i+1])/2
            #calculate average heat flux density on the
            # superficial area of the slice.
            heat_spike = spike_area * flux_average * 2
            # '*2' because there just half part is calculated.
            #calculate the heat flux power on the slice.
            heat_spike_together = heat_spike_together + heat_spike
            #add the heat flux power of the slice in total heat flux power.

        return flux_conv, flux_radi, flux_together,                           \
        heat_spike_together, heat_trans_coeff



    def anu_spike(self, x_spike, y_spike, x_exp, y_exp, r_throat, k,
                  viscosity, cp2, cp1, cp0, R, prandtl, h2o, co2):
        '''calculate the heat transfer on annular aerospike.
                Args:
                    x_spike, y_spike: x- and y-coordinates of the points
                                      on aerospike.
                    x_exp, y_exp: x- and y-coordinates of the expansion point.
                    r_throat: rounding radius at throat.
                    k: isentropic exponent.
                    viscosity: gas viscosity.
                    cp0, cp1, cp2: parameters of the 2 degreee polynomial for
                                   heat capacity of gas.
                    R: specific gas constant.
                    prandtl: prandtl number of gas.
                    h2o: mole fraction of h2o.
                    co2: mole fraction of co2.
                Returns:
                    flux_conv: convective heat flux density.
                    flux_radi: radiational heat flux denstiy.
                    flox_together: total heat flux density.
                    heat_spike_together: total heat flux power around
                                         the aerospike.
                    heat_trans_coeff: covective heat transfer coefficient.
                                                 '''
        Geo = Aerospike_Geometry(1)#use class 'Aerospike_Geometry' for annular.
        area_ratio = Geo.spike_area_ratio(x_exp, y_exp, x_spike, y_spike)
        # get the area-relationship along the aerospike.
        gasdy = Gasdynamic(k)#use class 'Gasdynamic' with isentropic exponent.
        mach = gasdy.ma_a(area_ratio)#get the Mach number along the aerospike.
        dt = list()# set list for hydraulic diameter.
        n = len(area_ratio)#get quantity of points on aerospike.
        for i in range(n):
            # calculate hydraulic diameter for every position on the aerospike.
            dt.append(Geo.anu_hydraulic_diameter(x_spike[i], y_spike[i],
                                                 x_exp, y_exp))

        if self.dt_throat is None:
            # if the hydraulic diameter at throat is not given.
            self.dt_throat = dt[0]
            # set the firs value of hydraulic diameter as hydraulic diameter
            # at throat.
        v_charact = gasdy.charact_velocity(R, self.Tc)
        # get characteristic velocity.
        T = gasdy.t_ma(self.Tc, mach)#get temperature at every positions.
        P = gasdy.p_ma(self.Pc, mach)#get pressure at every position.
        r_recovery = prandtl ** 0.33  # recovery factor for turbulent flow.
        T_re = list()# set list for recovery temperature of gas.
        heatcp = list()# set list for cp.
        for i in range(n):
            # calculate recovery temperature and cp value at every position.
            heatcp.append(cp2 * T[i] ** 2 + cp1 * T[i] + cp0)
            T_re.append(r_recovery * (self.Tc - T[i]) + T[i])

        heat_trans = Heat_transfer()#use class 'Heat_transfer'.
        heat_trans_coeff = list()
        # set list for convective heat transfer coeffizient.
        flux_conv = list()#set list for convective heat flux density.
        flux_radi = list()#set list for raditional heat flux density.
        flux_together = list()#set list for total heat flux density.
        heat_spike_together = 0#set the total heat flux power as 0.
        for i in range(n):
            # get the convective heat transfer coeffizient,
            # convective heat flux density, total heat flux density
            # at every position and set they in the corresponding lists.
            # for annular aerospike set constant in Bartz equation c=0.046
            coeffcient, conv_flux_density = heat_trans.convective(0.046,
                self.dt_throat, r_throat, viscosity, heatcp[i], prandtl,
                self.Pc, v_charact, area_ratio[i], self.Tc, T_re[i], self.Tw,
                mach[i], k)

            heat_trans_coeff.append(coeffcient)
            flux_conv.append(conv_flux_density)

            flux_radi.append(heat_trans.radiation(h2o, co2, P[i], dt[i],
                                                  T[i], self.Tw))
            # get raditional heat flux density at every position and set it in
            # the list for raditional heat flux density.
            flux_together.append(flux_conv[i] + flux_radi[i])
            # calculate the total heat flux density at every postion, and set
            # it in corresponding list.

        for i in range(n-1):
            # calculate the total heat flux power.
            spike_area = Geo.supfic_area_cone(x_spike[i], x_spike[i+1],
                                              y_spike[i], y_spike[i+1])
            # calculate the superficial area of a slice.
            flux_average = (flux_together[i]+flux_together[i+1])/2
            # calculate average heat flux density on the
            # superficial area of the slice.
            heat_spike = spike_area * flux_average
            # calculate the heat flux power on the slice.
            heat_spike_together = heat_spike_together + heat_spike
            # add the heat flux power of the slice in total heat flux power.

        return flux_conv, flux_radi, flux_together, heat_spike_together,      \
        heat_trans_coeff



    def lin_chamber(self, x_in_c, y_in_c, x_out_c, y_out_c, x_throat,
                    y_throat, x_exp, y_exp, depth, r_throat, k, viscosity,
                    cp2, cp1, cp0, R, prandtl, h2o, co2):
        '''calculate the heat transfer on linear aerospike combustion chamber.
                Args:
                    x_in_c, y_in_c: x- and y-coordinates of the points on
                                    inner contour of combustion chamber.
                    x_out_c, y_out_c: x- and y-coordinates of the points on
                                      outer contour of combustion chamber.
                    x_throat, y_throat: x- and y-coordinates of the
                                        throat point.
                    x_exp, y_exp: x- and y-coordinates of the expansion point.
                    depth: depth of the linear aerospike.
                    r_throat: rounding radius at throat.
                    k: isentropic exponent.
                    viscosity: gas viscosity.
                    cp0, cp1, cp2: parameters of the 2 degreee polynomial for
                                   heat capacity of gas.
                    R: specific gas constant.
                    prandtl: prandtl number of gas.
                    h2o: mole fraction of h2o.
                    co2: mole fraction of co2.
                Returns:
                    flux_conv: convective heat flux density.
                    flux_radi: radiational heat flux denstiy.
                    flox_together: total heat flux density.
                    heat_in_together: total heat flux power around the inner
                                      contour of combustion chamber.
                    heat_out_together: total heat flux power around the outer
                                       contour of combustion chamber.
                    heat_trans_coeff: covective heat transfer coefficient.
                                                '''
        Geo = Aerospike_Geometry(0)#use class 'Aerospike_Geometry' for linear.
        area_ratio = Geo.chamber_area_ratio(x_in_c, y_in_c, x_out_c, y_out_c,
                                            x_throat, y_throat, x_exp, y_exp)
        # get the area-relationship along the combustion chamber.
        n = len(area_ratio)
        #get quantity of points on contour of combustion chamber.
        gasdy = Gasdynamic(k)# use class 'Gasdynamic' with isentropic exponent.
        mach = gasdy.ma_a(area_ratio)#get the Mach number along the chamber.
        dt = list()# set list for hydraulic diameter.
        for i in range(n):
            # calculate hydraulic diameter for every position on the chamber.

            dt.append(Geo.lin_hydraulic_diameter(x_in_c[i], y_in_c[i],
                                                 x_out_c[i], y_out_c[i],
                                                 depth))

        self.dt_throat = dt[n-1]#set the value of hydraulic diameter at throat.
        v_charact = gasdy.charact_velocity(R, self.Tc)
        # get characteristic velocity.
        T = gasdy.t_ma(self.Tc, mach)#get temperature at every positions.
        P = gasdy.p_ma(self.Pc, mach)#get pressure at every position.
        r_recovery = prandtl ** 0.33  # recovery factor for turbulent flow.
        T_re = list() # set list for recovery temperature of gas.
        heatcp = list()# set list for cp.
        for i in range(n):
            # calculate recovery temperature and cp value at every position.
            heatcp.append(cp2 * T[i] ** 2 + cp1 * T[i] + cp0)
            T_re.append(r_recovery * (self.Tc - T[i]) + T[i])

        heat_trans = Heat_transfer()#use class 'Heat_transfer'.
        heat_trans_coeff = list()
        heat_trans_coeff_out = list()
        # set list for convective heat transfer coeffizient.
        flux_conv = list()#set list for convective heat flux density.
        flux_conv_out = list()
        flux_radi = list()#set list for raditional heat flux density.
        flux_together = list()#set list for total heat flux density.
        flux_together_out = list()
        heat_out_together = 0
        #set the total heat flux on outer contour power as 0.
        heat_in_together = 0
        # set the total heat flux on outer contour power as 0.
        for i in range(n):
            # get the convective heat transfer coeffizient,
            # convective heat flux density, total heat flux density
            # at every position and set they in the corresponding lists.
            # for linear aerospike set constant in Bartz equation c=0.035
            coeffcient, conv_flux_density = heat_trans.convective(0.032,
                self.dt_throat, r_throat, viscosity, heatcp[i], prandtl,
                self.Pc, v_charact, area_ratio[i], self.Tc, T_re[i], self.Tw,
                                                                  mach[i], k)

            coeffcient_out = coeffcient * (0.02/0.032)
            conv_flux_density_out = conv_flux_density * (0.02 / 0.032)
            # calculate heat transfer coeffizient for outer shell
            # with C = 0,021

            heat_trans_coeff.append(coeffcient)
            flux_conv.append(conv_flux_density)
            heat_trans_coeff_out.append(coeffcient_out)
            flux_conv_out.append(conv_flux_density_out)

            flux_radi.append(heat_trans.radiation(h2o, co2, P[i], dt[i],
                                                  T[i], self.Tw))
            # get raditional heat flux density at every position and set it in
            # the list for raditional heat flux density.
            flux_together.append(flux_conv[i] + flux_radi[i])
            flux_together_out.append(flux_conv_out[i] + flux_radi[i])
            # calculate the total heat flux density at every postion, and set
            # it in corresponding list.
        for i in range(n-1):
            # calculate the total heat flux power.
            super_area_out = depth * np.sqrt((x_out_c[i]-x_out_c[i+1])**2+
                                             (y_out_c[i]-y_out_c[i+1])**2)
            # calculate the outer superficial area of a slice.

            super_area_in = depth * np.sqrt((x_in_c[i]-x_in_c[i+1])**2+
                                            (y_in_c[i]-y_in_c[i+1])**2)
            # calculate the inner superficial area of a slice.

            flux_average = (flux_together[i] + flux_together[i+1])/2
            flux_average_out = (flux_together_out[i] +
                                flux_together_out[i + 1]) / 2
            # calculate average heat flux density on the
            # superficial area of the slice.
            heat_out = flux_average_out * super_area_out * 2
            heat_in = flux_average * super_area_in * 2
            # '*2' because there just half part is calculated.
            # calculate the heat flux power on the outer and inner slice.
            heat_out_together = heat_out_together + heat_out
            heat_in_together = heat_in_together + heat_in
            # add the heat flux power of the slice in total heat flux power.

        return flux_conv, flux_conv_out, flux_radi, flux_together,            \
        flux_together_out, heat_in_together, heat_out_together,               \
        heat_trans_coeff,                                                     \
        heat_trans_coeff_out



    def anu_chamber(self, x_in_c, y_in_c, x_out_c, y_out_c, x_throat,
                    y_throat, x_exp, y_exp, r_throat, k, viscosity, cp2,
                    cp1, cp0, R, prandtl, h2o, co2):
        '''calculate the heat transfer on linear aerospike combustion chamber.
                Args:
                    x_in_c, y_in_c: x- and y-coordinates of the points on
                                    inner contour of combustion chamber.
                    x_out_c, y_out_c: x- and y-coordinates of the points on
                                      outer contour of combustion chamber.
                    x_throat, y_throat: x- and y-coordinates of the
                                        throat point.
                    x_exp, y_exp: x- and y-coordinates of the expansion point.
                    r_throat: rounding radius at throat.
                    k: isentropic exponent.
                    viscosity: gas viscosity.
                    cp0, cp1, cp2: parameters of the 2 degreee polynomial for
                                   heat capacity of gas.
                    R: specific gas constant.
                    prandtl: prandtl number of gas.
                    h2o: mole fraction of h2o.
                    co2: mole fraction of co2.
                Returns:
                    flux_conv: convective heat flux density.
                    flux_radi: radiational heat flux denstiy.
                    flox_together: total heat flux density.
                    heat_in_together: total heat flux power around the inner
                                      contour of combustion chamber.
                    heat_out_together: total heat flux power around the outer
                                       contour of combustion chamber.
                    heat_trans_coeff: covective heat transfer coefficient.
                                                        '''
        Geo = Aerospike_Geometry(1)#use class 'Aerospike_Geometry' for annular.
        area_ratio = Geo.chamber_area_ratio(x_in_c, y_in_c, x_out_c, y_out_c,
                                            x_throat, y_throat, x_exp, y_exp)
        # get the area-relationship along the combustion chamber.
        n = len(area_ratio)
        # get quantity of points on contour of combustion chamber.
        gasdy = Gasdynamic(k)# use class 'Gasdynamic' with isentropic exponent.
        mach = gasdy.ma_a(area_ratio)#get the Mach number along the chamber.
        dt = list()# set list for hydraulic diameter.
        for i in range(n):
            # calculate hydraulic diameter for every position on the chamber.
            dt.append(Geo.anu_hydraulic_diameter(x_in_c[i], y_in_c[i],
                                                 x_out_c[i], y_out_c[i]))
        self.dt_throat = dt[n - 1]
        v_charact = gasdy.charact_velocity(R, self.Tc)
        # get characteristic velocity.
        T = gasdy.t_ma(self.Tc, mach)#get temperature at every positions.
        P = gasdy.p_ma(self.Pc, mach)#get pressure at every position.
        r_recovery = prandtl ** 0.33  # recovery factor for turbulent flow.
        T_re = list()  # list for recovery temperature of gas
        heatcp = list()# set list for cp.
        for i in range(n):
            # calculate recovery temperature and cp value at every position.
            heatcp.append(cp2 * T[i] ** 2 + cp1 * T[i] + cp0)
            T_re.append(r_recovery * (self.Tc - T[i]) + T[i])

        heat_trans = Heat_transfer()#use class 'Heat_transfer'.
        heat_trans_coeff = list()
        heat_trans_coeff_out = list()
        # set list for convective heat transfer coeffizient.
        flux_conv = list()#set list for convective heat flux density.
        flux_conv_out = list()
        flux_radi = list()#set list for raditional heat flux density.
        flux_together = list()#set list for total heat flux density.
        flux_together_out = list()
        heat_out_together = 0
        # set the total heat flux on outer contour power as 0.
        heat_in_together = 0
        # set the total heat flux on outer contour power as 0.
        for i in range(n):
            # get the convective heat transfer coeffizient,
            # convective heat flux density, total heat flux density
            # at every position and set they in the corresponding lists.
            # for annular aerospike set constant in Bartz equation c=0.046
            coeffcient, conv_flux_density = heat_trans.convective(0.046,
                self.dt_throat, r_throat, viscosity, heatcp[i], prandtl,
                self.Pc, v_charact, area_ratio[i], self.Tc, T_re[i], self.Tw,
                                                                  mach[i], k)

            coeffcient_out = coeffcient * (0.026 / 0.046)
            conv_flux_density_out = conv_flux_density * (0.026 / 0.046)
            # calculate heat transfer coeffizient for outer shell
            # with C = 0,026

            heat_trans_coeff.append(coeffcient)
            flux_conv.append(conv_flux_density)
            heat_trans_coeff_out.append(coeffcient_out)
            flux_conv_out.append(conv_flux_density_out)

            flux_radi.append(heat_trans.radiation(h2o, co2, P[i], dt[i],
                                                  T[i], self.Tw))
            # get raditional heat flux density at every position and set it in
            # the list for raditional heat flux density.
            flux_together.append(flux_conv[i] + flux_radi[i])
            flux_together_out.append(flux_conv_out[i] + flux_radi[i])
            # calculate the total heat flux density at every postion, and set
            # it in corresponding list.

        for i in range(n-1):
            # calculate the total heat flux power.
            super_area_out = Geo.supfic_area_cone(x_out_c[i], x_out_c[i+1],
                                                  y_out_c[i], y_out_c[i+1])
            # calculate the outer superficial area of a slice.

            super_area_in = Geo.supfic_area_cone(x_in_c[i], x_in_c[i+1],
                                                 y_in_c[i], y_in_c[i+1])
            # calculate the inner superficial area of a slice.

            flux_average = (flux_together[i] + flux_together[i+1])/2
            flux_average_out = (flux_together_out[i] +
                                flux_together_out[i + 1]) / 2
            # calculate average heat flux density on the
            # superficial area of the slice.
            heat_out = flux_average_out * super_area_out
            heat_in = flux_average * super_area_in
            # calculate the heat flux power on the outer and inner slice.
            heat_out_together = heat_out_together + heat_out
            heat_in_together = heat_in_together + heat_in
            # add the heat flux power of the slice in total heat flux power.

        return flux_conv, flux_conv_out, flux_radi, flux_together,            \
        flux_together_out ,heat_in_together,                                  \
        heat_out_together, heat_trans_coeff, heat_trans_coeff_out



