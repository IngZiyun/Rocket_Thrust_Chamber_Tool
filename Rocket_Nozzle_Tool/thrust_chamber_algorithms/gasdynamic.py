import numpy as np

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

    def critical_mach_from_mach(self, ma):

        ma_star = np.sqrt(
            (ma ** 2) / (1 + (((self.k - 1) / (self.k + 1)) * (ma ** 2 - 1))))

        return ma_star

    def mach_from_critical_mach(self, ma_star):

        ma = np.sqrt(
            (ma_star ** 2) / (1 - (((self.k - 1) * (ma_star ** 2 - 1)) / 2)))

        return ma

    def area_ratio_from_ma(self, ma):

        ma_k = self.critical_mach_from_mach(ma)
        k = self.k
        area_ratio = ma_k * (((1 - ((k - 1) / (k + 1)) * ma_k ** 2) /
                                (1 - ((k - 1) / (k + 1)))) ** (1 / (k - 1)))

        return area_ratio



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