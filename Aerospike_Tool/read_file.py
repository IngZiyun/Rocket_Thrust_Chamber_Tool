import numpy as np
import re

class Read_RPA():

    def __init__(self, address):
        '''get the address of RPA-file and read the parameters from it.
           '''
        self.add = address # get the address of file.

        name_product, self.Chamber_product, self.Throat_product =             \
        np.loadtxt(self.add, dtype=str, comments='#', skiprows=50,
                   delimiter="\t", usecols=(0, 2, 4), encoding="Latin-1",
                   unpack=True)
        # get the columns for names of products of propellants (0th column)
        # and the molar fraction in chamber (2th column),
        # at throat (4th column).

        name, self.Chamber, self.Throat, self.Exit =                          \
        np.loadtxt(self.add, dtype=str, comments='#', delimiter="\t",
                   usecols=(0, 2, 3, 4), encoding="Latin-1", unpack=True)
        # get the columns for names thermo parameters  (0th column)
        # and their values in chamber (2th column),
        # at throat (3th column), at exit (4th column).

        self.name = name.tolist()
        self.name_product = name_product.tolist()
        # translate the name and name_product to list-value.

    def h2o(self):
        '''get the index of desired parameter's name, the desired values of
           the parameter at different position can also be got with
           the same index. And translate the desired values to float.
           '''
        h2o = self.name_product.index('                 H2O')
        h2o_value_cham = float(self.Chamber_product[h2o])
        h2o_value_thro = float(self.Throat_product[h2o])
        return h2o_value_cham, h2o_value_thro

    def co2(self):
        '''get the index of desired parameter's name, the desired values of
                   the parameter at different position can also be got with
                   the same index. And translate the desired values to float.
                   '''
        co2 = self.name_product.index('                 CO2')
        co2_value_cham = float(self.Chamber_product[co2])
        co2_value_thro = float(self.Throat_product[co2])
        return co2_value_cham, co2_value_thro

    def pressure(self):
        '''get the index of desired parameter's name, the desired values of
                   the parameter at different position can also be got with
                   the same index. And translate the desired values to float.
                   '''
        pressure = self.name.index('            Pressure')
        pressure_value_cham = float(self.Chamber[pressure]) * 1000000
        pressure_value_thro = float(self.Throat[pressure]) * 1000000
        pressure_value_exit = float(self.Exit[pressure]) * 1000000
        # '*1000000' means: translate MPa to Pa.
        return pressure_value_cham, pressure_value_thro, pressure_value_exit

    def temperature(self):
        '''get the index of desired parameter's name, the desired values of
                   the parameter at different position can also be got with
                   the same index. And translate the desired values to float.
                   '''
        temperature = self.name.index('         Temperature')
        temperature_value_cham = float(self.Chamber[temperature])
        temperature_value_thro = float(self.Throat[temperature])
        temperature_value_exit = float(self.Exit[temperature])
        return temperature_value_cham, temperature_value_thro,                \
               temperature_value_exit

    def cp(self):
        '''get the index of desired parameter's name, the desired values of
                   the parameter at different position can also be got with
                   the same index. And translate the desired values to float.
                   '''
        cp = self.name.index('Specific heat (p=const)')
        cp_value_cham = float(self.Chamber[cp]) * 1000
        cp_value_thro = float(self.Throat[cp]) * 1000
        cp_value_exit = float(self.Exit[cp]) * 1000
        # '*1000' means: translate kj to j.
        return cp_value_cham, cp_value_thro, cp_value_exit

    def gamma(self):
        '''get the index of desired parameter's name, the desired values of
                   the parameter at different position can also be got with
                   the same index. And translate the desired values to float.
                   '''
        gamma = self.name.index('               Gamma')
        gamma_value_cham = float(self.Chamber[gamma])
        gamma_value_thro = float(self.Throat[gamma])
        gamma_value_exit = float(self.Exit[gamma])
        return gamma_value_cham, gamma_value_thro, gamma_value_exit

    def gas_constant(self):
        '''get the index of desired parameter's name, the desired values of
                   the parameter at different position can also be got with
                   the same index. And translate the desired values to float.
                   '''
        R = self.name.index('        Gas constant')
        R_value_cham = float(self.Chamber[R]) * 1000
        R_value_thro = float(self.Throat[R]) * 1000
        R_value_exit = float(self.Exit[R]) * 1000
        # '*1000' means: translate kj to j.
        return R_value_cham, R_value_thro, R_value_exit

    def prandtl_number(self):
        '''get the index of desired parameter's name, the desired values of
                   the parameter at different position can also be got with
                   the same index. And translate the desired values to float.
                   '''
        Pr = self.name.index('Prandtl number, effective')
        Pr_value_cham = float(self.Chamber[Pr])
        Pr_value_thro = float(self.Throat[Pr])
        Pr_value_exit = float(self.Exit[Pr])
        return Pr_value_cham, Pr_value_thro, Pr_value_exit

    def viscosity(self):
        '''get the index of desired parameter's name, the desired values of
                   the parameter at different position can also be got with
                   the same index. And translate the desired values to float.
                   '''
        vis = self.name.index('           Viscosity')
        vis_value_cham = float(self.Chamber[vis])
        vis_value_thro = float(self.Throat[vis])
        vis_value_exit = float(self.Exit[vis])
        return vis_value_cham, vis_value_thro, vis_value_exit

    def density(self):
        '''get the index of desired parameter's name, the desired values of
                   the parameter at different position can also be got with
                   the same index. And translate the desired values to float.
                   '''
        rho = self.name.index('             Density')
        rho_value_cham = float(self.Chamber[rho])
        rho_value_thro = float(self.Throat[rho])
        rho_value_exit = float(self.Exit[rho])
        return rho_value_cham, rho_value_thro, rho_value_exit

    def sonic_velocity(self):
        '''get the index of desired parameter's name, the desired values of
                   the parameter at different position can also be got with
                   the same index. And translate the desired values to float.
                   '''
        a = self.name.index('      Sonic velocity')
        a_value_cham = float(self.Chamber[a])
        a_value_thro = float(self.Throat[a])
        a_value_exit = float(self.Exit[a])
        return a_value_cham, a_value_thro, a_value_exit

    def mach_number(self):
        '''get the index of desired parameter's name, the desired values of
                   the parameter at different position can also be got with
                   the same index. And translate the desired values to float.
                   '''
        ma = self.name.index('         Mach number')
        ma_value_cham = float(self.Chamber[ma])
        ma_value_thro = float(self.Throat[ma])
        ma_value_exit = float(self.Exit[ma])
        return ma_value_cham, ma_value_thro, ma_value_exit


    def complete(self):
        '''get all important parameters from the RPA-file.
           And write them in a dictionary-value.'''

        values = dict()# set a dictionary-value.
        values['pressure chamber'], values['pressure throat'],                \
        values['pressure exit'] = self.pressure()

        values['temperature chamber'], values['temperature throat'],          \
        values['temperature exit'] = self.temperature()

        values['cp chamber'], values['cp throat'],                            \
        values['cp exit'] = self.cp()

        values['gamma chamber'], values['gamma throat'],                      \
        values['gamma exit'] = self.gamma()

        values['gas constant chamber'], values['gas constant throat'],        \
        values['gas constant exit'] = self.gas_constant()

        values['prandtl number chamber'], values['prandtl number throat'],    \
        values['prandtl number exit'] = self.prandtl_number()

        values['viscosity chamber'], values['viscosity throat'],              \
        values['viscosity exit'] = self.viscosity()

        values['mach chamber'], values['mach throat'],                        \
        values['mach exit'] = self.mach_number()

        values['density chamber'], values['density throat'],                  \
        values['density exit'] = self.density()

        values['sonic velocity chamber'], values['sonic velocity throat'],    \
        values['sonic velocity exit'] = self.sonic_velocity()

        values['h2o chamber'], values['h2o throat'] = self.h2o()
        values['co2 chamber'], values['co2 throat'] = self.co2()

        return values



class Read_CEA():

    def __init__(self, address):
        '''get the address of CEA file, and read it per line.
           '''
        self.add = address
        with open(self.add, 'r') as f:
            self.read = f.readlines()

    def pressure(self):
        '''translate every line to utf-8 code.
           get the first n code of every line.
           (n is length of codes of the desired parameter.)
           translate the name of desired parameter to uft-8 code.
           if the n codes and the codes of desired parameter's name
           are the same. Then get the first, second and third number
           of this line as the value of desired parameter in chamber,
           at throat and at exit.
           '''
        for i in self.read:
            a = i.encode('utf-8')
            b = a[0:7]
            c = ' P, BAR'
            c = c.encode('utf-8')
            if b == c:
                value_cham, value_throat, value_exit =                        \
                re.findall('\d+.?\d*', i) # get the numbers in this line.

        return float(value_cham) * 100000, float(value_throat) * 100000,      \
        float(value_exit) * 100000
        # '*100000' means: translate Bar in Pascal.

    def temperature(self):
        '''translate every line to utf-8 code.
                   get the first n code of every line.
                   (n is length of codes of the desired parameter.)
                   translate the name of desired parameter to uft-8 code.
                   if the n codes and the codes of desired parameter's name
                   are the same. Then get the first, second and third number
                   of this line as the value of desired parameter in chamber,
                   at throat and at exit.
                   '''
        for i in self.read:
            a = i.encode('utf-8')
            b = a[0:5]
            c = ' T, K'
            c = c.encode('utf-8')
            if b == c:
                value_cham, value_throat, value_exit =                        \
                re.findall('\d+.?\d*', i)# get the numbers in this line.

        return float(value_cham), float(value_throat), float(value_exit)

    def cp(self):
        '''translate every line to utf-8 code.
                   get the first n code of every line.
                   (n is length of codes of the desired parameter.)
                   translate the name of desired parameter to uft-8 code.
                   if the n codes and the codes of desired parameter's name
                   are the same. Then get the first, second and third number
                   of this line as the value of desired parameter in chamber,
                   at throat and at exit.
                   '''
        for i in self.read:
            a = i.encode('utf-8')
            b = a[0:15]
            c = ' Cp, KJ/(KG)(K)'
            c = c.encode('utf-8')
            if b == c:
                value_cham, value_throat, value_exit =                        \
                re.findall('\d+.?\d*', i)# get the numbers in this line.

        return float(value_cham) * 1000, float(value_throat) * 1000,          \
        float(value_exit) * 1000
        # '*1000' means: translate kj to j.

    def gamma(self):
        '''translate every line to utf-8 code.
                   get the first n code of every line.
                   (n is length of codes of the desired parameter.)
                   translate the name of desired parameter to uft-8 code.
                   if the n codes and the codes of desired parameter's name
                   are the same. Then get the first, second and third number
                   of this line as the value of desired parameter in chamber,
                   at throat and at exit.
                   '''
        for i in self.read:
            a = i.encode('utf-8')
            b = a[0:7]
            c = ' GAMMAs'
            c = c.encode('utf-8')
            if b == c:
                value_cham, value_throat, value_exit =                        \
                re.findall('\d+.?\d*', i)# get the numbers in this line.

        return float(value_cham), float(value_throat), float(value_exit)


    def sonic_velocity(self):
        '''translate every line to utf-8 code.
                   get the first n code of every line.
                   (n is length of codes of the desired parameter.)
                   translate the name of desired parameter to uft-8 code.
                   if the n codes and the codes of desired parameter's name
                   are the same. Then get the first, second and third number
                   of this line as the value of desired parameter in chamber,
                   at throat and at exit.
                   '''
        for i in self.read:
            a = i.encode('utf-8')
            b = a[0:14]
            c = ' SON VEL,M/SEC'
            c = c.encode('utf-8')
            if b == c:
                value_cham, value_throat, value_exit =                        \
                re.findall('\d+.?\d*', i)# get the numbers in this line.

        return float(value_cham), float(value_throat), float(value_exit)

    def gas_constant(self):
        '''gas constant is not given in CEA file, but the gas constant
           can be calculate through sonic velocity, isentropic exponent
           the temperature.
           so get their values in chamber, at throat and at exit
           to calculate the gas contant at the three position.
           '''
        a_c, a_t, a_e = self.sonic_velocity()
        k_c, k_t, k_e = self.gamma()
        t_c, t_t, t_e = self.temperature()
        r_c = (a_c ** 2) / (k_c * t_c)
        r_t = (a_t ** 2) / (k_t * t_t)
        r_e = (a_e ** 2) / (k_e * t_e)
        return r_c, r_t, r_e

    def prandtl_number(self):
        '''prandtl number is not given in CEA file, but it
           can be calculate through isentropic exponent.
           so get isentropic exponent in chamber, at throat and
           at exit to calculate the prandtl number at the three position.
                   '''
        k_c, k_t, k_e = self.gamma()
        pr_c = (4 * k_c) / (9 * k_c - 5)
        pr_t = (4 * k_t) / (9 * k_t - 5)
        pr_e = (4 * k_e) / (9 * k_e - 5)
        return pr_c, pr_t, pr_e

    def density(self):
        '''the value of density is not easy to read in CEA file, but density
           can be calculate through gas constant, pressure and temperature.
           so get their values in chamber, at throat and at exit
           to calculate the gas contant at the three position.
                   '''
        r_c, r_t, r_e = self.gas_constant()
        t_c, t_t, t_e = self.temperature()
        p_c, p_t, p_e = self.pressure()
        rho_c = p_c/(r_c * t_c)
        rho_t = p_t / (r_t * t_t)
        rho_e = p_e / (r_e * t_e)
        return rho_c, rho_t, rho_e

    def mach_number(self):
        '''translate every line to utf-8 code.
                   get the first n code of every line.
                   (n is length of codes of the desired parameter.)
                   translate the name of desired parameter to uft-8 code.
                   if the n codes and the codes of desired parameter's name
                   are the same. Then get the first, second and third number
                   of this line as the value of desired parameter in chamber,
                   at throat and at exit.
                   '''
        for i in self.read:
            a = i.encode('utf-8')
            b = a[0:12]
            c = ' MACH NUMBER'
            c = c.encode('utf-8')
            if b == c:
                value_cham, value_throat, value_exit =                        \
                re.findall('\d+.?\d*', i)# get the numbers in this line.

        return float(value_cham), float(value_throat), float(value_exit)

    def h2o(self):
        '''translate every line to utf-8 code.
                   get the first n code of every line.
                   (n is length of codes of the desired parameter.)
                   translate the name of desired parameter to uft-8 code.
                   if the n codes and the codes of desired parameter's name
                   are the same. Then get the first, second and third number
                   of this line as the value of desired parameter in chamber,
                   at throat and at exit.
                   '''
        for i in self.read:
            a = i.encode('utf-8')
            b = a[0:4]
            c = ' H2O'
            c = c.encode('utf-8')
            if b == c:
                _filter, value_cham, value_throat, value_exit =               \
                re.findall('\d+.?\d*', i)# get the numbers in this line.

        return float(value_cham), float(value_throat), float(value_exit)

    def co2(self):
        '''translate every line to utf-8 code.
                   get the first n code of every line.
                   (n is length of codes of the desired parameter.)
                   translate the name of desired parameter to uft-8 code.
                   if the n codes and the codes of desired parameter's name
                   are the same. Then get the first, second and third number
                   of this line as the value of desired parameter in chamber,
                   at throat and at exit.
                   '''
        for i in self.read:
            a = i.encode('utf-8')
            b = a[0:5]
            c = ' *CO2'
            c = c.encode('utf-8')
            if b == c:
                _filter, value_cham, value_throat, value_exit =               \
                re.findall('\d+.?\d*', i)# get the numbers in this line.

        return float(value_cham), float(value_throat), float(value_exit)

    def viscosity(self):
        '''the values of viscosity are not given in CEA file, but they can
           be calculated through specific gas constant, temperature and
           the molar masses. The molar masses can be taken through
           specific constant and universelle gas constant.
           '''
        r_c, r_t, r_e = self.gas_constant()
        t_c, t_t, t_e = self.temperature()
        R = 8.314462 # universelle gas constant.
        value_cham = (46.6 * (10 ** (-10))) * ((R/r_c) ** 0.5) * (t_c ** 0.6)
        value_throat = (46.6 * (10 ** (-10))) * ((R/r_t) ** 0.5) * (t_t ** 0.6)
        value_exit = (46.6 * (10 ** (-10))) * ((R/r_e) ** 0.5) * (t_e ** 0.6)
        return value_cham, value_throat, value_exit

    def complete(self):
        '''get all important parameters from the CEA-file.
                   And write them in a dictionary-value.'''

        values = dict()

        values['pressure chamber'], values['pressure throat'],                \
        values['pressure exit'] = self.pressure()

        values['temperature chamber'], values['temperature throat'],          \
        values['temperature exit'] = self.temperature()

        values['cp chamber'], values['cp throat'],                            \
        values['cp exit'] = self.cp()

        values['gamma chamber'], values['gamma throat'],                      \
        values['gamma exit'] = self.gamma()

        values['gas constant chamber'], values['gas constant throat'],        \
        values['gas constant exit'] = self.gas_constant()

        values['prandtl number chamber'], values['prandtl number throat'],    \
        values['prandtl number exit'] = self.prandtl_number()

        values['viscosity chamber'], values['viscosity throat'],              \
        values['viscosity exit'] = self.viscosity()

        values['mach chamber'], values['mach throat'],                        \
        values['mach exit'] = self.mach_number()

        values['h2o chamber'], values['h2o throat'],                          \
        values['h2o exit'] = self.h2o()

        values['co2 chamber'], values['co2 throat'],                          \
        values['co2 exit'] = self.co2()

        values['density chamber'], values['density throat'],                  \
        values['density exit'] = self.density()
        values['sonic velocity chamber'], values['sonic velocity throat'],    \
        values['sonic velocity exit'] = self.sonic_velocity()

        return values



