from PySide6.QtWidgets import QApplication, QMessageBox, QLabel
from PySide6.QtUiTools import QUiLoader
from PySide6 import QtWidgets
import sys
sys.path.append('./UI')
import Aerospike_P
import json
import Json
import Aerospike
import read_file
import numpy as np
import matplotlib.pyplot as plt

class User_interface:

    def __init__(self):
        '''load the interface-date of QT-Designer, and set the connections
           for Buttons.'''

        self.ui = QUiLoader().load('UI/EASE UI.ui')
        #load interface-data of QT-Designer.

        self.ui.file_load.clicked.connect(self.load_file)
        self.ui.select_ad.clicked.connect(self.select_file)
        self.ui.con_check.stateChanged.connect(self.check_contour)
        self.ui.select_in.clicked.connect(self.select_in)
        self.ui.select_out.clicked.connect(self.select_out)
        self.ui.linear.clicked.connect(self.lin_anu)
        self.ui.annular.clicked.connect(self.lin_anu)
        self.ui.lin2.clicked.connect(self.lin_anu_2)
        self.ui.anu2.clicked.connect(self.lin_anu_2)
        self.ui.input_contour.clicked.connect(self.contour_input)
        self.ui.heat_cal.clicked.connect(self.calculate_heat)
        self.ui.show_contour.clicked.connect(self.contour_show)
        self.ui.heat_analyse.clicked.connect(self.heatflux_show)
        self.ui.contour_cal.clicked.connect(self.calculate_contour)
        self.ui.contour_output.clicked.connect(self.output_contour)
        self.ui.heat_output.clicked.connect(self.output_heat)
        self.ui.area_analyse.clicked.connect(self.analyse_area)
        self.ui.chamber_analyse.clicked.connect(self.analyse_chamber)
        # set connections for Buttons.

        self.ui.setWindowTitle('Engineering AeroSpike Enignes Tool')
        # Window title.

        self.ui.linear.setChecked(True)
        # set a Button 'linear' as True.



    def select_file(self):
        '''select the RPA- or CEA-File.
           '''
        file, _filter = QtWidgets.QFileDialog.getOpenFileName(None,
                                        "Select", "./", "Text Files(*.txt)")
        # open the window for selection of files, and return the file address.

        self.file_ad = str(file)
        # get the file address, and translate it to string.

        self.ui.file_ad.setText(self.file_ad)
        # set the string of file address in the text field.



    def select_in(self):
        '''select the File for inner contour.
                   '''
        con_in, _filter = QtWidgets.QFileDialog.getOpenFileName(None,
                                        "Select", "./", "Text Filex(*txt)")
        # open the window for selection of files, and return the file address.

        self.con_in_ad = str(con_in)
        # get the file address, and translate it to string.

        self.ui.ad_in.setText(self.con_in_ad)
        # set the string of file address in the text field.



    def select_out(self):
        '''select the File for outer contour.
                           '''
        con_out, _filter = QtWidgets.QFileDialog.getOpenFileName(None,
                                        "Select", "./", "Text Filex(*txt)")
        # open the window for selection of files, and return the file address.

        self.con_out_ad = str(con_out)
        # get the file address, and translate it to string.

        self.ui.ad_out.setText(self.con_out_ad)
        # set the string of file address in the text field.



    def check_contour(self):
        '''check it out, the user want input external contour-data or not.
           And set situation of elements for the choice.'''
        if self.ui.con_check.isChecked() ==1:
            # if the user want input external contour-data.
            self.ui.ad_in.setEnabled(True)
            self.ui.ad_out.setEnabled(True)
            self.ui.select_in.setEnabled(True)
            self.ui.select_out.setEnabled(True)
            self.ui.lin2.setEnabled(True)
            self.ui.anu2.setEnabled(True)
            self.ui.depth2.setEnabled(True)
            self.ui.meter.setEnabled(True)
            self.ui.input_contour.setEnabled(True)
            self.ui.linear.setChecked(True)
            self.ui.lin2.setChecked(True)

            self.ui.annular.setEnabled(False)
            self.ui.linear.setEnabled(False)
            self.ui.dep.setEnabled(False)
            self.ui.contour_number.setEnabled(False)
            self.ui.y_exp.setEnabled(False)
            self.ui.shortening.setEnabled(False)
            self.ui.conv_ang.setEnabled(False)
            self.ui.in_r.setEnabled(False)
            self.ui.out_r.setEnabled(False)
            self.ui.l_cham.setEnabled(False)
            self.ui.l_conv.setEnabled(False)
            self.ui.distance.setEnabled(False)
            self.ui.contour_cal.setEnabled(False)
            # set the situation of all elements.

        else:
            # if the user don't want external contour-data.
            self.ui.ad_in.setEnabled(False)
            self.ui.ad_out.setEnabled(False)
            self.ui.select_in.setEnabled(False)
            self.ui.select_out.setEnabled(False)
            self.ui.lin2.setEnabled(False)
            self.ui.lin2.setChecked(True)
            self.ui.anu2.setEnabled(False)
            self.ui.depth2.setEnabled(False)
            self.ui.meter.setEnabled(False)
            self.ui.input_contour.setEnabled(False)

            self.ui.annular.setEnabled(True)
            self.ui.linear.setEnabled(True)
            self.ui.dep.setEnabled(True)
            self.ui.contour_number.setEnabled(True)
            self.ui.y_exp.setEnabled(True)
            self.ui.shortening.setEnabled(True)
            self.ui.conv_ang.setEnabled(True)
            self.ui.in_r.setEnabled(True)
            self.ui.out_r.setEnabled(True)
            self.ui.l_cham.setEnabled(True)
            self.ui.l_conv.setEnabled(True)
            self.ui.distance.setEnabled(True)
            self.ui.contour_cal.setEnabled(True)
            self.ui.linear.setChecked(True)
            # set the situation of all elements.



    def lin_anu(self):
        '''check it out, the user want linear or annular aerospike.
                   And set situation of elements for the choice.
                   '''
        if self.ui.linear.isChecked() == 1: # if linear aerospike.
            self.ui.dep.setEnabled(True)
            self.ui.dep.setText('100')
        else:# if annular aerospike.
            self.ui.dep.setEnabled(False)
            self.ui.dep.setText('0')


    def lin_anu_2(self):
        '''check it out, the user want linear or annular aerospike.
                           And set situation of elements for the choice.
                           '''
        if self.ui.lin2.isChecked() == 1: # if linear aerospike.
            self.ui.depth2.setEnabled(True)
            self.ui.depth2.setText('100')
        else:# if annular aerospike.
            self.ui.depth2.setEnabled(False)
            self.ui.depth2.setText('0')



    def load_file(self):
        '''load the wanted parameters from RPA- or CEA-File.
           '''
        if self.ui.rpa.isChecked() == 1:# if RPA.
            read = read_file.Read_RPA(self.ui.file_ad.text())
            # use the class for reading the data from RPA
            data = read.complete()
            # use the function for reading all parameters, and save them
            # in value 'data'.
            json_ad, _filter = QtWidgets.QFileDialog.getSaveFileName(None,
                                        "Save", "./", "Json Files(*.json)")
            # open the window for creation of a JSON-File, and get the
            # address the of JSON-File.
            self.json_data = Json.Json_Edit(json_ad)
            # use the class for Edition of JSON-File.
            self.json_data.write(data)
            # write the data in JSON-File.

        else:# if CEA.
            read = read_file.Read_CEA(self.ui.file_ad.text())
            # use the class for reading the data from CEA.
            data = read.complete()
            # use the function for reading all parameters, and save them
            # in value 'data'.
            json_ad, _filter = QtWidgets.QFileDialog.getSaveFileName(None,
                                        "Save", "./", "Json Files(*.json)")
            # open the window for creation of a JSON-File, and get the
            # address the of JSON-File.
            self.json_data = Json.Json_Edit(json_ad)
            # use the class for Edition of JSON-File.
            self.json_data.write(data)
            # write the data in JSON-File.

        msg = QMessageBox.information(None,
                                      "Messagebox", "File load is finished")
        print(msg)
        # give a Message.



    def calculate_contour(self):
        '''generate the contour-data with wanted geometric parameters.
           '''
        depth = float(self.ui.dep.text()) / 1000
        contour_N = float(self.ui.contour_number.text())
        y_exp = float(self.ui.y_exp.text()) / 1000
        shortening = float(self.ui.shortening.text())
        inner_radius = float(self.ui.in_r.text()) / 1000
        outer_radius = float(self.ui.out_r.text()) / 1000
        l_chamber = float(self.ui.l_cham.text()) / 1000
        l_convegence = float(self.ui.l_conv.text()) / 1000
        delta_s = float(self.ui.distance.text()) / 1000
        gas_parameter = self.json_data.read(0)
        # get the necessary parameters for the calculation.

        if self.ui.linear.isChecked() == 1: # if linear aerospike is desired.
            form = 0 # set the form value as 0
        else: # if annular aerospike is desired.
            form = 1 # set the form value as 1

        contour = Aerospike.Aerospike_Geometry(form)
        # use the class for aerospike geometry with form value.

        x_spike_r, y_spike_r, convergence_angle =                             \
        contour.spike_contour(gas_parameter['gamma throat'], contour_N,
                              gas_parameter['mach exit'], shortening)
        # get the relationship between the coordinates of points on
        # aerospike-contour and the Y-coordinate of expansion point.
        # get the convergence angle(total turning angle).

        convergence_angle_degree = 180 * (convergence_angle/np.pi)
        # translate the angle from degree to radial.

        x_out, y_out = contour.outer_chamber_contour(0, y_exp,
            convergence_angle_degree, outer_radius, l_chamber, l_convegence,
                                                     delta_s)
        # get the coordinates of the points on outer contour of
        # combustion chamber.

        x_end = min(x_out)
        # get the x-coordinate of first point on the parallel part
        # of combustion chamber.

        x_spike = [ix * y_exp for ix in x_spike_r]
        y_spike = [iy * y_exp for iy in y_spike_r]
        # get the coordinates of the points on aerospike-contour.

        faktor =                                                              \
        ((inner_radius * convergence_angle / (2 * np.pi)) / l_convegence)
        # get the factor for refinement of the inner part.

        x_cham_in, y_cham_in =                                                \
        contour.inner_chamber_contour(x_spike[0], y_spike[0],
                convergence_angle_degree, inner_radius, x_end, delta_s*faktor)
        # get the coordinates of the points on inner contour of
        # combustion chamber.

        x_in = x_cham_in + x_spike
        y_in = y_cham_in + y_spike
        # assemble the all points on inner contour.
        # (inner contour of combustion chamber + aerospike-contour.)

        contour_data = dict()
        # set a dictionary value for saving the coordinates.
        contour_data['x in'] = x_in
        contour_data['y in'] = y_in
        contour_data['x out'] = x_out
        contour_data['y out'] = y_out
        contour_data['depth'] = depth
        # set the coordinates in corresponding groups.

        n = self.json_data.check_line_number()  # check json file.
        if n > 1:  # if the second line is not blank, then use json.revise
            self.json_data.revise(1, contour_data)
        else:  # if the second line is blank, then use json.write
            self.json_data.write(contour_data)

        self.ui.conv_ang.setText(str(round(convergence_angle_degree, 2)))
        # writ the dictionary 'contour_data' in the second line
        # of the JSON-File.
        msg =                                                                 \
        QMessageBox.information(None, "Messagebox",
                                "Calculate of contour is finished")
        print(msg)
        # give a Message.



    def contour_input(self):
        '''read the geometric data, if users want input external contour.
           '''
        address_in = self.ui.ad_in.text()
        address_out = self.ui.ad_out.text()
        # get the addresses of contour-data.
        x_in, y_in =                                                          \
        np.loadtxt(address_in, dtype=float, usecols=(0, 1), unpack=True)
        # read the coordinates of points on inner contour.
        x_out, y_out =                                                        \
        np.loadtxt(address_out, dtype=float, usecols=(0, 1), unpack=True)
        # read the coordinates of points on outer contour.
        x_in = x_in.tolist()
        y_in = y_in.tolist()
        x_out = x_out.tolist()
        y_out = y_out.tolist()
        # translate the coordinates-data to list value.
        if self.ui.meter.isChecked() == 0:
            # if unit is not meter, then translate to meter
            n = len(x_in)
            m = len(x_out)
            for i in range(n):
                x_in[i] = x_in[i] / 1000
                y_in[i] = y_in[i] / 1000
            for im in range(m):
                x_out[im] = x_out[im] / 1000
                y_out[im] = y_out[im] / 1000
        contour = dict()
        contour['x in'] = x_in
        contour['y in'] = y_in
        contour['x out'] = x_out
        contour['y out'] = y_out
        contour['depth'] = float(self.ui.depth2.text()) / 1000
        # set a dictionary value for saving the coordinates, and set
        # the coordinates in corresponding groups.

        n = self.json_data.check_line_number()  # check json file.
        if n > 1:  # if the second line is not blank, then use 'json.revise'.
            self.json_data.revise(1, contour)
        else:  # if the second line is blank, then use 'json.write'.
            self.json_data.write(contour)

        msg = QMessageBox.information(None,
                                      "Messagebox", "File load is finished")
        print(msg)
        # give a Message.



    def calculate_heat(self):
        '''calculate the thermo load.
           '''
        Twall = float(self.ui.wall_temp.text())
        # get the desired wall temperature.
        gas_parameter = self.json_data.read(0)
        # read the gasdynamic parameters.
        contour_data = self.json_data.read(1)
        # read the geometric data.
        lin_anu = self.ui.lin2.isChecked() + self.ui.linear.isChecked()
        #check the form, if one of the two checkbox for linear is checked,
        # then use linear methode

        Tem = list()
        cp = list()
        Tem.append(gas_parameter['temperature chamber'])
        Tem.append(gas_parameter['temperature throat'])
        Tem.append(gas_parameter['temperature exit'])
        # get 3 values of temperature.
        cp.append(gas_parameter['cp chamber'])
        cp.append(gas_parameter['cp throat'])
        cp.append(gas_parameter['cp exit'])
        # get 3 values of cp.
        cp_koe = np.polyfit(Tem, cp, 2)
        # get the 2 degree polynomial of cp about Temperature.

        if lin_anu > 0:
            # at last one of the two checkbox for linear is checked.
            # so use linear methode.
            Geo = Aerospike.Aerospike_Geometry(0)
            # use class for aerospike geometry with linear form.

            x_cham_in, y_cham_in, x_throat, y_throat, x_spike, y_spike,       \
            x_exp, y_exp, x_cham_out, y_cham_out, x_r_throat, y_r_throat,     \
            r_throat = Geo.points_arrange(
                contour_data['x in'], contour_data['y in'],
                contour_data['x out'], contour_data['y out'])
            # enter the geometric data,
            # and get the coordinates of points on
            # inner contour (x_cham_in, y_cham_in) of combustion chamber,
            # outer contour (x_cham_out, y_cham_out) of combustion chamber,
            # points on aerospike (x_spike, y_spike),
            # throat point (x_throat, y_throat),
            # expansion point (x_exp, y_exp),
            # center of rounding circle arc at throat (x_r_throat, y_r_throat)
            # radius of rounding circle arc at throat (r_throat)

            thermo_load = Aerospike.Thermo_load(
                gas_parameter['pressure chamber'],
                gas_parameter['temperature chamber'], Twall)
            # use class for calculation of thermo load with parameters :
            # chamber pressure, chamber temperature, desired wall temperature.

            x_spike.insert(0, x_throat)
            y_spike.insert(0, y_throat)
            # insert the coordinate of throat point as the first point
            # of aerospike.

            chamber_conv, chamber_conv_out, chamber_radi, chamber_flux,       \
            chamber_flux_out, chamber_in_load, chamber_out_load,                                \
            chamber_heat_transfer_coefficient,                                \
            chamber_heat_transfer_coefficient_out =                           \
            thermo_load.lin_chamber(
                x_cham_in, y_cham_in, x_cham_out, y_cham_out,
                x_throat, y_throat, x_exp, y_exp, contour_data['depth'],
                r_throat, gas_parameter['gamma chamber'],
                gas_parameter['viscosity chamber'],
                cp_koe[0], cp_koe[1], cp_koe[2],
                gas_parameter['gas constant chamber'],
                gas_parameter['prandtl number chamber'],
                gas_parameter['h2o chamber'], gas_parameter['co2 chamber'])
            # enter the geometric parameters and gasdynamic parameters,
            # get the convective heat flux density in chamber(chamber_conv),
            # radiation heat flux density in chamber(chamber_radi),
            # total heat flux density in chamber(chamber_flux),
            # total heat transfer power on inner part of chamber
            # (chamber_in_load),
            # total heat transfer power on outer part of chamber
            # (chamber_out_load),
            # convective heat transfer coefficient in chamber
            # (chamber_heat_transfer_coefficient).

            spike_conv, spike_radi, spike_flux, spike_load,                   \
            spike_heat_transfer_coefficient =                                 \
            thermo_load.lin_spike(x_spike, y_spike, x_exp, y_exp,
                                  contour_data['depth'], r_throat,
                                  gas_parameter['gamma throat'],
                                  gas_parameter['viscosity throat'],
                                  cp_koe[0], cp_koe[1], cp_koe[2],
                                  gas_parameter['gas constant throat'],
                                  gas_parameter['prandtl number throat'],
                                  gas_parameter['h2o throat'],
                                  gas_parameter['co2 throat'])
            # enter the geometric parameters and gasdynamic parameters,
            # get the convective heat flux density on aerospike(spike_conv),
            # radiation heat flux density on aerospike(spike_radi),
            # total heat flux density on aerospike(spike_flux),
            # total heat transfer power on aerospike
            # (spike_load),
            # convective heat transfer coefficient on aerospike.
            # (spike_heat_transfer_coefficient).

        else: # use annular methode
            Geo = Aerospike.Aerospike_Geometry(1)
            # use class for aerospike geometry with annular form.
            x_cham_in, y_cham_in, x_throat, y_throat, x_spike, y_spike,       \
            x_exp, y_exp, x_cham_out, y_cham_out, x_r_throat, y_r_throat,     \
            r_throat = Geo.points_arrange(
                contour_data['x in'], contour_data['y in'],
                contour_data['x out'], contour_data['y out'])
            # enter the geometric data,
            # and get the coordinates of points on
            # inner contour (x_cham_in, y_cham_in) of combustion chamber,
            # outer contour (x_cham_out, y_cham_out) of combustion chamber,
            # points on aerospike (x_spike, y_spike),
            # throat point (x_throat, y_throat),
            # expansion point (x_exp, y_exp),
            # center of rounding circle arc at throat (x_r_throat, y_r_throat)
            # radius of rounding circle arc at throat (r_throat)

            x_spike.insert(0, x_throat)
            y_spike.insert(0, y_throat)
            # insert the coordinate of throat point as the first point
            # of aerospike.

            thermo_load =                                                     \
            Aerospike.Thermo_load(gas_parameter['pressure chamber'],
                                  gas_parameter['temperature chamber'], Twall)
            # use class for calculation of thermo load with parameters:
            # chamber pressure, chamber temperature, desired wall temperature.

            chamber_conv, chamber_conv_out, chamber_radi, chamber_flux,       \
            chamber_flux_out, chamber_in_load, chamber_out_load,              \
            chamber_heat_transfer_coefficient,                                \
            chamber_heat_transfer_coefficient_out =                           \
            thermo_load.anu_chamber(
                x_cham_in, y_cham_in, x_cham_out, y_cham_out,
                x_throat, y_throat, x_exp, y_exp, r_throat,
                gas_parameter['gamma chamber'],
                gas_parameter['viscosity chamber'],
                cp_koe[0], cp_koe[1], cp_koe[2],
                gas_parameter['gas constant chamber'],
                gas_parameter['prandtl number chamber'],
                gas_parameter['h2o chamber'], gas_parameter['co2 chamber'])
            # enter the geometric parameters and gasdynamic parameters,
            # get the convective heat flux density in chamber(chamber_conv),
            # radiation heat flux density in chamber(chamber_radi),
            # total heat flux density in chamber(chamber_flux),
            # total heat transfer power on inner part of chamber
            # (chamber_in_load),
            # total heat transfer power on outer part of chamber
            # (chamber_out_load),
            # convective heat transfer coefficient in chamber
            # (chamber_heat_transfer_coefficient).

            spike_conv, spike_radi, spike_flux, spike_load,                   \
            spike_heat_transfer_coefficient =                                 \
            thermo_load.anu_spike(x_spike, y_spike, x_exp, y_exp, r_throat,
                                  gas_parameter['gamma throat'],
                                  gas_parameter['viscosity throat'],
                                  cp_koe[0], cp_koe[1], cp_koe[2],
                                  gas_parameter['gas constant throat'],
                                  gas_parameter['prandtl number throat'],
                                  gas_parameter['h2o throat'],
                                  gas_parameter['co2 throat'])
            # enter the geometric parameters and gasdynamic parameters,
            # get the convective heat flux density on aerospike(spike_conv),
            # radiation heat flux density on aerospike(spike_radi),
            # total heat flux density on aerospike(spike_flux),
            # total heat transfer power on aerospike
            # (spike_load),
            # convective heat transfer coefficient on aerospike.
            # (spike_heat_transfer_coefficient).

        thermo_flux_engine = chamber_flux + spike_flux
        thermo_load_engine = (chamber_out_load + chamber_in_load + spike_load)
        # assemble the thermo values.
        x_points = x_cham_in + x_spike
        # assemble the x-coordinates.
        thermo_load_values = dict()
        # set a dictionary for saving the thermo data.
        thermo_load_values['Convective Thermo flux density'] =                \
        chamber_conv + spike_conv
        thermo_load_values['Radiation Thermo flux density'] =                 \
        chamber_radi + spike_radi
        thermo_load_values['Total Thermo flux density'] =                     \
        thermo_flux_engine
        thermo_load_values['Inner x coordinate'] = x_points
        thermo_load_values['Outer x coordinate'] = x_cham_out
        thermo_load_values['Convective Thermo flux density in Chamber'] =     \
        chamber_conv
        thermo_load_values['Convective Thermo flux density on Chamber Shell']=\
        chamber_conv_out
        thermo_load_values['Radiation Thermo flux density in Chamber'] =      \
        chamber_radi
        thermo_load_values['Total Thermo flux density on Chamber Shell'] =    \
        chamber_flux_out
        thermo_load_values['Convective Heat transfer coefficient'] =          \
        chamber_heat_transfer_coefficient + spike_heat_transfer_coefficient
        thermo_load_values['Convective Heat transfer coefficient in Chamber'] \
        = chamber_heat_transfer_coefficient
        thermo_load_values[
            'Convective Heat transfer coefficient on Chamber Shell']          \
            = chamber_heat_transfer_coefficient_out
        thermo_load_values['Total Thermo Load in kW'] =                       \
        thermo_load_engine / 1000
        thermo_load_values['Thermo Load on inner core in kW'] =               \
        (spike_load + chamber_in_load) / 1000
        thermo_load_values['Thermo Load on outer shell in kW'] =              \
        chamber_out_load / 1000
        # set the thermo values in corresponding positions.

        n = self.json_data.check_line_number()
        if n > 2:  # if the third line is not blank, then use json.revise
            self.json_data.revise(2, thermo_load_values)
        else:  # if the third line is blank, then use json.write
            self.json_data.write(thermo_load_values)
        msg = QMessageBox.information(None,
                    "Messagebox", "Calculation for Heat Transfer is finished")
        print(msg)
        #give a Message.



    def contour_show(self):
        '''draw the diagram for thrust chamber contour.
           '''

        contour = self.json_data.read(1)
        n_in = len(contour['x in'])
        n_out = len(contour['x out'])
        x_in = list()
        y_in = list()
        x_out = list()
        y_out = list()
        for i in range(n_in):
            x_in.append(contour['x in'][i] * 1000)
            y_in.append(contour['y in'][i] *1000)
        for i in range(n_out):
            x_out.append(contour['x out'][i] * 1000)
            y_out.append(contour['y out'][i] * 1000)

        plt.plot(x_in, y_in, 'k')
        plt.plot(x_out, y_out, 'k')
        plt.xlabel('X [Milimeter]')
        plt.ylabel('Y [Milimeter]')
        plt.show()


    def analyse_area(self):
        '''calculate the relationship of cross-sectional areas,
           and draw diagrams.
           '''
        contour_data = self.json_data.read(1) # read geometric data.
        if contour_data['depth'] == 0: # depth = 0 means annular aerospike
            Geo = Aerospike.Aerospike_Geometry(1)
            # use the class 'Aerospike_Geometry' for annular aerospike.
        else:
            Geo = Aerospike.Aerospike_Geometry(0)
            # use the class 'Aerospike_Geometry' for annular aerospike.

        x_cham_in, y_cham_in, x_throat, y_throat, x_spike, y_spike,           \
        x_exp, y_exp, x_cham_out, y_cham_out, x_r_throat,                     \
        y_r_throat, r_throat                                                  \
        = Geo.points_arrange(
            contour_data['x in'], contour_data['y in'], contour_data['x out'],
            contour_data['y out'])
        # enter the geometric data,
        # and get the coordinates of points on
        # inner contour (x_cham_in, y_cham_in) of combustion chamber,
        # outer contour (x_cham_out, y_cham_out) of combustion chamber,
        # points on aerospike (x_spike, y_spike),
        # throat point (x_throat, y_throat),
        # expansion point (x_exp, y_exp),
        # center of rounding circle arc at throat (x_r_throat, y_r_throat)
        # radius of rounding circle arc at throat (r_throat)

        x_point_in = x_cham_in + x_spike # assemble x-coordinates.
        y_point_in = y_cham_in + y_spike
        #x_spike.insert(0, x_throat)
        #y_spike.insert(0, y_throat)
        area_ratio = Geo.full_area_ratio(contour_data['x in'],
        contour_data['y in'], contour_data['x out'], contour_data['y out'])
        # entry geometric data, get the relationship of all
        # cross-sectional areas.

        m = len(x_cham_out)
        o = len(x_spike)
        fig1, ax1 = plt.subplots(num="Cross Section", figsize=(8, 5))
        ax1.set_xlabel('X [Meter]', fontsize=30)
        ax1.set_ylabel('Y [Meter]', fontsize=30)
        ax1.plot(contour_data['x in'], contour_data['y in'], 'k')
        ax1.plot(contour_data['x out'], contour_data['y out'], 'k')
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        for i in range(m):
            ax1.plot([x_cham_out[i], x_cham_in[i]],
                     [y_cham_out[i], y_cham_in[i]], 'r', linewidth=1)
        for i in range(o):
            ax1.plot([x_exp, x_spike[i]],
                     [y_exp, y_spike[i]], 'r', linewidth=1,)
        ax2 = ax1.twinx()
        ax2.plot(x_point_in, area_ratio, 'b', linewidth=2)
        ax2.set_xlabel('X [Meter]', fontsize=35)
        ax2.set_ylabel('Area Ratio', color='b', fontsize=35)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.show()
        # draw the diagrams for contour and area-relationship.



    def analyse_chamber(self):
        gas_parameter = self.json_data.read(0)
        contour_data = self.json_data.read(1)
        if contour_data['depth'] == 0:  # depth = 0 means annular aerospike
            Geo = Aerospike.Aerospike_Geometry(1)
            # use the class 'Aerospike_Geometry' for annular aerospike.
            x_cham_in, y_cham_in, x_throat, y_throat, x_spike, y_spike,       \
            x_exp, y_exp, x_cham_out, y_cham_out, x_r_throat,                 \
            y_r_throat, r_throat = Geo.points_arrange(
                contour_data['x in'], contour_data['y in'],
                contour_data['x out'], contour_data['y out'])
            # enter the geometric data,
            # and get the coordinates of points on
            # inner contour (x_cham_in, y_cham_in) of combustion chamber,
            # outer contour (x_cham_out, y_cham_out) of combustion chamber,
            # points on aerospike (x_spike, y_spike),
            # throat point (x_throat, y_throat),
            # expansion point (x_exp, y_exp),
            # center of rounding circle arc at throat (x_r_throat, y_r_throat)
            # radius of rounding circle arc at throat (r_throat)

            x_cham_in.append(x_throat)
            y_cham_in.append(y_throat)
            x_cham_out.append(x_exp)
            y_cham_out.append(y_exp)
            # set the throat point and expansion point
            # at the end of points group for combustion chamber.

            chamber_volume =                                                  \
            Geo.anu_chamber_volume(x_cham_in, y_cham_in,
                                   x_cham_out, y_cham_out)
            # input the points on combustion chamber, get the volume of
            # combustion chamber.

            a_throat = Geo.supfic_area_cone(x_throat, x_exp, y_throat, y_exp)
            # get the A* (cross-sectional area at throat.)
        else:
            Geo = Aerospike.Aerospike_Geometry(0)
            # use the class 'Aerospike_Geometry' for linear aerospike.

            x_cham_in, y_cham_in, x_throat, y_throat, x_spike, y_spike,       \
            x_exp, y_exp, x_cham_out, y_cham_out, x_r_throat, y_r_throat,     \
            r_throat = Geo.points_arrange(
                contour_data['x in'], contour_data['y in'],
                contour_data['x out'], contour_data['y out'])
            # enter the geometric data,
            # and get the coordinates of points on
            # inner contour (x_cham_in, y_cham_in) of combustion chamber,
            # outer contour (x_cham_out, y_cham_out) of combustion chamber,
            # points on aerospike (x_spike, y_spike),
            # throat point (x_throat, y_throat),
            # expansion point (x_exp, y_exp),
            # center of rounding circle arc at throat (x_r_throat, y_r_throat)
            # radius of rounding circle arc at throat (r_throat)

            x_cham_in.append(x_throat)
            y_cham_in.append(y_throat)
            x_cham_out.append(x_exp)
            y_cham_out.append(y_exp)
            # set the throat point and expansion point
            # at the end of points group for combustion chamber.

            chamber_volume = Geo.lin_chamber_volume(x_cham_in, y_cham_in,
                                x_cham_out, y_cham_out, contour_data['depth'])
            # input the points on combustion chamber, get the volume of
            # combustion chamber.

            a_throat = np.sqrt((x_throat - x_exp) ** 2 +
                               (y_throat - y_exp) ** 2) *                     \
            contour_data['depth']
            # get the A* (cross-sectional area at throat.)

        l_cha = chamber_volume / a_throat
        # calculate the characteristic length of combustion chamber.

        t_stay = chamber_volume / ((1 / gas_parameter['density chamber']) *
                                   (gas_parameter['sonic velocity throat'] *
                                    gas_parameter['density throat'] *
                                    a_throat))
        # calculate the stay time of propellants in combustion chamber.

        self.ui.c_v.setText(str(chamber_volume))
        self.ui.g_t.setText(str(t_stay))
        self.ui.l_c.setText(str(l_cha))
        # set the results in text field.

        contour_data['chamber volume'] = chamber_volume
        contour_data['charactic chamber length'] = l_cha
        contour_data['stay time of propellant gases'] = t_stay
        self.json_data.revise(1, contour_data)
        # write the results in JSON file in the line for geometric data.

        msg = QMessageBox.information(None, "Messagebox",
                            "Analyse for Combustion Chamber is finished")
        print(msg)



    def heatflux_show(self):
        '''draw the diagrams for heat transfer.
           '''
        heat_data = self.json_data.read(2)
        # read the thermo data from JSON file.
        self.ui.heat_shell.setText(str(
            heat_data['Thermo Load on outer shell in kW']))

        self.ui.heat_spike.setText(str(
            heat_data['Thermo Load on inner core in kW']))

        self.ui.heat_total.setText(str(
            heat_data['Total Thermo Load in kW']))
        # set the values in text fields.

        n_x = len(heat_data['Inner x coordinate'])
        n_x_c = len(heat_data['Outer x coordinate'])
        x_cor = list()
        x_cor_c = list()
        for i in range(n_x):
            x_cor.append(heat_data['Inner x coordinate'][i] * 1000)
        for i in range(n_x_c):
            x_cor_c.append(heat_data['Outer x coordinate'][i] * 1000)
        # translate meter to millimeter.
        contour = self.json_data.read(1)
        n_in = len(contour['x in'])
        n_out = len(contour['x out'])
        x_in = list()
        y_in = list()
        x_out = list()
        y_out = list()
        for i in range(n_in):
            x_in.append(contour['x in'][i] * 1000)
            y_in.append(contour['y in'][i] * 1000)
        for i in range(n_out):
            x_out.append(contour['x out'][i] * 1000)
            y_out.append(contour['y out'][i] * 1000)

        figt_c, ax1t_c = plt.subplots(
            num="Convective Heat transfer coefficient", figsize=(8, 5))

        ax1t_c.plot(x_cor,
                    heat_data['Convective Heat transfer coefficient'], 'r')

        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax1t_c.set_xlabel('X [Milimeter]', fontsize=30)

        ax1t_c.set_ylabel(
            'Convective Heat transfer coefficient [W/k * m^2]',
            color='r', fontsize=20)

        ax2t_c = ax1t_c.twinx()
        ax2t_c.plot(x_in, y_in, 'k--')
        ax2t_c.plot(x_out, y_out, 'k--')
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax2t_c.set_ylabel('Y [Milimeter]', color='k', fontsize=30)


        figt, ax1t = plt.subplots(
            num="Total Heat Flux Density on Spike", figsize=(8, 5))

        ax1t.plot(x_cor, heat_data['Total Thermo flux density'], 'r')
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax1t.set_xlabel('X [Milimeter]', fontsize=30)
        ax1t.set_ylabel('Total Heat flux [W/m^2]', color = 'r', fontsize=30)
        ax2t = ax1t.twinx()
        ax2t.plot(x_in, y_in, 'k--')
        ax2t.plot(x_out, y_out, 'k--')
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax2t.set_ylabel('Y [Milimeter]', color='k', fontsize=30)


        figc, ax1c = plt.subplots(
            num="Convective Heat Flux Density on Spike", figsize=(8, 5))

        ax1c.plot(x_cor, heat_data['Convective Thermo flux density'], 'r')
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax1c.set_xlabel('X [Milimeter]', fontsize=30)
        ax1c.set_ylabel('Convective Heat flux [W/m^2]', color='r', fontsize=30)
        ax2c = ax1c.twinx()
        ax2c.plot(x_in, y_in, 'k--')
        ax2c.plot(x_out, y_out, 'k--')
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax2c.set_ylabel('Y [Milimeter]', color='k', fontsize=30)


        figr, ax1r = plt.subplots(
            num="Radition Heat Flux Density on Spike", figsize=(8, 5))

        ax1r.plot(x_cor, heat_data['Radiation Thermo flux density'], 'r')
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax1r.set_xlabel('X [Milimeter]', fontsize=30)
        ax1r.set_ylabel('Radition Heat flux [W/m^2]', color='r', fontsize=30)
        ax2r = ax1r.twinx()
        ax2r.plot(x_in, y_in, 'k--')
        ax2r.plot(x_out, y_out, 'k--')
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax2r.set_ylabel('Y [Milimeter]', color='k', fontsize=30)


        figtc, ax1tc = plt.subplots(
            num="Total Heat Flux Density on Shell", figsize=(8, 5))

        ax1tc.plot(x_cor_c,
            heat_data['Total Thermo flux density on Chamber Shell'], 'r')

        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax1tc.set_xlabel('X [Milimeter]', fontsize=30)
        ax1tc.set_ylabel('Total Heat flux [W/m^2]', color='r', fontsize=30)
        ax2tc = ax1tc.twinx()
        ax2tc.plot(x_in, y_in, 'k--')
        ax2tc.plot(x_out, y_out, 'k--')
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax2tc.set_ylabel('Y [Milimeter]', color='k', fontsize=30)

        figcc, ax1cc = plt.subplots(
            num="Convective Heat Flux Density on Shell", figsize=(8, 5))

        ax1cc.plot(x_cor_c,
            heat_data['Convective Thermo flux density on Chamber Shell'], 'r')

        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax1cc.set_xlabel('X [Milimeter]', fontsize=30)
        ax1cc.set_ylabel(
            'Convective Heat flux [W/m^2]', color='r', fontsize=30)

        ax2cc = ax1cc.twinx()
        ax2cc.plot(x_in, y_in, 'k--')
        ax2cc.plot(x_out, y_out, 'k--')
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax2cc.set_ylabel('Y [Milimeter]', color='k', fontsize=30)

        figrc, ax1rc = plt.subplots(
            num="Radition Heat Flux Density on Shell", figsize=(8, 5))

        ax1rc.plot(x_cor_c,
                   heat_data['Radiation Thermo flux density in Chamber'], 'r')

        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax1rc.set_xlabel('X [Milimeter]', fontsize=30)
        ax1rc.set_ylabel('Radition Heat flux [W/m^2]', color='r', fontsize=30)
        ax2rc = ax1rc.twinx()
        ax2rc.plot(x_in, y_in, 'k--')
        ax2rc.plot(x_out, y_out, 'k--')
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax2rc.set_ylabel('Y [Milimeter]', color='k', fontsize=30)

        figrc_c, ax1rc_c = plt.subplots(
            num="Convective Heat transfer coefficient on Shell",
            figsize=(8, 5))

        ax1rc_c.plot(x_cor_c,
        heat_data['Convective Heat transfer coefficient on Chamber Shell'],
                     'r')

        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax1rc_c.set_xlabel('X [Milimeter]', fontsize=30)

        ax1rc_c.set_ylabel(
            'Convective Heat transfer coefficient [W/k * m^2]',
            color='r', fontsize=20)

        ax2rc_c = ax1rc_c.twinx()
        ax2rc_c.plot(x_in, y_in, 'k--')
        ax2rc_c.plot(x_out, y_out, 'k--')
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        ax2rc_c.set_ylabel('Y [Milimeter]', color='k', fontsize=30)

        plt.show()


    def output_contour(self):
        '''output the contour data in txt file.
            '''
        contour_ad, _filter =                                                 \
        QtWidgets.QFileDialog.getSaveFileName(None,
                                            "Save", "./", "txt Files(*.txt)")
        contour_data = self.json_data.read(1)

        with open(contour_ad, 'w') as f:
            n = len(contour_data['x in'])
            m = len(contour_data['x out'])
            f.write('########Inner contour##########')
            f.write('\n')
            for i in range(n):
                f.write(str(contour_data['x in'][i]))
                f.write('\t')
                f.write(str(contour_data['y in'][i]))
                f.write('\n')
            f.write('########outer contour#######')
            f.write('\n')
            for i in range(m):
                f.write(str(contour_data['x out'][i]))
                f.write('\t')
                f.write(str(contour_data['y out'][i]))
                f.write('\n')
            f.write(
            '########Depth(if Depth = 0 means annular Aerospike)########')
            f.write('\n')
            f.write(str(contour_data['depth']))
            if contour_data.has_key('chamber volume'):
                f.write('\n')
                f.write('#chamber volume[cubic metre]')
                f.write('\t')
                f.write('#charactic chamber length[Meter]')
                f.write('\t')
                f.write('#stay time of propellant gases[sec]')
                f.write('\n')
                f.write(str(contour_data['chamber volume']))
                f.write('\t')
                f.write(str(contour_data['charactic chamber length']))
                f.write('\t')
                f.write(str(contour_data['stay time of propellant gases']))
                f.write('\n')
        msg = QMessageBox.information(None, "Messagebox", "The file is saved")
        print(msg)


    def output_heat(self):
        '''output the thermo data in txt file.
                    '''
        heat_ad, _filter = QtWidgets.QFileDialog.getSaveFileName(None,
                                            "Save", "./", "txt Files(*.txt)")
        heat_data = self.json_data.read(2)

        with open(heat_ad, 'w') as f:
            n = len(heat_data['Inner x coordinate'])
            m = len(heat_data['Outer x coordinate'])
            f.write('#Total Heat Transfer [kW]')
            f.write('\t')
            f.write('#Total Heat Transfer around the Spike [kW]')
            f.write('\t')
            f.write('#Total Heat Transfer around the Shell [kW]')
            f.write('\n')
            f.write(str(heat_data['Total Thermo Load in kW']))
            f.write('\t')
            f.write(str(heat_data['Thermo Load on inner core in kW']))
            f.write('\t')
            f.write(str(heat_data['Thermo Load on outer shell in kW']))
            f.write('\n')
            f.write('########Thermo flux density on inner part##########')
            f.write('\n')
            f.write('#Total Thermo flux density [W/m^2]')
            f.write('\t')
            f.write('Convective Thermo flux density [W/m^2]')
            f.write('\t')
            f.write('Radiation Thermo flux density [W/m^2]')
            f.write('\t')
            f.write('Convective Heat transfer coefficient [W/K * m^2]')
            f.write('\t')
            f.write('x coordinate [Meter]')
            f.write('\n')
            for i in range(n):
                f.write(str(heat_data['Total Thermo flux density'][i]))
                f.write('\t')
                f.write(str(heat_data['Convective Thermo flux density'][i]))
                f.write('\t')
                f.write(str(heat_data['Radiation Thermo flux density'][i]))
                f.write('\t')

                f.write(str(
                    heat_data['Convective Heat transfer coefficient'][i]))

                f.write('\t')
                f.write(str(heat_data['Inner x coordinate'][i]))
                f.write('\n')
            f.write('########Thermo flux density on outer part##########')
            f.write('\n')
            f.write('#Total Thermo flux density [W/m^2]')
            f.write('\t')
            f.write('Convective flux density [W/m^2]')
            f.write('\t')
            f.write('Radiation Thermo flux density [W/m^2]')
            f.write('\t')
            f.write('Radiation Thermo flux density [W/m^2]')
            f.write('\t')
            f.write('Convective Heat transfer coefficient [W/K * m^2]')
            f.write('\t')
            f.write('x coordinate [Meter]')
            f.write('\n')
            for i in range(m):
                f.write(str(
                    heat_data['Total Thermo flux density in Chamber'][i]))
                f.write('\t')
                f.write(str(
                    heat_data['Convective Thermo flux density in Chamber'][i]))
                f.write('\t')
                f.write(str(
                    heat_data['Radiation Thermo flux density in Chamber'][i]))
                f.write('\t')
                f.write(str(
                    heat_data['Convective Heat transfer coefficient           \
                              in Chamber'][i]))
                f.write('\t')
                f.write(str(heat_data['Outer x coordinate'][i]))
                f.write('\n')
        msg = QMessageBox.information(None, "Messagebox", "The file is saved")
        print(msg)


app = QApplication([])
stats = User_interface()
stats.ui.show()
app.exec()