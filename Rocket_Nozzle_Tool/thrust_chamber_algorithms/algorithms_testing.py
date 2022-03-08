import matplotlib.pyplot as plt
import numpy as np
import thrust_optimised_aerospike
import ideal_nozzle
import ideal_aerospike_nozzle
import conical_nozzle
import aerospike_combustion_chamber
import tictop_nozzle
import parabola_nozzle
import combustion_chamber

ideal = ideal_nozzle.IdealNozzle()
x_ideal_R, y_ideal_R = ideal.ideal_nozzle_with_R(1.4, 5, 1, 100, 0.5)
x_ideal_lambda, y_ideal_lambda = ideal.ideal_nozzle_with_lambda(1.4, 5, 1, 100, 1.5)
x_conical, y_conical = conical_nozzle.conical_nozzle_with_ratio(0.2618, 1, 0.5, 25, 50)
x_top_L, y_top_L, theta_top = parabola_nozzle.parabola_nozzle_with_length(0.5236, 1, 10, 25, 100)
print('top exit angle:', theta_top * 180 / np.pi)
x_top_theta, y_top_theta, L_top = parabola_nozzle.parabola_nozzle_with_exit_angle(0.5236, 1, 0.08727, 25, 100)
print('top_length', L_top)
x_tictop_L, y_tictop_L, theta_tictop = tictop_nozzle.tictop_nozzle_with_length(1.4, 8, 3, 0.5, 10, 1, 25, 100)
print('tictop exit angle', theta_tictop * 180 /np.pi)
x_tictop_theta, y_tictop_theta, L_tictop = tictop_nozzle.tictop_nozzle_with_exit_angle(1.4, 8, 3, 0.5, 0.08727, 1, 25, 100)
print('tictop length:', L_tictop)
x_chamber, y_chamber, v_chamber = combustion_chamber.cylinder_combustion_chamber(1, 1, 1, 2.5, 10, 1.04615)
print('chamber volume', v_chamber)



x_ideal_aerospike_lin, y_ideal_aerospike_lin, nu_lin = ideal_aerospike_nozzle.ideal_aerospike_contour(1.23, 100, 3.4, 0, 0)
print('linear total turning angle', nu_lin * 180 / np.pi)
x_ideal_aerospike_ann, y_ideal_aerospike_ann, nu_ann = ideal_aerospike_nozzle.ideal_aerospike_contour(1.23, 100, 3.4, 0, 1)
print('annular total turning angle', nu_ann * 180 / np.pi)
x_optimised_aerospike, y_optimised_aerospike = thrust_optimised_aerospike.thrust_optimised_aerospike_nozzle(3.2, -0.1047, 1.23, 80, 1)



'''plt.plot(x_ideal_R, y_ideal_R)
plt.plot(x_ideal_lambda, y_ideal_lambda)
plt.plot(x_conical, y_conical)
plt.plot(x_top_L, y_top_L)
plt.plot(x_top_theta, y_top_theta)
plt.plot(x_tictop_L, y_tictop_L)
plt.plot(x_tictop_theta, y_tictop_theta)
plt.plot(x_chamber, y_chamber)'''

plt.plot(x_ideal_aerospike_lin, y_ideal_aerospike_lin)
plt.plot(x_ideal_aerospike_ann, y_ideal_aerospike_ann)
plt.plot(x_optimised_aerospike, y_optimised_aerospike)
plt.show()