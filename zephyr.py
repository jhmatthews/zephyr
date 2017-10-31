import py_read_output as r 
import sv_sub
import sys
import geo
import numpy as np 

RADIAN=57.29598
C=2.997925e10
MSOL=1.979e33
G=6.670e-8
YR=3.1556925e7
EPSILON=1e-6
PI=3.1416
RHO2NH=4.217851e23

class compact_object:
	def __init__(self, mass, radius):
		self.mstar = mass
		self.rstar = radius

class grid():
	
	def __init__(self, root):
		'''
		Initialise a grid class for a python Radiative Transfer model.

		Parameters:
			root 	string
					parameter filename
		'''

		# read a root filename 
		self.params = r.read_pf(root)

		# get the compact object parameters
		mass = float(self.params["mstar"]) * MSOL
		rstar = float(self.params["rstar"]) 
		self.compact_object = compact_object(mass, rstar)

		# initialise a Shlosman and Vitello wind model
		self.sv = sv_sub.sv_wind(self.compact_object)

		# set the actual SV parameters from the file 
		self.sv.r_scale = self.params["sv.acceleration_length"]
		self.sv.alpha = self.params["sv.acceleration_exponent"]
		self.sv.rmin = self.params["sv.diskmin"] * self.compact_object.rstar
		self.sv.rmax = self.params["sv.diskmax"] * self.compact_object.rstar
		self.sv.thetamin = float(self.params["sv.thetamin"]) / RADIAN
		self.sv.thetamax = float(self.params["sv.thetamax"]) / RADIAN
		self.sv.wind_rmax = float(self.params["wind.radmax"])
		self.sv.v_zero = 6e5 		# hardwired value of 6 km/s
		self.sv.v_infinity = float(self.params["wind.radmax"]) 		# hardwired value of 6 km/s
		self.sv.sv_lambda = float(self.params["sv.mdot_r_exponent"])
		self.sv.gamma = 1.0
		self.sv.wind_mdot = float(self.params["wind.mdot"]) * MSOL / YR

		self.sv.get_mdot_norm()



	def setup(self):

		self.nx = int(self.params["Wind.dim.in.x_or_r.direction"])
		self.nz = int(self.params["Wind.dim.in.z_or_theta.direction"])

		# hardwired x and z scales at the moment, as well as "cen_or_not"
		self.x, self.z = geo.get_grid(self.compact_object.rstar/10.0, self.sv.rmin, self.nz, self.nx, self.sv.wind_rmax, 0)

		print (self.z)

	def populate_primitive_variables(self):
		'''
		populate the grid with primitive variables (density, velocity and so on)

		IMPROVE: make more pythonic, less for looping please!
		'''
		self.nh = np.zeros([len(self.x), len(self.z)])
		self.density = np.zeros_like(self.nh)
		self.velocity = np.zeros_like(self.nh)

		for  j in range(len(self.z)):
			for i in range(len(self.x)):
				rho_min = self.sv.rmin + self.z[j] * np.tan (self.sv.thetamin)
				rho_max = self.sv.rmax + self.z[j] * np.tan (self.sv.thetamax)

				if (self.x[i] <= rho_max and self.x[i] >= rho_min):

					# 3d vector, x, theta, z 
					xtest = [self.x[i], 0.0, self.z[j]]

					# get density from sv subroutines, convert to hydrogen density 
					self.density[i,j] = self.sv.sv_rho(xtest)
					n_h = self.density[i,j] * RHO2NH

					# copy to array
					self.nh[i,j] = np.log10(n_h)

					# need to add velocity here

			else:
				# improve: out of wind, the values should be masked.
				self.density[i,j] = -99
				self.nh[i,j] = -99


if __name__ == "__main__":

	root = sys.argv[1]

	mod = grid(root)
	mod.setup()
	mod.populate_primitive_variables()

	import matplotlib.pyplot as plt 

	print (mod.x, mod.z, mod.nh)
	plt.pcolormesh(mod.x, mod.z, mod.nh.T, vmin = 4, vmax=14)
	plt.loglog()
	plt.xlim(1e15,1e20)
	plt.ylim(1e15,1e20)
	plt.show()




























