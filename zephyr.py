#!/usr/bin/env python
'''
				zephyr 

Zephyr is a program that sets up a biconical wind
model by setting up a series of streamlines with an
analytical velocity law and then integrating continuity
equations to get the density at a given point.

Usage:
	python zephyr.py [filename]
	where filename is a radiative transfer style parameter filename
'''
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
					parameter filename to read 
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
		self.sv.v_infinity = float(self.params["sv.v_infinity"]) 		# hardwired value of 6 km/s
		self.sv.sv_lambda = float(self.params["sv.mdot_r_exponent"])
		self.sv.gamma = 1.0
		self.sv.wind_mdot = float(self.params["wind.mdot"]) * MSOL / YR
		self.fill_factor = float(self.params["filling_factor"]) 

		self.sv.get_mdot_norm()



	def setup(self):

		'''
		setup some grid arrays 
		'''

		self.nx = int(self.params["Wind.dim.in.x_or_r.direction"])
		self.nz = int(self.params["Wind.dim.in.z_or_theta.direction"])

		# hardwired x and z scales at the moment, as well as "cen_or_not"
		self.x, self.z = geo.get_grid(self.compact_object.rstar/10.0, self.sv.rmin, self.nz, self.nx, self.sv.wind_rmax, 0)

		self.grid_shape = [len(self.x), len(self.z)]
		indices = np.indices(self.grid_shape)
		self.x2d = self.x[indices[0]]
		self.z2d = self.z[indices[1]]

		#print (self.x2d.shape)

		#print (self.z)

	def populate_primitive_variables(self):
		'''
		populate the grid with primitive variables (density, velocity and so on)

		IMPROVE: make more pythonic, less for-looping please!
		'''

		# these help construct a mask for in and out of wind
		rho_min = self.sv.rmin + self.z2d * np.tan (self.sv.thetamin)
		rho_max = self.sv.rmax + self.z2d * np.tan (self.sv.thetamax)

		# this masks values out of the wind
		inwind_mask = (self.x2d > rho_max) + (self.x2d < rho_min)

		print (self.x2d)

		# create masked arrays
		shape = self.x2d.shape
		self.nh = np.ma.masked_array(np.zeros(shape), mask=inwind_mask)
		self.density = np.ma.masked_array(np.zeros(shape), mask=inwind_mask)
		self.velocity = np.ma.masked_array(np.zeros(shape), mask=inwind_mask)

		for  j in range(len(self.z)):
			for i in range(len(self.x)):

				if not (inwind_mask[i,j]):

					# 3d vector, x, theta, z 
					xtest = [self.x[i], 0.0, self.z[j]]

					# get density from sv subroutines, convert to hydrogen density 
					self.density[i,j] = self.sv.sv_rho(xtest) / self.fill_factor

					n_h = self.density[i,j] * RHO2NH

					# copy to array
					self.nh[i,j] = np.log10(n_h)

					# need to add velocity here
					print ( self.density[i,j], n_h)


if __name__ == "__main__":

	if len(sys.argv) > 1:
		root = sys.argv[1]

		mod = grid(root)
		mod.setup()
		mod.populate_primitive_variables()

		import matplotlib.pyplot as plt 

		#print (mod.x, mod.z, mod.nh)
		plt.pcolormesh(mod.x2d, mod.z2d, mod.nh)
		plt.colorbar()
		plt.loglog()
		plt.xlim(1e15,1e20)
		plt.ylim(1e15,1e20)
		plt.show()

	else:
		print (__doc__)




