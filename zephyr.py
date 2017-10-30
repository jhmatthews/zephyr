import py_read_output as r 
import sv_sub as sv
import geo

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
		self.sv_thetamin = self.params["sv.diskmin"]
		self.sv_thetamax = self.params["sv.diskmax"]
		self.sv_rmin = float(self.params["sv.thetamin"]) / RADIAN
		self.sv_rmax = float(self.params["sv.thetamax"]) / RADIAN
		self.wind_rmax = float(self.params["wind.radmax"])

	def setup(self):

		self.nx = self.params["Wind.dim.in.x_or_r.direction"]
		self.nz = self.params["Wind.dim.in.z_or_theta.direction"]

		# hardwired x and z scales at the moment, as well as "cen_or_not"
		self.x, self.z = geo.get_grid(self.rstar/10.0, self.sv_rmin, nz, nx, self.rmax, 0)

		self.density = np.empty([len(x),len(z)])
		self.nh = np.empty([len(x),len(z)])


	def populate_primitive_variables(self):
		'''
		populate the grid with primitive variables (density, velocity and so on)

		IMPROVE: make more pythonic, less for looping please!
		'''
		self.nh = np.zeros([len(x), len(z)])
		self.density = np.zeros_like(self.nh)
		self.velocity = np.zeros(self.nh)

		for  j in range(len(self.z)):
			for i in range(len(self.x)):
				rho_min = self.sv_rmin + z[j] * np.tan (self.sv_thetamin)
				rho_max = self.sv_rmax + z[j] * np.tan (self.sv_thetamax)

				if (x[i] <= rho_max and x[i] >= rho_min):

					# 3d vector, x, theta, z 
					xtest = [x[i],0.0,z[j]]

					# get density from sv subroutines, convert to hydrogen density 
					self.density[i,j] = sv.sv_rho(xtest)
					n_h = density[i,j] * sv.RHO2NH

					# copy to array
					self.nh[i,j] = np.log10(n_h)

					# need to add velocity here

			else:
				# improve: out of wind, the values should be masked.
				self.density[i,j] = -99
				self.nh[i,j] = -99


if __name__ == "__main__":

	mod = grid(root)
	mod.setup()
	mod.populate_primitive_variables()

	import matplotlib.pyplot as plt 
	plt.pcolormesh(mod.x, mod.z, mod.nh)
	plt.loglog()
	plt.show()




























