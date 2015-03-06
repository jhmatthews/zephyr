########################################################
#
#				JM/NSH 10/03/13: GEO.PY
#
#	gives geo structure and defines the grid for zephyr
#
########################################################

import numpy as np
import math
import sv_sub as sv
import scipy
from collections import OrderedDict
import disk_sub

RADIAN=57.2959
C=2.997925e10
MSOL=1.979e33
G=6.670e-8
RHO2NH=4.217851e23
PI=3.1415 
HEV=4.13620e-15
H=6.6262e-27
WIEN=5.879e10
YR=3.1556925e7
THOMPSON=0.66524e-24
STEFAN_BOLTZMANN=5.6696e-5

sv_v_zero = 6e5
sv_gamma = 1.0
disk_radmax=0.0
xlog_scale=1.0e7 #set some default values for log scales
zlog_scale=1.0e7



class model(object):

	def __init__(self, pf_dict):

		'''
		Initialise the variables by passing a dictionary containing parameters
		'''

		self.mstar = pf_dict["mstar(msol)"] * MSOL 
		self.rstar = pf_dict["rstar(cm)"]
		self.disk_mdot = pf_dict["disk.mdot(msol/yr)"]*MSOL/YR
		self.wind_mdot = pf_dict["wind.mdot(msol/yr)"]*MSOL/YR
		self.fill_factor = pf_dict["wind.filling_factor(1=smooth,<1=clumped)"]


		# geometrical quantities
		self.wind_rmax = pf_dict["wind.radmax(cm)"]
		self.disk_radmax = pf_dict["disk.radmax(cm)"]
		self.sv_rmin=pf_dict["sv.diskmin(wd_rad)"]*self.rstar
		self.sv_rmax=pf_dict["sv.diskmax(wd_rad)"]*self.rstar
		self.sv_thetamin=pf_dict["sv.thetamin(deg)"] / RADIAN
		self.sv_thetamax=pf_dict["sv.thetamax(deg)"] / RADIAN 

		# kinematic quantities                   
		self.sv_lambda=pf_dict["sv.mdot_r_exponent"]
		self.sv_v_infinity=pf_dict["sv.v_infinity(in_units_of_vescape"]
		self.sv_r_scale=pf_dict["sv.acceleration_length(cm)"]
		self.sv_alpha=pf_dict["sv.acceleration_exponent"]  
		

		# cell dimensions
		self.nx = int(pf_dict["Wind.dim.in.x_or_r.direction"])
		self.nz = int(pf_dict["Wind.dim.in.z_or_theta.direction"])

		# some hardwired values
		self.mdot_norm = 0.0
		self.sv_v_zero = 6e5
		self.sv_gamma = 1.0

		#set some default values for log scales
		self.xlog_scale = self.sv_rmin 
		self.zlog_scale = self.rstar / 10.0

		# agn specific stuff - lx, spectral index etc.
		self.lum_agn = pf_dict["lum_agn(ergs/s)"]
		self.alpha_agn = pf_dict["agn_power_law_index"]
		self.QSO = 1

		# some inferred quantities- Eddington quantites, covering factor
		self.ledd = 1.3e38 * self.mstar / MSOL
		self.medd = self.ledd/ (0.1 * C * C)
		self.edd_frac = (self.disk_mdot)/self.medd
		self.omega = np.cos(self.sv_thetamin)-np.cos(self.sv_thetamax)

  		self.mdot_norm =scipy.integrate.romberg (sv.sv_wind_mdot_integral, self.sv_rmin, self.sv_rmax, rtol=1e-6)
  		self.sightlines = np.array(pf_dict["angle(0=pole)"]) / RADIAN 

  	def init_disk(self):

  		self.disk = disk(self)

  		return 0



  	def init_plasma(self):

  		'''
  		initialise a plasma instance
  		'''

  		self.plasma = plasma(self)

  		return 0

  	def change_sightlines(self,sightlines):

  		self.sightlines = sightlines

	def print_lums(self):

		'''
		print some information
		'''

		print "accretion rate", self.disk_mdot
		print "mass", self.mstar
		print "The eddington luminosity for this black hole is ", self.ledd
		print "Which implies an accretion rate of ", self.medd/(MSOL/YR)
		print "This geometry is accreting at an eddington fraction of ", self.edd_frac 
		print "With these angles, the wind will subtend a fraction of ", self.omega, " possible viewing angles"
		return 0

	

	def compute_sightline_taus(self):

		'''
		compute electron scattering optical depths along sightline by integrating
		through cells. Probably a smarter way to do this
		'''

		start = (0,0,0)

		for itheta in self.sightlines:

			print itheta

		return 0

	def get_rad_properties(self):

		'''
		get radiation properties- for example, the ionizing luminosity,
		number of ionizing photons, etc 
		'''

		fev=13.62/HEV
		c4ev=64.4/HEV
		f2500=C/(2500*1e-8)
		f2=C/(800*1e-8)
		f1=C/(2000*1e-8)


		fev_bool = (self.disk.f_2d > fev)
		fuv_bool = (self.disk.f_2d > f1) * (self.disk.f_2d < f2)

		self.lum_ioniz = 2 * PI * self.disk.area_2d * fev_bool * self.disk.bnus 
		self.n_ioniz = 2 * PI * self.disk.area_2d * fev_bool * self.disk.bnus / H / self.disk.f_2d 
		self.lum_uv = 2 * PI * self.disk.area_2d * fuv_bool * self.disk.bnus

		self.tot_ioniz = np.sum(self.lum_ioniz)
		self.tot_nioniz = np.sum(self.n_ioniz)

		self.L_2500_array = 2 * PI * disk_sub.planck_nu (self.disk.tdisk,f2500) * self.disk.area
		self.L_2500 = np.sum(self.L_2500_array)

		return 0





class plasma(object):

	def __init__(self, model):

		'''
		initialise the grid
		'''

		xmaxscaled= model.wind_rmax / model.xlog_scale
		zmaxscaled= model.wind_rmax / model.zlog_scale
		dlogx = (np.log10(xmaxscaled)) / (model.nx - 3)
		dlogz = (np.log10(zmaxscaled)) / (model.nz - 3)
		x,z = np.zeros(model.nx), np.zeros(model.nz)
		xcen,zcen = np.zeros(model.nx), np.zeros(model.nz)

		for ix in range(model.nx):
			xcen[ix]= model.xlog_scale * (10.0** (dlogx * (ix)) )
			if ix>0: x[ix]= model.xlog_scale * (10.0** (dlogx * (ix-1)) )	 	

		for iz in range(model.nz):
			zcen[iz]= model.zlog_scale * (10.0** (dlogz * (iz)) )	
			if ix>0: model.zlog_scale * (10.0** (dlogz * (iz-1)) )	

		self.xcen, self.zcen = np.meshgrid(xcen,zcen)
		self.x, self.z = np.meshgrid(xcen,zcen)

		# access indexed things with nn easily one d array
		self.nn = np.indices((model.nx,model.nz)).T.reshape((model.nx*model.nz,2))

		self.rho_min = model.sv_rmin + (z * np.tan (model.sv_thetamin))
		self.rho_max = model.sv_rmax + (z * np.tan (model.sv_thetamax))

		self.inwind = (self.x > self.rho_min) * (self.x < self.rho_max)
		# if (x[i]<=rho_max and x[i]>=rho_min):


	def get_cell_props(self):

		'''
		calculate velocities and densities in cell.
		Assume Hydrogen is ionized. ne = nh
		'''

		for ix in range(len(self.x)):
			for iz in range(len(self.z)):

				if self.inwind[ix][iz]:
					pos = (self.x[ix][iz], 0, self.z[ix][iz])
					self.v = sv.sv_velocity(pos)
					self.rho = sv.sv_rho(pos)
					self.nh = RHO2NH * self.rho
				else:
					self.h = None
					self.rho = None
					self.speed,self.v = None, None

		return 0

	def get_IP(self, model):

		'''
		get the logarithm of non-attenuated ionization parameter 

		U = Q_H / (4 * PI * r * r) / nh / C
		'''

		self.IP = np.log10(model.tot_nioniz/((4 * PI * (self.x+self.z) * ((10.0**self.nh)) * C)))

		#self.IP_attenuated = 
		return 0

	def get_tau_theta(self, theta):

		'''
		get electron scattering optical depth along a sightline theta from the origin to 
		max radius of wind 
		'''

		# go to wind.radmax
		dl = model.wind_rmax / 1e5

		ray = np.arange(0, model.wind_rmax, dl)

		tau = get_tau(ray, theta, dl)

		return tau

			

	def get_tau_location(self, position):

		'''
		get electron scattering optical depth to a position in the wind from the origin
		'''

		r = np.sqrt((position[0]*position[0]) + (position[2]*position[2]))
		theta = np.arctan(position[0])

		dl = r / 1e5

		ray = np.arange(0, r, dl)

		get_tau(ray, )

		return tau

	def get_tau(self, r, theta, dl):

		'''
		get electron scattering optical depth along positions (r,theta)
		in steps dl 
		'''

		col = np.zeros(len(r))

		z = r * np.cos(theta)
		x = r * np.cos(theta)

		inwind = (x > self.rho_min) * (x < self.rho_max)

		for i in range(len(x)):

			rho_min = model.sv_rmin + (z[i] * np.tan (model.sv_thetamin))
			rho_max = model.sv_rmax + (z[i] * np.tan (model.sv_thetamax))

			if (x[i] > rho_min) and (x[i] < rho_max):

				pos = [x[i],0,z[i]]
				nh = sv_rho(pos) * RHO2NH
				col[i] = nh * dl

		column = np.sum(col)
		tau = column * THOMPSON 

		return tau





class disk(object):


	def __init__(self, model):

		'''
		initialise the disk class 

		uses physical quantities from the model to calculate the properties and 
		spectrum of the disk
		'''

		# now do some disk stuff
		NRINGS = 1000

		# get the reference temperature of the disk 
		self.tref = self.t_disk(model.mstar, model.disk_mdot, model.rstar)

		# create the annuli
		rtemp = np.linspace (model.rstar,model.disk_radmax,num=NRINGS)

		# radii of disk annuli - has length len(rtemp) - 1
		rmins = rtemp[:-1]
		rmaxes = rtemp[1:]
		self.r = 0.5 * (rmins + rmaxes)

		# units of stellar radii
		self.r_over_star = self.r / model.rstar

		# get the disk spectrum (specific intensity? TODO)
		self.spec_freq,self.disk_spec=disk_sub.spec_disk2 (1e14,1e18,model.mstar,model.disk_mdot,model.rstar,model.disk_radmax)

		# area of annuli
		self.area = PI * (rmaxes*rmaxes - rmins*rmins)

		self.tdisk= self.teff(self.tref,self.r_over_star)

		# peak frequency of disk
		self.fmax = self.tdisk * WIEN
		self.f = np.linspace(1e15,1e4/HEV, num=1e3)

		# 2d meshgrid. Needed for computed luminosities quickly and easily
		self.t_2d, self.f_2d = np.meshgrid(self.tdisk, self.f)

		# likewise with areas
		self.area_2d, self.f_2d = np.meshgrid(self.area, self.f)

		# create an array of bnu for each annulus
		# will have shape (len(tdisk), len(f))
		self.bnus = disk_sub.planck_nu (self.tdisk, self.f)




	def t_disk (self, m, mdot, r):
		t = 3. * G / (8. * PI * STEFAN_BOLTZMANN) * m * mdot / (r * r * r)
		t = pow (t, 0.25)
		return (t)

	def teff (self, t, x):
		q = (1.e0 -  (x ** -0.5e0)) / (x * x * x)
		q = t * (q ** 0.25e0)
		return (q)


