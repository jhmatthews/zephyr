########################################################
#
#				JM/NSH 10/03/13: GEO.PY
#
#	gives geo structure and defines the grid for zephyr
#
########################################################

import numpy as np
import math
import sv_sub
import scipy

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


mstar=0.0
rstar=0.0
disk_mdot=0.0
wind_rmax=0.0
wind_mdot=0.0
sv_rmin=0.0
sv_rmax=0.0
sv_thetamin=0.0
sv_thetamax=0.0                     
sv_lambda=0.0
sv_v_infinity=0.0
sv_r_scale=0.0
sv_alpha=0.0   
mdot_norm = 0.0
sv_v_zero = 6e5
sv_gamma = 1.0
disk_radmax=0.0
xlog_scale=1.0e7 #set some default values for log scales
zlog_scale=1.0e7
QSO=0

def get_grid(zscale, xscale, nz, nx, rmax, cen_or_not):

	'''defines your grid in x and z coordinates- being the cen/upperboundaries of cells
	depending on the variable cen_or_not. if not 0: upperbound, if 0: centre'''

	xmaxscaled= rmax / xscale
	zmaxscaled= rmax / zscale
	dlogx = (np.log10(xmaxscaled)) / (nx - 3)
	dlogz = (np.log10(zmaxscaled)) / (nz - 3)
	x,z = np.zeros(nx), np.zeros(nz)

	for ix in range(nx):
		if cen_or_not!=0: 
			x[ix]= xscale * (10.0** (dlogx * (ix)) )	
		else:
			if ix>0: x[ix]= xscale * (10.0** (dlogx * (ix-1)) )	 	

	for iz in range(nz):
		if cen_or_not!=0:
			z[iz]= zscale * (10.0** (dlogz * (iz)) )	
		else:
			if ix>0: zscale * (10.0** (dlogz * (iz-1)) )	

	return x,z


class model:

	def __init__(self, pf_dict):

		self.mstar = pf_dict["mstar(msol)"] * MSOL 
		self.rstar = pf_dict["rstar(cm)"]
		self.disk_mdot = pf_dict["disk.mdot(msol/yr)"]*MSOL/YR
		self.wind_rmax = pf_dict["wind.radmax(cm)"]
		self.wind_mdot = pf_dict["wind.mdot(msol/yr)"]*MSOL/YR
		self.sv_rmin=pf_dict["sv.diskmin(wd_rad)"]*self.rstar
		self.sv_rmax=pf_dict["sv.diskmax(wd_rad)"]*self.rstar
		self.sv_thetamin=pf_dict["sv.thetamin(deg)"]/RADIAN
		self.sv_thetamax=pf_dict["sv.thetamax(deg)"]/RADIAN                    
		self.sv_lambda=pf_dict["sv.mdot_r_exponent"]
		self.sv_v_infinity=pf_dict["sv.v_infinity(in_units_of_vescape"]
		self.sv_r_scale=pf_dict["sv.acceleration_length(cm)"]
		self.sv_alpha=pf_dict["sv.acceleration_exponent"]  
		self.fill_factor = pf_dict["wind.filling_factor(1=smooth,<1=clumpted)"]

		# some hardwired values
		self.mdot_norm = 0.0
		self.sv_v_zero = 6e5
		self.sv_gamma = 1.0
		self.disk_radmax = pf_dict["disk.radmax(cm)"]
		self.xlog_scale = 1.0e7 #set some default values for log scales
		self.zlog_scale = 1.0e7

		self.lum_agn = pf_dict["lum_agn(ergs/s)"]
		self.alpha_agn = pf_dict["agn_power_law_index"]
		self.QSO = 1

		self.ledd=1.3e38*self.mstar/MSOL
		self.medd=self.ledd/(0.1*C*C)
		self.edd_frac=(self.disk_mdot)/self.medd
		self.omega=np.cos(self.sv_thetamin)-np.cos(self.sv_thetamax)


  		self.mdot_norm =scipy.integrate.romberg (sv_sub.sv_wind_mdot_integral, self.sv_rmin, self.sv_rmax, rtol=1e-6);



	def print_lums(self):

		print "accretion rate", self.disk_mdot
		print "mass", self.mstar
		print "The eddington luminosity for this black hole is ",self.ledd

		print "Which implies an accretion rate of ",self.medd/(MSOL/YR)

		print "This geometry is accreting at an eddington fraction of ",self.edd_frac 

		print "With these angles, the wind will subtend a fraction of ",self.omega," possible viewing angles"

		import sys
		return 0









		
	
