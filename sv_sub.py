#!/usr/bin/env python 
########################################################
#
#				SV_SUB.PY
#
#	NSH subroutines related to Shlosman Vitello wind
#
########################################################

import sys
import numpy as np
import pylab
import matplotlib.pyplot as plt
import scipy.integrate 
import scipy.optimize
from collections import namedtuple
import py_read_output as r 

RADIAN=57.29598
C=2.997925e10
MSOL=1.979e33
G=6.670e-8
YR=3.1556925e7
EPSILON=1e-6
PI=3.1416
RHO2NH=4.217851e23

class sv_wind:
	def __init__(self, my_object):
		self.object = my_object

	def get_mdot_norm(self):
		self.mdot_norm = scipy.integrate.romberg (self.sv_wind_mdot_integral, self.rmin, self.rmax, rtol=1e-6)

	def sv_velocity (self, x):
		'''
		Get the wind velocity at a given position vector x 
		'''
		v = np.zeros(3)
		zzz = v_escape = -99.
		rzero = self.sv_find_wind_rzero (x)
		theta = self.sv_theta_wind (rzero)
		r = np.sqrt(x[0] * x[0] + x[1] * x[1])
		ldist = np.sqrt ((r - rzero) * (r - rzero) + x[2] * x[2])
		vl = self.v_zero
		if (ldist > 0):
			zzz = pow (ldist / self.r_scale, self.alpha)
			if (rzero < self.object.rstar):
				v_escape = np.sqrt (2. * G * self.object.mstar / self.object.rstar);
			else:
				v_escape = np.sqrt (2. * G * self.object.mstar / rzero);
			vl = self.v_zero + (self.v_infinity * v_escape - self.v_zero) * zzz / (1. + zzz)
		
		v[0] = vl * np.sin (theta)
		if (r > 0):
			v[1] = np.sqrt (G * self.object.mstar * rzero) / r
		else:
			v[1] = 0

		v[2] = vl * np.cos (theta)

		if (x[2] < 0):		
			v[2] *= (-1)

		if (x[1] != 0.0):
			self.project_from_cyl_xyz (x, v, xtest)
			v = self.stuff_v (xtest)

		speed = (np.sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))

		return (speed, v)



	def sv_rho (self, x):

		speed,v = self.sv_velocity (x)
		rzero = self.sv_find_wind_rzero (x)
		theta = self.sv_theta_wind (rzero)
		r = np.sqrt (x[0] * x[0] + x[1] * x[1])
		ldist = np.sqrt ((r - rzero) * (r - rzero) + x[2] * x[2])
		dmdot_da = self.wind_mdot * (rzero**(self.sv_lambda)) * np.cos (theta) / self.mdot_norm / 2.
		if (rzero > self.rmax):
			rzero = self.rmax;
		dtheta_drzero =(self.sv_theta_wind (rzero) - self.sv_theta_wind ((1. - EPSILON) * rzero)) / (EPSILON * rzero)
		dr_drzero = 1. + ldist * dtheta_drzero / np.cos (theta)

		rho = rzero * dmdot_da / (dr_drzero * r * v[2])
		return (rho)



	def sv_find_wind_rzero (self, p):

		z = np.fabs (p[2])
		if (z == 0):
			x = (np.sqrt (p[0] * p[0] + p[1] * p[1]))
			return (x)
		zero_p = self.sv_zero_init (p)	
		rho_min = self.rmin + z * np.tan (self.thetamin)
		rho_max = self.rmax + z * np.tan (self.thetamax)
		rho = np.sqrt (p[0] * p[0] + p[1] * p[1])
		if (rho <= rho_min):
			x = self.rmin * rho / rho_min
			return (x)
		if (rho >= rho_max):
			x = self.rmax + rho - rho_max
			return (x)
		x = scipy.optimize.brentq (self.sv_zero_r, self.rmin, self.rmax, args=zero_p, xtol=100.)
		return (x)



	def sv_zero_init (self, p):
		zero_p= self.stuff_v (p)
		zero_p[2] = np.fabs (zero_p[2])
		return (zero_p)


	def stuff_v (self, vin):
		vout=np.zeros(3)
		vout[0] = vin[0];
		vout[1] = vin[1];
		vout[2] = vin[2];
		return (vout)

	def sv_theta_wind (self, r):
		if (r <= self.rmin):
			return (np.arctan (np.tan (self.thetamin * r / self.rmin)));
		if (r >= self.rmax):
			return (self.thetamax);
		theta = self.thetamin + (self.thetamax - self.thetamin) *  ((r - self.rmin) / (self.rmax - self.rmin)**self.gamma)
		return (theta);

	def sv_wind_mdot_integral (self, r):
		x = 2 * PI *  ( r**( self.sv_lambda + 1.)) * np.cos (self.sv_theta_wind (r))
		return (x)

	def sv_zero_r (self, r, zero_p):  
		theta = self.sv_theta_wind (r);
		rho = np.sqrt (zero_p[0] * zero_p[0] + zero_p[1] * zero_p[1])
		rho_guess = r + np.tan (theta) * zero_p[2]
		return (rho_guess - rho)


