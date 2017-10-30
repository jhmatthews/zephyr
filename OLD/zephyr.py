#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python -i

########################################################
#
#				JM/NSH 10/03/13: ZEPHYR.PY
#
#	models SV wind. written by NSH with some additions
#
#
########################################################
import csv, sys, os, array, warnings
import matplotlib.pyplot as plt
import numpy as np
import sv_sub as sv
import disk_sub as disk
from collections import namedtuple
import astro_help as ah
import geo
import astroconst as consts
import py_read_output as r
from pylab import *

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



nogrid_read=True

if (len(sys.argv)>1):
	fname=(sys.argv[1])
	em_fname=fname
	for i in range(len(sys.argv)):
		if sys.argv[i]=="-grid":
			zscale=float(sys.argv[i+1])
			nz=int(sys.argv[i+2])
			nogrid_read=False
			em_fname=fname+'_scale1e'+str(np.log10(zscale))
else:
	fname="test"

fname_to_open=fname+".pf"
specfname=fname+".spec"

print'ZEPHYR: models SV93 wind for parameter file', fname_to_open

sv.input_params(fname_to_open)



#################################################
# print out some simple luminosities


# frequency boundaries
# print 'frequency boundaries:',fev,f1,f2
fev=13.62/HEV
c4ev=64.4/HEV
f2500=C/(2500*1e-8)
f2=C/(800*1e-8)
f1=C/(2000*1e-8)

model_dictionaries = [r.read_pf(fname)]

# diagnostic quantities
class diags:

	def __init__(self):
		self.alpha_ox = []
		self.L_x_bal = []
		self.IP = []
		self.EM = []

diag = diags()

L_x = 10.0 ** np.arange(43,45,0.1)		# x-ray luminosity
fill = 10.0 ** np.arange(-3,1,1)		# filling factor	
alpha = np.array([0.5,1.5,3,4,10])		# acceleration exponent
omegas = np.arange(0.2,0.8,0.1)			# covering factor of wind

nparams  = len(L_x) * len(fill) * len(alpha) * len(omegas)


for ix in range(1):

	clf()

	mod = geo.model(model_dictionaries[0])

	#mod.print_lums()

	mod.lum_agn = L_x[ix]
	mod.init_disk()

	mod.init_plasma()
	mod.plasma.get_cell_props()
	mod.get_rad_properties()

	mod.plasma.get_IP(mod)

	



	#sys.exit()

	# TO DO CUMUL_IONIZ
	print "The total ionizing luminosity is ",mod.tot_ioniz," which represents ",mod.tot_nioniz," ionizing photons"

	mod.init_plasma()

	#print "The peak occurs at a radius of ",mod.r_peak,"cm from the origin"

	if mod.QSO==1:
		const_agn=mod.lum_agn / ((( (2.42e18**( mod.alpha_agn + 1.))) -  (4.84e17**( mod.alpha_agn + 1.0))) /(mod.alpha_agn + 1.0))
		QSO_ioniz=const_agn*((((1e20** (mod.alpha_agn + 1.))) -  (fev** (mod.alpha_agn + 1.))) /(mod.alpha_agn + 1.0))
		QSO_nioniz=(const_agn/H)*((((1e20** ((mod.alpha_agn-1.) + 1.))) -  (fev** ((mod.alpha_agn-1.) + 1.))) /((mod.alpha_agn-1.) + 1.0))
		L_2kev=const_agn*4.84e17**mod.alpha_agn
		Xray_Lum=const_agn*((( (1e20**( mod.alpha_agn + 1.))) -  (1e14**( mod.alpha_agn + 1.0))) /(mod.alpha_agn + 1.0))
		agn_spec=np.zeros(len(mod.disk.spec_freq))
		tot_spec=np.zeros(len(mod.disk.spec_freq))
		for i in range (len(mod.disk.spec_freq)):
			agn_spec[i]=const_agn*mod.disk.spec_freq[i]**mod.alpha_agn
			tot_spec[i]=agn_spec[i]+mod.disk.disk_spec[i]
	else:
		QSO_ioniz=0.0
		QSO_nioniz=0.0
		L_2kev=0.0

	print "Disk_ioniz=", mod.tot_ioniz
	print "QSO_ioniz=", QSO_ioniz
	print "QSO_nioniz=", QSO_nioniz
	#print "F_2500=", f2500
	print "L_2500=", mod.L_2500, disk.lnu_disk(C/(2500*1e-8),mod.mstar,mod.disk_mdot,mod.rstar,mod.disk_radmax)
	print "L_1549=", disk.lnu_disk(C/(1549*1e-8),mod.mstar,mod.disk_mdot,mod.rstar,mod.disk_radmax)
	print "L_2kev=", L_2kev
	print "L_bol=", mod.ledd*mod.edd_frac
	print "Xray_lum=", Xray_Lum

	alpha_ox = -0.3838 * np.log10(mod.L_2500/L_2kev)
	print "ALPHA_OX=",alpha_ox
	f_3560=disk.llamb_disk(3560.0*1e-8,mod.mstar,mod.disk_mdot,mod.rstar,mod.disk_radmax)*1e-8
	# print "L_3560",f_3560
	# f_3560=f_3560/(4*consts.PI*consts.PC*consts.PC*100.0) 
	# print "F_3560 (10pc)=",f_3560
	# print "Mu=",2.5*np.log10(419.6e-11/f_3560)
	# print "Mbol=",4.75-2.5*np.log10((mod.ledd*mod.edd_frac)/3.82e33)

	# first, lets 

	

'''
The total ionizing luminosity is  6.25180875082e+45  which represents  2.06971801329e+56  ionizing photons
The peak occurs at a radius of  1.23291441291e+15 cm from the origin
Disk_ioniz= 6.25180875082e+45
QSO_ioniz= 6.28260658126e+43
QSO_nioniz= 1.7703490215e+53
F_2500= 1.19917e+15
L_2500= 6.34656912965e+30 6.34373260559e+30
L_1549= 5.41404138628e+30
L_2kev= 1.18321395488e+25
L_bol= 2.81814054615e+46
Xray_lum= 7.30796621787e+43
ALPHA_OX= -2.19897278588
L_3560 1.5363311042e+42
F_3560 (10pc)= 128.375784146
Mu= -26.2141190746
Mbol= -27.4197482138
Plotting various spectra of e.g. disk, censrc...
You are defining a 100 by 100 grid.
zscale: 8.85667e+13 xscale: 4.428335e+16 rmax: 1e+19
'''

# dl = 
# slx = 


# sightline_array=[71.0,72.0,73.0,74.0,75.0,76.0,77.0,78.0,79.0,80.0,81.0,82.0,83.0,84.0,85.0,86.0,87.0,88.0,89.0]
# total_column=[]
# for sightline in sightline_array:
# 	sightline=sightline/RADIAN
# 	test_sightline=(PI/2.0)-sightline
# 	launch=1e15
# 	test_thetamin=(PI/2.0)-geo.sv_thetamin
# 	test_thetamax=(PI/2.0)-geo.sv_thetamax
# 	if (test_sightline < test_thetamin):
# 		z1=(geo.sv_rmin-launch)/(1/np.tan(test_sightline)-1/np.tan(test_thetamin))
# 		x1=launch+z1/np.tan(test_sightline)
# 		if (test_sightline < test_thetamax):
# 			z2=(geo.sv_rmax-launch)/(1/np.tan(test_sightline)-1/np.tan(test_thetamax))
# 			x2=launch+z2/np.tan(test_sightline)
# 			if (z2 > geo.wind_rmax):
# 				z2=geo.wind_rmax
# 			if (x2 > geo.wind_rmax):
# 				x2=geo.wind_rmax		
# 		else:
#  			x2=geo.wind_rmax
# 			z2=geo.wind_rmax
# 	print 'At sightline ',np.degrees(sightline),' we go from',x1,z1,' to ',x2,z2
# 	xsl=[]
# 	zsl=[]
# 	lxsl=[]
# 	lzsl=[]
# 	lsl=[]
# 	rhosl=[]
# 	col=0.0
# 	column=[]
# 	dl=np.sqrt(((z2-z1)*(z2-z1))+((x2-x1)*(x2-x1)))
# 	dl=dl/10000.0
# 	for i in range(10000):
# 		xsl.append(x1+i*dl*np.cos(test_sightline))
# 		lxsl.append(np.log10(xsl[i]))
# 		zsl.append(z1+i*dl*np.sin(test_sightline))
# 		lzsl.append(np.log10(zsl[i]))
# 		lsl.append(np.sqrt(xsl[i]*xsl[i]+zsl[i]*zsl[i]))
# 		rhosl.append((sv.sv_rho([xsl[i],0.0,zsl[i]])*RHO2NH))
# 		col=col+(rhosl[i]*dl)
# 		column.append(col)
# 	print "Cumulative hydrogen density along sightline ",np.degrees(sightline)," to origin is ",column[-1]
# 	total_column.append(column[-1])

# tot_col=np.array(total_column)
# tau=tot_col*consts.THOMPSON


