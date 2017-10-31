########################################################
#
#				JM/NSH 10/03/13: GEO.PY
#
#	gives geo structure and defines the grid for zephyr
#
########################################################

import numpy as np

def get_grid(zscale, xscale, nz, nx, rmax, cen_or_not):

	'''defines your grid in x and z coordinates- being the cen/upperboundaries of cells
	depending on the variable cen_or_not. if not 0: upperbound, if 0: centre'''

	xmaxscaled = rmax / xscale
	zmaxscaled = rmax / zscale
	dlogx = (np.log10(xmaxscaled)) / (nx - 3)
	dlogz = (np.log10(zmaxscaled)) / (nz - 3)
	x,z = np.zeros(nx), np.zeros(nz)

	print (zscale, xscale)

	for ix in range(nx):
		if cen_or_not != 0: 
			x[ix]= xscale * (10.0** (dlogx * (ix)) )	
		else:
			if ix > 0: x[ix]= xscale * (10.0** (dlogx * (ix-1)) )	 	

	for iz in range(nz):
		if cen_or_not != 0:
			z[iz] = zscale * (10.0 ** (dlogz * (iz)) )	
		else:
			if iz > 0: z[iz] = zscale * (10.0 ** (dlogz * (iz-1)) )	

	print z
	return x,z






		
	
