# -*- coding: iso-8859-1 -*-
"""
Purpose of this code is to provide functions for (And to execute) deriving optical parameters of cloud particles (q_abs, w_0, g, and maaaybe sigma), for given particle sizes and a particle size distribution.
"""

########################
###Import useful libraries
########################
import numpy as np
import matplotlib.pyplot as plt
import pdb
import scipy.integrate
from matplotlib.pyplot import cm
from scipy import interpolate as interp
import mie_calcs as mie#Pierrehumbert PPC port of bhmie.c
import scipy.special
import time
import cPickle as pickle

########################
###Define useful constants, all in CGS (via http://www.astro.wisc.edu/~dolan/constants.html)
########################

#Unit conversions
km2m=1.e3 #1 km in m
km2cm=1.e5 #1 km in cm
cm2km=1.e-5 #1 cm in km
cm2m=1.e-2 #1 cm in m
cm2inch=1./2.54 #1 cm in inches
amu2g=1.66054e-24 #1 amu in g
g2kg=1.e-3 #1 gram in kg
bar2atm=0.9869 #1 bar in atm
atm2bar=1./bar2atm #1 atm in bar
Pascal2bar=1.e-5 #1 Pascal in bar
bar2Pa=1.e5 #1 bar in Pascal
Pa2bar=1.e-5 #1 Pascal in bar
deg2rad=np.pi/180.
bar2barye=1.e6 #1 Bar in Barye (the cgs unit of pressure)
barye2bar=1.e-6 #1 Barye in Bar
micron2m=1.e-6 #1 micron in m
micron2cm=1.e-4 #1 micron in cm
cm2micron=1.e4 #1 cm in micron
nm2cm=1.e-7 #1 nm in cm
nm2micron=1.e-3 #1 nm in um
metricton2kg=1000. #1 metric ton in kg

#Fundamental constants
c=2.997924e10 #speed of light, cm/s
h=6.6260755e-27 #planck constant, erg/s
k=1.380658e-16 #boltzmann constant, erg/K
sigma=5.67051e-5 #Stefan-Boltzmann constant, erg/(cm^2 K^4 s)
R_earth=6371.*km2m#radius of earth in m
R_sun=69.63e9 #radius of sun in cm
AU=1.496e13#1AU in cm

#Mean molecular masses
m_co2=44.01*amu2g #co2, in g
m_h2o=18.02*amu2g #h2o, in g


##########################################################################
###Define useful functions
##########################################################################

########################
###Size parameter
########################
def size_parameter(radius, wav):
	"""
	This function computes the size parameter for Mie scattering calculations.
	It takes particle radius in microns, and wavelength of light in nm
	It returns the dimensionless size parameter 2*pi*r/lambda (Taken from Hansen & Travis 1974, pg 546)
	NOTE that this implementation of the size parameter is DIFFERENT from that done in ES237, where the relative index of refraction was also included.
	"""
	wav_um=wav*nm2micron#convert wavelength from nm to cm
	return 2.*np.pi*radius/wav_um #size parameter, dimensionless

########################
### Indices of refraction
########################

###Water ice, using Warren+1984
def indref_water_ice_warren_1984(wav):
	"""
	This function returns the index of refraction for water ice.
	The data used are taken from Warren (1984), via the compilation made available by the NASA Ames Climate Modelling Group at: http://spacescience.arc.nasa.gov/mars-climate-modeling-group/documents/Waterice_Refractive_Indicies.txt. This group's calculation of optical properties is also available at: http://spacescience.arc.nasa.gov/mars-climate-modeling-group/documents/Waterice_Scattering_Properties.txt. This forms a useful check (Evaluated for r_eff=4.0 microns, sigma=0.1, log-normal distribution)
	
	Input: wavelength in nm
	
	Output: complex index of refraction, calculated using interpolation
	For values outside the range encompassed by the data (i.e. <195 nm or >167 microns) the constant value at 195 nm and 167 um respectively is used.
	"""
	refwav_um,N_r, N_i=np.genfromtxt('./Raw_Data/ComplexRefractionIndices/Waterice_Refractive_Indicies.txt', skip_header=1, skip_footer=0, usecols=(0, 1, 2), unpack=True)#wavelength in microns, real refractive index, imaginary refractive index.
	
	#index of refraction at 195 nm, for padding
	refwav_um_minval=refwav_um[0]
	N_r_minval=N_r[0]
	N_i_minval=N_i[0]
	
	#index of refraction at 167 um, for padding
	refwav_um_maxval=refwav_um[-1]
	N_r_maxval=N_r[-1]
	N_i_maxval=N_i[-1]
	
	#Create functionalized forms of indices for linear interpolation
	N_r_func=interp.interp1d(refwav_um, N_r, kind='linear')
	N_i_func=interp.interp1d(refwav_um, N_i, kind='linear')
	
	wav_um=wav*nm2micron
	if (type(wav)==float or type(wav)==np.float64): #branch: have to treat floats and arrays separately.
		if wav_um<refwav_um_minval:
			return N_r_minval+1.j*N_i_minval
		elif wav_um>refwav_um_maxval:
			return N_r_maxval+1.j*N_i_maxval
		else:
			return N_r_func(wav_um)+1.j*N_i_func(wav_um)
	else:
		result=np.zeros(np.shape(wav_um), dtype=np.complex)
		result[wav_um<refwav_um_minval]=N_r_minval+1.j*N_i_minval
		result[wav_um>refwav_um_maxval]=N_r_maxval+1.j*N_i_maxval
		
		inds=np.where((wav_um>=refwav_um_minval) & (wav_um<=refwav_um_maxval))
		result[inds]=N_r_func(wav_um[inds])+1.j*N_i_func(wav_um[inds])
		return result

###Water ice, using Warren and Brandt 2008
def indref_water_ice_warren_2008(wav):
	"""
	This function returns the index of refraction for water ice.
	The data used are taken from Warren and Brandt (2008), which forms an update to Warren 1984. They report significant changes at every wavelength bin...though the biggest changes seem to be in weak imaginary index of refraction at 200-500 nm. Data taken from http://www.atmos.washington.edu/ice_optical_constants/
	
	Input: wavelength in nm
	
	Output: complex index of refraction, calculated using interpolation
	For values outside the range encompassed by the data (i.e. <44.3 nm or >2 meters) the constant value at 44.3 nm and 2 m respectively is used.
	"""
	refwav_um,N_r, N_i=np.genfromtxt('./Raw_Data/ComplexRefractionIndices/IOP_2008_ASCIItable.dat', skip_header=0, skip_footer=0, usecols=(0, 1, 2), unpack=True)#wavelength in microns, real refractive index, imaginary refractive index.
	
	#index of refraction at 195 nm, for padding
	refwav_um_minval=refwav_um[0]
	N_r_minval=N_r[0]
	N_i_minval=N_i[0]
	
	#index of refraction at 167 um, for padding
	refwav_um_maxval=refwav_um[-1]
	N_r_maxval=N_r[-1]
	N_i_maxval=N_i[-1]
	
	#Create functionalized forms of indices for linear interpolation
	N_r_func=interp.interp1d(refwav_um, N_r, kind='linear')
	N_i_func=interp.interp1d(refwav_um, N_i, kind='linear')
	
	wav_um=wav*nm2micron
	if (type(wav)==float or type(wav)==np.float64): #branch: have to treat floats and arrays separately.
		if wav_um<refwav_um_minval:
			return N_r_minval+1.j*N_i_minval
		elif wav_um>refwav_um_maxval:
			return N_r_maxval+1.j*N_i_maxval
		else:
			return N_r_func(wav_um)+1.j*N_i_func(wav_um)
	else:
		result=np.zeros(np.shape(wav_um), dtype=np.complex)
		result[wav_um<refwav_um_minval]=N_r_minval+1.j*N_i_minval
		result[wav_um>refwav_um_maxval]=N_r_maxval+1.j*N_i_maxval
		
		inds=np.where((wav_um>=refwav_um_minval) & (wav_um<=refwav_um_maxval))
		result[inds]=N_r_func(wav_um[inds])+1.j*N_i_func(wav_um[inds])
		return result

###CO2 Ice, using Hansen (2005)...er, compilation of Pierrehumbert
def indref_co2_ice_ppc(wav):
	"""
	This function returns the index of refraction for water ice.
	The data used are nominally taken from Hansen (2005), which extends from 174 nm to 1.8 microns. However, those data are absorption cross-sections only, and do not include real index of refraction data.
	
	For the moment, we use the compendium provided by Pierrehumbert in Principles of Planetary Climate at http://geosci.uchicago.edu/~rtp1/PrinciplesPlanetaryClimate/Data/WorkbookDatasets/Chapter5Data/co2i4a.rfi.txt. These appear to extend from 52 nm to 0.2 m
	
	Input: wavelength in nm
	Output: complex index of refraction, calculated using interpolation
	
	For values outside the range encompassed by the data (i.e. <52 nm or >0.2 meters) the constant value at 52 nm and 0.2 m respectively is used.
	"""
	refwav_um,N_r, N_i=np.genfromtxt('./Raw_Data/ComplexRefractionIndices/co2i4a.rfi.txt', skip_header=4, skip_footer=0, usecols=(0, 1, 2), unpack=True)#wavelength in microns, real refractive index, imaginary refractive index.
	
	#index of refraction at 195 nm, for padding
	refwav_um_minval=refwav_um[0]
	N_r_minval=N_r[0]
	N_i_minval=N_i[0]
	
	#index of refraction at 167 um, for padding
	refwav_um_maxval=refwav_um[-1]
	N_r_maxval=N_r[-1]
	N_i_maxval=N_i[-1]
	
	#Create functionalized forms of indices for linear interpolation
	N_r_func=interp.interp1d(refwav_um, N_r, kind='linear')
	N_i_func=interp.interp1d(refwav_um, N_i, kind='linear')
	
	wav_um=wav*nm2micron
	
	if (type(wav)==float or type(wav)==np.float64): #branch: have to treat floats and arrays separately.
		if wav_um<refwav_um_minval:
			return N_r_minval+1.j*N_i_minval
		elif wav_um>refwav_um_maxval:
			return N_r_maxval+1.j*N_i_maxval
		else:
			return N_r_func(wav_um)+1.j*N_i_func(wav_um)
	else:
		result=np.zeros(np.shape(wav_um), dtype=np.complex)
		result[wav_um<refwav_um_minval]=N_r_minval+1.j*N_i_minval
		result[wav_um>refwav_um_maxval]=N_r_maxval+1.j*N_i_maxval
		
		inds=np.where((wav_um>=refwav_um_minval) & (wav_um<=refwav_um_maxval))
		result[inds]=N_r_func(wav_um[inds])+1.j*N_i_func(wav_um[inds])
		return result

###Dust using data from Wolff+2009 via http://spacescience.arc.nasa.gov/mars-climate-modeling-group/brief.html
def indref_dust_wolff_2009(wav):
	"""
	This function returns the index of refraction for dust.
	The data used are nominally taken from Wolff (2009), which extends from 236 nm to 98.5423 micron. There are only 5 points in our wavelength region of interest.
	
	For the moment, we use the compendium provided by the Ames GCM group at http://spacescience.arc.nasa.gov/mars-climate-modeling-group/brief.html. They assume an r_eff of 1.5 micron and a var_eff of 0.5, lognormal.
	
	Input: wavelength in nm
	Output: complex index of refraction, calculated using interpolation
	
	For values outside the range encompassed by the data (i.e. <236 nm or >98.5423 micron) the constant value at 52 nm and 98.5423 micron respectively is used.
	"""
	refwav_um,N_r, N_i=np.genfromtxt('./Raw_Data/ComplexRefractionIndices/Dust_Refractive_Indicies.txt', skip_header=1, skip_footer=0, usecols=(0, 1, 2), unpack=True)#wavelength in microns, real refractive index, imaginary refractive index.
	
	#index of refraction at 263 nm, for padding
	refwav_um_minval=refwav_um[0]
	N_r_minval=N_r[0]
	N_i_minval=N_i[0]
	
	#index of refraction at 98.5423 um, for padding
	refwav_um_maxval=refwav_um[-1]
	N_r_maxval=N_r[-1]
	N_i_maxval=N_i[-1]
	
	#Create functionalized forms of indices for linear interpolation
	N_r_func=interp.interp1d(refwav_um, N_r, kind='linear')
	N_i_func=interp.interp1d(refwav_um, N_i, kind='linear')
	
	wav_um=wav*nm2micron
	
	if (type(wav)==float or type(wav)==np.float64): #branch: have to treat floats and arrays separately.
		if wav_um<refwav_um_minval:
			return N_r_minval+1.j*N_i_minval
		elif wav_um>refwav_um_maxval:
			return N_r_maxval+1.j*N_i_maxval
		else:
			return N_r_func(wav_um)+1.j*N_i_func(wav_um)
	else:
		result=np.zeros(np.shape(wav_um), dtype=np.complex)
		result[wav_um<refwav_um_minval]=N_r_minval+1.j*N_i_minval
		result[wav_um>refwav_um_maxval]=N_r_maxval+1.j*N_i_maxval
		
		inds=np.where((wav_um>=refwav_um_minval) & (wav_um<=refwav_um_maxval))
		result[inds]=N_r_func(wav_um[inds])+1.j*N_i_func(wav_um[inds])
		return result


###Dust using data from Wolff+2009 via http://spacescience.arc.nasa.gov/mars-climate-modeling-group/brief.html, plus Pang+Ajello 1977
def indref_dust_wolff_pang_ajello(wav):
	"""
	This function returns the index of refraction for dust.
	The data used are nominally taken from Wolff (2009), which extends from 236 nm to 98.5423 micron. There are only 5 points in our wavelength region of interest.
	The shortwave data are padded with n_i from Pang and Ajello (1977). Zureck (1978) uses n_r that are flat but ~1.8-1.9. But, these are systematically ~.4 above the Wolff values at area of overlap between the two datasets. To join them, we offset these n_r down by 0.4
	
	The compendium provided by the Ames GCM group at http://spacescience.arc.nasa.gov/mars-climate-modeling-group/brief.html. They assume an r_eff of 1.5 micron and a var_eff of 0.5, lognormal.
	
	Input: wavelength in nm
	Output: complex index of refraction, calculated using interpolation
	
	For values outside the range encompassed by the data (i.e. <194 nm or >98.5423 micron) the constant value at 194 nm and 98.5423 micron respectively is used.
	"""
	refwav_um,N_r, N_i=np.genfromtxt('./Raw_Data/ComplexRefractionIndices/Dust_Refractive_Indicies_Wolff_Pang_Ajello.txt', skip_header=1, skip_footer=0, usecols=(0, 1, 2), unpack=True)#wavelength in microns, real refractive index, imaginary refractive index.
	
	#index of refraction at 194 nm, for padding
	refwav_um_minval=refwav_um[0]
	N_r_minval=N_r[0]
	N_i_minval=N_i[0]
	
	#index of refraction at 98.5423 um, for padding
	refwav_um_maxval=refwav_um[-1]
	N_r_maxval=N_r[-1]
	N_i_maxval=N_i[-1]
	
	#Create functionalized forms of indices for linear interpolation
	N_r_func=interp.interp1d(refwav_um, N_r, kind='linear')
	N_i_func=interp.interp1d(refwav_um, N_i, kind='linear')
	
	wav_um=wav*nm2micron
	
	if (type(wav)==float or type(wav)==np.float64): #branch: have to treat floats and arrays separately.
		if wav_um<refwav_um_minval:
			return N_r_minval+1.j*N_i_minval
		elif wav_um>refwav_um_maxval:
			return N_r_maxval+1.j*N_i_maxval
		else:
			return N_r_func(wav_um)+1.j*N_i_func(wav_um)
	else:
		result=np.zeros(np.shape(wav_um), dtype=np.complex)
		result[wav_um<refwav_um_minval]=N_r_minval+1.j*N_i_minval
		result[wav_um>refwav_um_maxval]=N_r_maxval+1.j*N_i_maxval
		
		inds=np.where((wav_um>=refwav_um_minval) & (wav_um<=refwav_um_maxval))
		result[inds]=N_r_func(wav_um[inds])+1.j*N_i_func(wav_um[inds])
		return result

###Dust using data from Ockert-Bell+1997
def indref_dust_ockertbell_1997(wav):
	"""
	This function returns the index of refraction for dust.
	The data used are transcribed from Table 2 of Ockert-Bell et al (1997). These data were used by, e.g., Patel+2002, but some authors, e.g. Korablev+2005, consider them unreliable in the UV/blue. 
	
	These data were formally derived assuming a mix of randomly oriented nonspherical particles for large size parameter (i.e. the regime we work in, normally). We don't have this capability, so we just assume Mie theory.\
	
	The data go up to 4.15 micron, but we only transcribe to 1.015 because we don't need any more. However, if need be the additional data can be trivially transcribed.
	
	Input: wavelength in nm
	Output: complex index of refraction, calculated using interpolation
	
	For values outside the range encompassed by the data (i.e. <210 nm or >4.15 micron) the constant value at 210 nm and 4.15 micron respectively is used.
	"""
	refwav_um=np.array([0.21, 0.30, 0.35, 0.40, 0.50, 0.60, 0.67, 0.70, 0.80, 1.015]) #wavelength in microns
	N_r=np.array([1.47, 1.48, 1.50, 1.51, 1.52, 1.51, 1.51, 1.51, 1.50, 1.50]) #real refractive index
	N_i=np.array([0.008, 0.038, 0.039, 0.034, 0.011, 0.004, 0.003, 0.003, 0.003, 0.003]) #imaginary refractive index.
	
	#index of refraction at 210 nm, for padding
	refwav_um_minval=refwav_um[0]
	N_r_minval=N_r[0]
	N_i_minval=N_i[0]
	
	#index of refraction at 4.15 um, for padding
	refwav_um_maxval=refwav_um[-1]
	N_r_maxval=N_r[-1]
	N_i_maxval=N_i[-1]
	
	#Create functionalized forms of indices for linear interpolation
	N_r_func=interp.interp1d(refwav_um, N_r, kind='linear')
	N_i_func=interp.interp1d(refwav_um, N_i, kind='linear')
	
	wav_um=wav*nm2micron
	
	if (type(wav)==float or type(wav)==np.float64): #branch: have to treat floats and arrays separately.
		if wav_um<refwav_um_minval:
			return N_r_minval+1.j*N_i_minval
		elif wav_um>refwav_um_maxval:
			return N_r_maxval+1.j*N_i_maxval
		else:
			return N_r_func(wav_um)+1.j*N_i_func(wav_um)
	else:
		result=np.zeros(np.shape(wav_um), dtype=np.complex)
		result[wav_um<refwav_um_minval]=N_r_minval+1.j*N_i_minval
		result[wav_um>refwav_um_maxval]=N_r_maxval+1.j*N_i_maxval
		
		inds=np.where((wav_um>=refwav_um_minval) & (wav_um<=refwav_um_maxval))
		result[inds]=N_r_func(wav_um[inds])+1.j*N_i_func(wav_um[inds])
		return result


########################
###Compute optical parameters (sigma, w_0, g) as a function of wavelength, particle size, and index of refraction
########################

def get_optical_parameters_mie_ppc(radius, wav, n_medium, n_particle, nterm, parameter):
	"""
	Calculates optical parameters sigma, w_0, g as a function of particle radius (micron), wavelength (nm), index of refraction of the surrounding medium (assumed close to 1, real), and index of refraction of the particle (complex), and number of terms in expansion
	
	Uses Mie theory, and especially the Mie scattering implementation of Pierrehumbert from Principles of Planetary Climate, which is a port of bhmie.c
	"""
	sizeparams=size_parameter(radius, wav)
	x=sizeparams*n_medium #The Pierrehumbert code requires size paramters that have the index of refraction of the medium included.
	refrel=n_particle/n_medium #relative index of refraction
	
	if (type(x)==float or type(x)==np.float64): #branch: have to treat floats and arrays separately.
		stuff=mie.bhmie(x, refrel, nterm) #only need 2 phase functions because only going to 2-stream approx
		qabs=stuff[0]
		qsca=stuff[1]
		g=stuff[2]
	else:
		qabs=np.zeros(np.shape(x))
		qsca=np.zeros(np.shape(x))
		g=np.zeros(np.shape(x))
		
		for ind in range(0, len(x)):
			stuff=mie.bhmie(x[ind], refrel[ind], nterm)
			qabs[ind]=stuff[0]
			qsca[ind]=stuff[1]
			g[ind]=stuff[2]			
	
	w_0=qsca/(qabs+qsca) #single-scattering albedo
	#return qsca, w_0, g
	if parameter=='qsca':
		return qsca
	if parameter=='w_0':
		return w_0
	if parameter=='g':
		return g
	if parameter=='all':
		return qsca, w_0, g

########################
###Compute different distributions (modified gamma and lognormal)
########################

def modified_gamma(r, r_eff, var_eff):
	"""
	Returns shape of modified gamma distribution from Hansen & Travis 1974 eqn 2.56. NOTE: Not normalized! If used directly, needs to be in the denominator. Requires var_eff<0.5
	"""
	constant=(r_eff*var_eff)**((2.*var_eff-1.)/var_eff)/scipy.special.gamma((1.-2*var_eff)/var_eff)
	return constant*r**((1.-3.*var_eff)/var_eff)*np.exp(-r/(r_eff*var_eff))

##Check that this normalizes to 1...huh, seems to!
#print scipy.integrate.quad(modified_gamma, 1., 100., args=(10., 0.1), epsabs=0, epsrel=1.e-5)[0]

def lognormal(r, r_eff, var_eff):
	"""
	Returns shape of log normal distribution from Hansen & Travis 1974 eqn 2.60, with material from pgs 557-8.
	"""
	
	r_g=r_eff/(1.+var_eff)**2.5
	sigma_g=np.sqrt(np.log(1.+var_eff))
	
	return (1./(sigma_g*(2.*np.pi)**0.5))*(1./r)*np.exp(-(np.log(r)-np.log(r_g))**2/(2.*sigma_g**2))

##Check that this normalizes to 1...and, seems to!
#print scipy.integrate.quad(lognormal, 1., 100., args=(10., 0.1), epsabs=0, epsrel=1.e-5)[0]

########################
###Put it all together: calculate optical parameters as a function of wavelength for a given distribution with specified r_eff, var_eff
########################

def get_ice_particle_optical_parameters(wav, r_eff, var_eff, indreffunc, dist):
	"""
	Compute the average values using numerical quadrature, and weighting by n(r)*pi*r**2, as shown for Q_sca in equation 2.58 of Hansen and Travis
	"""
	#indreffunc=indref_water_ice_warren_2008
	#dist=modified_gamma #set distribution to use
	
	def denom(r, r_eff, var_eff):
		return np.pi*r**2.*dist(r, r_eff, var_eff)
	def num_r(r, r_eff, var_eff):
		return r*np.pi*r**2.*dist(r, r_eff, var_eff)
	def num_qsca(r, r_eff, var_eff, wavelength, n_particle):
		return get_optical_parameters_mie_ppc(r, wavelength, n_medium, n_particle, nterm, 'qsca')*np.pi*r**2.*dist(r, r_eff, var_eff)
	def num_w_0(r, r_eff, var_eff,wavelength, n_particle):
		return get_optical_parameters_mie_ppc(r, wavelength, n_medium, n_particle, nterm, 'w_0')*np.pi*r**2.*dist(r, r_eff, var_eff)
	def num_g(r, r_eff, var_eff, wavelength, n_particle):
		return get_optical_parameters_mie_ppc(r, wavelength, n_medium, n_particle, nterm, 'g')*np.pi*r**2.*dist(r, r_eff, var_eff)
	
	sigma_mean=np.zeros(np.shape(wav))
	qsca_mean=np.zeros(np.shape(wav))
	w_0_mean=np.zeros(np.shape(wav))
	g_mean=np.zeros(np.shape(wav))
	
	lowerbound=r_eff*10.**(-10.*var_eff) #this is a very arbitrary formulation for the lower bound...
	upperbound=r_eff*10.**(4.*var_eff) #also super arbitrary for the lower bound...
	
	n_medium=1. #Set index of refraction of external medium, approximate to vaccuum for initial tests
	nterm=2 #Use 2 angles in Mie sum
	
	r_mean=scipy.integrate.quad(num_r, lowerbound, upperbound, args=(r_eff, var_eff), epsabs=0, epsrel=1.e-6)[0]/scipy.integrate.quad(denom, lowerbound, upperbound, args=(r_eff, var_eff), epsabs=0, epsrel=1.e-6)[0]
	print 'r_mean', r_mean
	G=scipy.integrate.quad(denom, lowerbound, upperbound, args=(r_eff, var_eff), epsabs=0, epsrel=1.e-6)[0]
	#print G/(np.pi*r_eff**2.) 
	start_time=time.time()
	for ind in range(0, len(wav)):
		wavelength=wav[ind]
		n_particle=indreffunc(wavelength) #wav in nm here...
		
		qsca_mean[ind]=scipy.integrate.quad(num_qsca, lowerbound, upperbound, args=(r_eff, var_eff, wavelength, n_particle), epsabs=0, epsrel=1.e-2, limit=200)[0]/G
		w_0_mean[ind]=scipy.integrate.quad(num_w_0, lowerbound, upperbound, args=(r_eff, var_eff, wavelength, n_particle), epsabs=0, epsrel=1.e-2, limit=200)[0]/G
		g_mean[ind]=scipy.integrate.quad(num_g, lowerbound, upperbound, args=(r_eff, var_eff, wavelength, n_particle), epsabs=0, epsrel=1.e-2, limit=200)[0]/G
		
	sigma_mean=qsca_mean*G/w_0_mean #convert the scattering cross-section to the total extinction cross-section (used to calculate optical depth)
	####sigma_mean[ind]=qsca_mean[ind]*G #this is really the scattering cross-section only...will eventually need to fix and rerun with the above code.

		#print 'Done with', wavelength
	print time.time()-start_time, 'seconds' #Runs in about .1 second for 12 x 60.
	return sigma_mean, w_0_mean, g_mean, qsca_mean
	
#########################
####Function to save output
#########################

def generate_cloud_optical_params(wav, r_eff, var_eff, indreffunc, dist, filename):
	sigma, w_0, g, qsca=get_ice_particle_optical_parameters(wav, r_eff, var_eff, indreffunc, dist)
	
	#save results in pickle
	f=open(filename+'.pickle', 'w')
	pickle.dump([wav, sigma, w_0, g, qsca], f)
	f.close()

	#save results in txt file
	toprint=np.zeros([len(sigma),5])
	toprint[:,0]=wav
	toprint[:,1]=sigma
	toprint[:,2]=w_0
	toprint[:,3]=g
	toprint[:,4]=qsca
	np.savetxt(filename+'.dat', toprint, delimiter='		', fmt='%1.7e', newline='\n', header='Wavelength (nm)	Total Extinction Cross-Section (micron**2)	w_0	g	Q_Scattering')
	
	
#########################
####Generate output
#########################
wav=np.arange(99.5, 500.5, step=0.1)
#wav=np.arange(90., 510., step=10.)

#generate_cloud_optical_params(wav, 1., 0.1, indref_co2_ice_ppc, lognormal, './ParticulateOpticalParameters/cloud_co2_reff1_vareff0p1_lognormal')
#generate_cloud_optical_params(wav, 1., 0.1, indref_water_ice_warren_2008, lognormal, './ParticulateOpticalParameters/cloud_h2o_reff1_vareff0p1_lognormal')
#generate_cloud_optical_params(wav, 10., 0.1, indref_co2_ice_ppc, lognormal, './ParticulateOpticalParameters/cloud_co2_reff10_vareff0p1_lognormal')
#generate_cloud_optical_params(wav, 10., 0.1, indref_water_ice_warren_2008, lognormal, './ParticulateOpticalParameters/cloud_h2o_reff10_vareff0p1_lognormal')
#generate_cloud_optical_params(wav, 100., 0.1, indref_co2_ice_ppc, lognormal, './ParticulateOpticalParameters/cloud_co2_reff100_vareff0p1_lognormal')
#generate_cloud_optical_params(wav, 100., 0.1, indref_water_ice_warren_2008, lognormal, './ParticulateOpticalParameters/cloud_h2o_reff100_vareff0p1_lognormal')

generate_cloud_optical_params(wav, 1.5, 0.5, indref_dust_wolff_2009, lognormal, './ParticulateOpticalParameters/dust_wolff_reff1p5_vareff0p5_lognormal') #they technically allowed for nonspherical particles; we ignore this because we're not doing a detailed retrieval, we just want OOM estimates
generate_cloud_optical_params(wav, 1.5, 0.5, indref_dust_wolff_pang_ajello, lognormal, './ParticulateOpticalParameters/dust_wolff_pangajello_reff1p5_vareff0p5_lognormal') #they technically allowed for nonspherical particles; we ignore this because we're not doing a detailed retrieval, we just want OOM estimates
#####generate_cloud_optical_params(wav, 1.85, 0.51, indref_dust_ockertbell_1997, lognormal, './ParticulateOpticalParameters/dust_ockertbell_reff1p85_vareff0p51_lognormal') #they technically allowed for nonspherical particles; we ignore this because we're not doing a detailed retrieval, we just want OOM estimates


plt.show()
