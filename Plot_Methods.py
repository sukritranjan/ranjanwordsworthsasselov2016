# -*- coding: iso-8859-1 -*-
"""
Purpose of this code is to plot the figures from the Methods section of our paper. 
"""
########################
###Import useful libraries
########################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pdb
import cookbook
from matplotlib.pyplot import cm
import cPickle as pickle

########################
###Define useful constants, all in CGS (via http://www.astro.wisc.edu/~dolan/constants.html)
########################

#Unit conversions
km2m=1.e3 #1 km in m
km2cm=1.e5 #1 km in cm
cm2km=1.e-5 #1 cm in km
cm2inch=1./2.54 #1 cm in inches
amu2g=1.66054e-24 #1 amu in g
bar2atm=0.9869 #1 bar in atm
Pascal2bar=1.e-5 #1 Pascal in bar
Pa2bar=1.e-5 #1 Pascal in bar
bar2Pa=1.e5 #1 bar in Pascal
deg2rad=np.pi/180.
bar2barye=1.e6 #1 Bar in Barye (the cgs unit of pressure)
barye2bar=1.e-6 #1 Barye in Bar
micron2m=1.e-6 #1 micron in m
micron2cm=1.e-4 #1 micron in cm
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

#Mars parameters
g=371. #surface gravity of Mars, cm/s**2, from: http://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html

deg2rad=np.pi/180. #1 degree in radian


########################
###Which plots to generate?
########################
plot_mie_parameters=True #plot the Mie parameters used in our calculation
plot_atmprofile_z=True #Plot the sample TP profiles used in our calculation
plot_atmprofile_tp=True #Plot the sample TP profiles used in our calculation

########################
###
########################
if plot_mie_parameters:
	"""
	Plot the Mie parameters used in our study
	"""
	########################
	###Read in Models
	########################
	wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	qsca_dict={} #initialize dict to hold Q_sca
	w_0_dict={} #initialize dict to hold w_0 (single scattering albedo)
	g_dict={} #initialize dict to hold asymmetry parameters
	file_list=np.array(['./ParticulateOpticalParameters/cloud_co2_reff1_vareff0p1_lognormal.pickle','./ParticulateOpticalParameters/cloud_co2_reff10_vareff0p1_lognormal.pickle','./ParticulateOpticalParameters/cloud_co2_reff100_vareff0p1_lognormal.pickle','./ParticulateOpticalParameters/cloud_h2o_reff1_vareff0p1_lognormal.pickle','./ParticulateOpticalParameters/cloud_h2o_reff10_vareff0p1_lognormal.pickle','./ParticulateOpticalParameters/cloud_h2o_reff100_vareff0p1_lognormal.pickle','./ParticulateOpticalParameters/dust_wolff_reff1p5_vareff0p5_lognormal.pickle'])
	
	label_list=np.array(['$r_{eff}=1 \mu m$\n$var_{eff}=0.1$', '$r_{eff}=10 \mu m$\n$var_{eff}=0.1$', '$r_{eff}=100 \mu m$\n$var_{eff}=0.1$', '$r_{eff}=1 \mu m$\n$var_{eff}=0.1$', '$r_{eff}=10 \mu m$\n$var_{eff}=0.1$', '$r_{eff}=100 \mu m$\n$var_{eff}=0.1$', '$r_{eff}=1.5 \mu m$\n$var_{eff}=0.5$'])
	for picklefile in file_list:
		f=open(picklefile, 'r')
		wav, sigma, w_0, g, qsca=pickle.load(f) #units: nm, microns**2, dimless, dimless, dimless
		sigma_cgs=sigma*(micron2cm)**2 #convert XC from microns**2 to cm**2		
		
		wav_dict[picklefile]=wav
		qsca_dict[picklefile]=qsca
		w_0_dict[picklefile]=w_0
		g_dict[picklefile]=g
	
	########################
	###Plot
	########################
	fig, ax=plt.subplots(3,3, figsize=(8,7.), sharex=True, sharey=False)
	markersizeval=5.
	colorlist=np.array(['red', 'black', 'blue'])
	
	for ind in range(0, 3):
		label=label_list[ind]
		elt=file_list[ind]
		ax[0,0].set_title('CO$_2$ Ice')
		ax[0,0].plot(wav_dict[elt], qsca_dict[elt]/w_0_dict[elt], linewidth=1, color=colorlist[ind], label=label)
		ax[1,0].plot(wav_dict[elt], w_0_dict[elt], linewidth=1, color=colorlist[ind], label=label)
		ax[2,0].plot(wav_dict[elt], g_dict[elt],  linewidth=1, color=colorlist[ind], label=label)
		
	for ind in range(3, 6):
		label=label_list[ind]
		elt=file_list[ind]
		ax[0,1].set_title('H$_2$O Ice')
		ax[0,1].plot(wav_dict[elt], qsca_dict[elt]/w_0_dict[elt],linewidth=1, color=colorlist[ind-3], label=label)
		ax[1,1].plot(wav_dict[elt], w_0_dict[elt], linewidth=1, color=colorlist[ind-3], label=label)
		ax[2,1].plot(wav_dict[elt], g_dict[elt],linewidth=1, color=colorlist[ind-3], label=label)

	for ind in range(6,7):
		label=label_list[ind]
		elt=file_list[ind]
		ax[0,2].set_title('Dust')
		ax[0,2].plot(wav_dict[elt], qsca_dict[elt]/w_0_dict[elt], linewidth=1, color='black', label=label)
		ax[1,2].plot(wav_dict[elt], w_0_dict[elt], linewidth=1, color='black', label=label)
		ax[2,2].plot(wav_dict[elt], g_dict[elt], linewidth=1, color='black', label=label)

	ax[0,0].set_ylabel('$Q_{ext}$', fontsize=12)
	ax[1,0].set_ylabel('$\omega_0$', fontsize=12)
	ax[2,0].set_ylabel('$g$', fontsize=12)
	
	ax[0,0].set_ylim([1., 2.5])
	ax[0,1].set_ylim([1., 2.5])
	ax[0,2].set_ylim([1., 2.5])

	ax[1,0].set_ylim([0.5, 1.01])
	ax[1,1].set_ylim([0.5, 1.01])
	ax[1,2].set_ylim([0.5, 1.01])

	ax[2,0].set_ylim([0.6, 1.01])
	ax[2,1].set_ylim([0.6, 1.01])
	ax[2,2].set_ylim([0.6, 1.01])

	ax[1,0].legend(loc='lower right', fontsize=8)
	ax[1,1].legend(loc='lower right', fontsize=8)
	ax[1,2].legend(loc='upper right', fontsize=8)

	ax[2,2].set_xlim([100, 500])
	ax[2,2].xaxis.set_ticks([100, 300, 500])
	ax[2,1].set_xlabel('Wavelength (nm)')
	plt.tight_layout()
	plt.subplots_adjust(wspace=0.2, hspace=0.1)
	plt.savefig('./Plots/methods_mieparams.pdf', orientation='portrait',papertype='letter', format='pdf')
	
	plt.show()


########################
###
########################
if plot_atmprofile_z:
	"""
	Plot some sample atmospheric profiles
	"""
	########################
	###Read in Models
	########################
	z_tp_dict={}#initialize dict to hold altitude from TP file (cm)
	p_dict={} #initialize dict to hold pressure (bar)
	t_dict={} #initialize dict to hold temperature (K)
	z_conc_dict={}#initialize dict to hold altitude from molar concentration file (cm)
	conc_co2_dict={} #initialize dict to hold molar concentrations of CO2
	conc_h2o_dict={} #initialize dict to hold molar concentrations of H2O

	file_list=np.array(['0.02bar_250K','0.2bar_250K','2bar_250K'])
	label_list=np.array(['pCO$_2=0.02$ bar, $T_0=250$K', 'pCO$_2=0.2$ bar, $T_0=250$K', 'pCO$_2=2$ bar, $T_0=250$K'])
	
	for elt in file_list:
		z_tp_dict[elt], t_dict[elt], p_dict[elt]=np.genfromtxt('./TPProfiles/colddrymars_'+elt+'_tpprofile.dat', skip_header=2, skip_footer=0,usecols=(0,1,3), unpack=True) #cm, K, bar
		
		z_conc_dict[elt], conc_co2_dict[elt], conc_h2o_dict[elt]=np.genfromtxt('./MolarConcentrations/colddrymars_'+elt+'_molarconcentrations.dat', skip_header=2, skip_footer=0,usecols=(0,2,3), unpack=True) #K, bar

	########################
	###Plot
	########################
	fig, ax=plt.subplots(3,figsize=(16.5*cm2inch,10.), sharex=False, sharey=True)
	colors=cm.rainbow(np.linspace(0,1,len(file_list)))
	for ind in range(0, len(file_list)):
		label=label_list[ind]
		elt=file_list[ind]
		ax[0].plot(p_dict[elt], z_tp_dict[elt]*1.e-5, linewidth=1, color=colors[ind], label=label)
		ax[1].plot(t_dict[elt], z_tp_dict[elt]*1.e-5, linewidth=1, color=colors[ind], label=label)
		ax[2].plot(conc_h2o_dict[elt], z_conc_dict[elt]*1.e-5, linewidth=1, color=colors[ind], label=label)
	#ax[0].set_title('CO$_2$-H$_2$O Atmosphere Vertical Profiles')
	ax[0].set_xlabel('Pressure (bar)', fontsize=12)
	ax[0].set_xscale('log')
	ax[0].set_xlim([2.e-5, 2.])
	ax[1].set_xlabel('Temperature (K)', fontsize=12)
	ax[1].set_xlim([165, 255])
	ax[2].set_xlabel('[H$_2$O]', fontsize=12)
	ax[2].set_xscale('log')
	ax[2].set_xlim([2.e-8, 2.])
	
	ax[0].set_ylabel('Altitude (km)', fontsize=12)
	ax[1].set_ylabel('Altitude (km)', fontsize=12)
	ax[2].set_ylabel('Altitude (km)', fontsize=12)
	ax[2].set_ylim([0., 64.])
	
	ax[0].legend(loc='best', fontsize=10)
	plt.tight_layout()
	plt.savefig('./Plots/methods_atmo_prof_z.pdf', orientation='portrait',papertype='letter', format='pdf')

	plt.show()

########################
###
########################
if plot_atmprofile_tp:
	"""
	Plot some sample atmospheric profiles
	"""
	########################
	###Read in Models
	########################
	P_0_list=np.array(['0.02', '0.2', '2'])
	T_0_list=np.array(['210', '250', '300'])
	
	p_dict={} #initialize dict to hold pressure (bar)
	t_dict={} #initialize dict to hold temperature (K)

	#file_list=np.array(['0.02bar_250K','0.2bar_250K','2bar_250K'])
	#label_list=np.array(['$P_0=0.02$ bar, $T_0=250$K', '$P_0=0.2$ bar, $T_0=250$K', '$P_0=2$ bar, $T_0=250$K'])
	
	for P_0 in P_0_list:
		for T_0 in T_0_list:
			label=P_0+'bar_'+T_0+'K'
			
			t_dict[label], p_dict[label]=np.genfromtxt('./TPProfiles/colddrymars_'+label+'_tpprofile.dat', skip_header=2, skip_footer=0,usecols=(1,3), unpack=True) #cm, K, bar

	########################
	###Plot
	########################	
	fig, ax=plt.subplots(3,figsize=(16.5*cm2inch,10.), sharex=True, sharey=False)
	colors=cm.rainbow(np.linspace(0,1,len(T_0_list)))
	
	for P_ind in range(0, len(P_0_list)):
		P_0=P_0_list[P_ind]
		ax[P_ind].set_title('pCO$_2$='+P_0+' (bar)')
		for T_ind in range(0, len(T_0_list)):
			T_0=T_0_list[T_ind]
			label=P_0+'bar_'+T_0+'K'

			ax[P_ind].plot(t_dict[label], p_dict[label], linewidth=1, color=colors[T_ind], label='$T_0=$'+T_0+'K')
		
		ax[P_ind].set_ylabel('Pressure (bar)', fontsize=12)
		ax[P_ind].set_yscale('log')
		ax[P_ind].set_ylim([np.min(p_dict[label]), np.float(P_0)])
		ax[P_ind].invert_yaxis()
	
	ax[len(P_0_list)-1].set_xlabel('Temperature (K)', fontsize=12)
	ax[len(P_0_list)-1].set_xlim([165, 305])
	ax[len(P_0_list)-1].legend(loc='best', fontsize=12)

	plt.tight_layout()
	plt.savefig('./Plots/methods_atmo_prof_tp.pdf', orientation='portrait',papertype='letter', format='pdf')
	plt.show()
