# -*- coding: iso-8859-1 -*-
"""
Purpose of this code is to plot the figures from the Results section of our paper. 
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
plot_tpdep=True #plot the dependence of surface fluence on clear-sky H2O-CO2 atmospheres as a function of temperature and pressure.
#plot_td_dep=False #plot the dependence of surface fluence on whether or not CO2 temperature dependence is included, for clear-sky H2O-CO2 atmospheres as a function of pressure.
plot_td_dep_lowpressure=True #plot the dependence of surface fluence on whether or not CO2 temperature dependence is included, for clear-sky H2O-CO2 atmospheres, just for low pressure (pCO2=2e-5 bar)

plot_clouds=True #plot the impact of varying levels of CO2 and H2O clouds in a CO2-H2O atmosphere.

plot_pso2_pco2=True #plot the impact of varying levels of SO2 in a CO2-H2O-SO2 atmosphere
plot_ph2s_pco2=True #plot the impact of varying levels of SO2 in a CO2-H2O-SO2 atmosphere
plot_pso2_clouds=True #plot the impact of varying levels of SO2 a CO2-H2O-SO2 atmosphere with varying levels of CO2 clouds.
plot_ph2s_clouds=True #plot the impact of varying levels of H2S a CO2-H2O-H2S atmosphere with varying levels of CO2 clouds.

plot_dust_pco2=True #plot the impact of varying levels of dust in a CO2-H2O atmosphere
plot_dust_clouds=True #plot the impact of varying levels of dust a CO2-H2O atmosphere with varying levels of CO2 clouds.



########################
###
########################
if plot_dust_clouds:
	"""
	The purpose of this script is to plot the UV surface fluence for Mars with a dusty, CO2-H2O atmosphere.
	pCO2=0.02 bar
	dust OD=0.1-10
	cloud OD=1-1000
	
	Conditions: T_0=250, A=desert, z=0, no TD XCs, DeltaScaling
	"""
	########################
	###Read in Models
	########################	
	###Read in TOA input
	toa_wavelengths, toa_intensity=np.genfromtxt('./Solar_Input/general_youngsun_mars_spectral_input.dat', skip_header=1, skip_footer=0,usecols=(2,3), unpack=True) #nm, erg/s/cm2/nm
	
	###Read in files programatically
	#cloudOD=1
	od1_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	od1_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	#cloudOD=10
	od10_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	od10_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	#cloudOD=100
	od100_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	od100_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	#cloudOD=1000
	od1000_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	od1000_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)

	elt_list=np.array(['dustod=0.1', 'dustod=1', 'dustod=10'])

	for elt in elt_list:
		od1_wav_dict[elt], od1_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/colddrymars_0.02bar_250K_z=0_A=desert_noTD_DS_'+elt+'_co2cloudod=1_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
		
		od10_wav_dict[elt], od10_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/colddrymars_0.02bar_250K_z=0_A=desert_noTD_DS_'+elt+'_co2cloudod=10_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)

		od100_wav_dict[elt], od100_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/colddrymars_0.02bar_250K_z=0_A=desert_noTD_DS_'+elt+'_co2cloudod=100_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
		
		od1000_wav_dict[elt], od1000_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/colddrymars_0.02bar_250K_z=0_A=desert_noTD_DS_'+elt+'_co2cloudod=1000_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)



	########################
	###Plot
	########################
	numelt=len(elt_list)
	fig, ax=plt.subplots(4, figsize=(16.5*cm2inch,10), sharex=True, sharey=True)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,numelt))
	conclist=np.array([r'$\tau_d=0.1$',r'$\tau_d=1$',r'$\tau_d=10$'])
	
	ax[0].set_title(r'$\tau_{cloud}=1$')
	ax[1].set_title(r'$\tau_{cloud}=10$')
	ax[2].set_title(r'$\tau_{cloud}=100$')
	ax[3].set_title(r'$\tau_{cloud}=1000$')
	
	ax[0].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[1].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[2].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[3].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')

	for ind in range(0, numelt):
		elt=elt_list[ind]
		
		ax[0].plot(od1_wav_dict[elt], od1_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[1].plot(od10_wav_dict[elt], od10_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[2].plot(od100_wav_dict[elt], od100_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[3].plot(od1000_wav_dict[elt], od1000_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])


	#ax[0].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	ax[1].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	#ax[2].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')

	ax[3].set_ylim([1.e-6, 1.e4])
	ax[3].set_yscale('log')
	ax[3].set_xlim([100, 500])
	ax[3].set_xlabel('Wavelength (nm)')
	
	ax[0].legend(bbox_to_anchor=[0, 1.13, 1., .152], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=12)
	plt.tight_layout(rect=(0,0,1,0.94))
	plt.subplots_adjust(wspace=0., hspace=0.2)
	plt.savefig('./Plots/results_dust_clouds.pdf', orientation='portrait',papertype='letter', format='pdf')
	
	plt.show()

########################
###
########################
if plot_dust_pco2:
	"""
	The purpose of this script is to plot the UV surface fluence for Mars with a dusty, CO2-H2O atmosphere.
	pCO2=0.02-2 bar
	dust OD=0.1-10 (100--> absorption optical epth of >20 -->just no flux anywhere)
	
	Conditions: T_0=250, A=desert, z=0, no TD XCs, DeltaScaling, no clouds
	"""
	########################
	###Read in Models
	########################
	###Read in TOA input
	toa_wavelengths, toa_intensity=np.genfromtxt('./Solar_Input/general_youngsun_mars_spectral_input.dat', skip_header=1, skip_footer=0,usecols=(2,3), unpack=True) #nm, erg/s/cm2/nm
	
	###Read in files programatically
	#pCO2=0.02 bar
	pco2_0p02_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	pco2_0p02_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)

	#pCO2=0.02 bar
	pco2_0p2_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	pco2_0p2_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	
	#pCO2=0.02 bar
	pco2_2_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	pco2_2_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)

	elt_list=np.array(['dustod=0.1', 'dustod=1', 'dustod=10'])


	for elt in elt_list:
		pco2_0p02_wav_dict[elt], pco2_0p02_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/colddrymars_0.02bar_250K_z=0_A=desert_noTD_DS_'+elt+'.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)

		pco2_0p2_wav_dict[elt], pco2_0p2_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/colddrymars_0.2bar_250K_z=0_A=desert_noTD_DS_'+elt+'.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)

		pco2_2_wav_dict[elt], pco2_2_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/colddrymars_2bar_250K_z=0_A=desert_noTD_DS_'+elt+'.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
	########################
	###Plot
	########################
	numelt=len(elt_list)
	fig, ax=plt.subplots(3, figsize=(16.5*cm2inch,10), sharex=True, sharey=True)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,numelt))
	conclist=np.array([r'$\tau_d=0.1$',r'$\tau_d=1$',r'$\tau_d=10$'])
	
	ax[0].set_title('pCO2=0.02 bar')
	ax[1].set_title('pCO2=0.2 bar')
	ax[2].set_title('pCO2=2 bar')
	
	ax[0].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[1].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[2].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')

	for ind in range(0, numelt):
		elt=elt_list[ind]

		
		ax[0].plot(pco2_0p02_wav_dict[elt], pco2_0p02_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[1].plot(pco2_0p2_wav_dict[elt], pco2_0p2_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[2].plot(pco2_2_wav_dict[elt], pco2_2_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])

	#ax[0].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	ax[1].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	#ax[2].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')

	ax[2].set_ylim([1.e-5, 1.e4])
	ax[2].set_yscale('log')
	ax[2].set_xlim([100, 500])
	ax[2].set_xlabel('Wavelength (nm)')
	
	ax[0].legend(bbox_to_anchor=[0, 1.13, 1., .152], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10)
	plt.tight_layout(rect=(0,0,1,0.94))
	plt.savefig('./Plots/results_dust_pco2.pdf', orientation='portrait',papertype='letter', format='pdf')
	
	plt.show()

########################
###
########################
if plot_ph2s_clouds:
	"""
	The purpose of this script is to plot the UV surface fluence for Mars with a CO2-H2O-SO2/H2S atmosphere, with CO2 cloud decks of varying thickness emplaced at 20.5 km (20-21 km). 
	pCO2=0.02 bar (optically thin in gas scattering)
	pH2S=2e-9 -- 2e-4 bar
	CO2 cloud OD varies from 1-1000.
	
	Conditions: T_0=250, A=desert, z=0, no TD XCs, DeltaScaling
	"""
	########################
	###Read in Models
	########################	
	###Read in TOA input
	toa_wavelengths, toa_intensity=np.genfromtxt('./Solar_Input/general_youngsun_mars_spectral_input.dat', skip_header=1, skip_footer=0,usecols=(2,3), unpack=True) #nm, erg/s/cm2/nm
	
	###Read in files programatically
	#cloudOD=1
	od1_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	od1_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	#cloudOD=10
	od10_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	od10_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	#cloudOD=100
	od100_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	od100_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	#cloudOD=1000
	od1000_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	od1000_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)

	elt_list=np.array(['0.02bar_250K_0so2_100ppbh2s','0.02bar_250K_0so2_1ppmh2s','0.02bar_250K_0so2_10ppmh2s', '0.02bar_250K_0so2_100ppmh2s', '0.02bar_250K_0so2_1000ppmh2s', '0.02bar_250K_0so2_10000ppmh2s'])

	for elt in elt_list:
		od1_wav_dict[elt], od1_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_DS_co2cloudod=1_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
		
		od10_wav_dict[elt], od10_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_DS_co2cloudod=10_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)

		od100_wav_dict[elt], od100_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_DS_co2cloudod=100_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
		
		od1000_wav_dict[elt], od1000_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_DS_co2cloudod=1000_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)



	########################
	###Plot
	########################
	numelt=len(elt_list)
	fig, ax=plt.subplots(4, figsize=(16.5*cm2inch,10), sharex=True, sharey=True)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,numelt))
	conclist=np.array([r'pH$_2$S$=2\times10^{-9}$ bar',r'pH$_2$S$=2\times10^{-8}$ bar',r'pH$_2$S$=2\times10^{-7}$ bar',r'pH$_2$S$=2\times10^{-6}$ bar',r'pH$_2$S$=2\times10^{-5}$ bar',r'pH$_2$S$=2\times10^{-4}$ bar'])
	
	ax[0].set_title(r'$\tau_{cloud}=1$')
	ax[1].set_title(r'$\tau_{cloud}=10$')
	ax[2].set_title(r'$\tau_{cloud}=100$')
	ax[3].set_title(r'$\tau_{cloud}=1000$')
	
	ax[0].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[1].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[2].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[3].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')

	for ind in range(0, numelt):
		elt=elt_list[ind]
		
		ax[0].plot(od1_wav_dict[elt], od1_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[1].plot(od10_wav_dict[elt], od10_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[2].plot(od100_wav_dict[elt], od100_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[3].plot(od1000_wav_dict[elt], od1000_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])


	#ax[0].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	ax[1].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	#ax[2].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')

	ax[3].set_ylim([1.e-6, 1.e4])
	ax[3].set_yscale('log')
	ax[3].set_xlim([100, 500])
	ax[3].set_xlabel('Wavelength (nm)')
	
	ax[0].legend(bbox_to_anchor=[0, 1.13, 1., .152], loc=3, ncol=3, mode='expand', borderaxespad=0., fontsize=10)
	plt.tight_layout(rect=(0,0,1,0.92))
	plt.subplots_adjust(wspace=0., hspace=0.2)
	plt.savefig('./Plots/results_ph2s_clouds.pdf', orientation='portrait',papertype='letter', format='pdf')
	
	plt.show()

########################
###
########################
if plot_pso2_clouds:
	"""
	The purpose of this script is to plot the UV surface fluence for Mars with a CO2-H2O-SO2/H2S atmosphere, with CO2 cloud decks of varying thickness emplaced at 20.5 km (20-21 km). 
	pCO2=0.02 bar (optically thin in gas scattering)
	pSO2=2e-9-2e-5 bar
	CO2 cloud OD varies from 1-1000.
	
	Conditions: T_0=250, A=desert, z=0, no TD XCs, DeltaScaling
	"""
	########################
	###Read in Models
	########################
	albedo='desert' #specify albedo condition
	
	###Read in TOA input
	toa_wavelengths, toa_intensity=np.genfromtxt('./Solar_Input/general_youngsun_mars_spectral_input.dat', skip_header=1, skip_footer=0,usecols=(2,3), unpack=True) #nm, erg/s/cm2/nm
	
	###Read in files programatically
	#cloudOD=1
	od1_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	od1_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	#cloudOD=10
	od10_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	od10_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	#cloudOD=100
	od100_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	od100_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	#cloudOD=1000
	od1000_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	od1000_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)

	elt_list=np.array(['0.02bar_250K_100ppbso2_0h2s','0.02bar_250K_1ppmso2_0h2s','0.02bar_250K_10ppmso2_0h2s', '0.02bar_250K_100ppmso2_0h2s', '0.02bar_250K_1000ppmso2_0h2s'])

	for elt in elt_list:
		od1_wav_dict[elt], od1_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_DS_co2cloudod=1_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
		
		od10_wav_dict[elt], od10_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_DS_co2cloudod=10_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)

		od100_wav_dict[elt], od100_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_DS_co2cloudod=100_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
		
		od1000_wav_dict[elt], od1000_surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_DS_co2cloudod=1000_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)



	########################
	###Plot
	########################
	numelt=len(elt_list)
	fig, ax=plt.subplots(4, figsize=(16.5*cm2inch,10), sharex=True, sharey=True)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,numelt))
	conclist=np.array([r'pSO$_2=2\times10^{-9}$ bar',r'pSO$_2=2\times10^{-8}$ bar',r'pSO$_2=2\times10^{-7}$ bar',r'pSO$_2=2\times10^{-6}$ bar',r'pSO$_2=2\times10^{-5}$ bar'])
	
	ax[0].set_title(r'$\tau_{cloud}=1$')
	ax[1].set_title(r'$\tau_{cloud}=10$')
	ax[2].set_title(r'$\tau_{cloud}=100$')
	ax[3].set_title(r'$\tau_{cloud}=1000$')
	
	ax[0].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[1].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[2].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[3].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')

	for ind in range(0, numelt):
		elt=elt_list[ind]
		
		ax[0].plot(od1_wav_dict[elt], od1_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[1].plot(od10_wav_dict[elt], od10_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[2].plot(od100_wav_dict[elt], od100_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[3].plot(od1000_wav_dict[elt], od1000_surfint_dict[elt], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])


	#ax[0].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	ax[1].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	#ax[2].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')

	ax[3].set_ylim([1.e-6, 1.e4])
	ax[3].set_yscale('log')
	ax[3].set_xlim([100, 500])
	ax[3].set_xlabel('Wavelength (nm)')
	
	ax[0].legend(bbox_to_anchor=[0, 1.13, 1., .152], loc=3, ncol=3, mode='expand', borderaxespad=0., fontsize=10)
	plt.tight_layout(rect=(0,0,1,0.94))
	plt.subplots_adjust(wspace=0., hspace=0.2)
	plt.savefig('./Plots/results_pso2_clouds.pdf', orientation='portrait',papertype='letter', format='pdf')
	
	plt.show()

########################
###
########################
if plot_ph2s_pco2:
	"""
	The purpose of this script is to plot the UV surface fluence for Mars with a CO2-H2O-H2S atmosphere. The H2S inventory varies. pCO2=0.02-2 bar is considered (to take into account CO2 multi-scattering effects). 
	pCO2=0.02-2 bar
	pH2S=2e-9 -- 2e-4 bar
	
	Conditions: T_0=250, A=desert, z=0, no TD XCs, no DeltaScaling, no clouds
	"""
	########################
	###Read in Models
	########################
	###Read in TOA input
	toa_wavelengths, toa_intensity=np.genfromtxt('./Solar_Input/general_youngsun_mars_spectral_input.dat', skip_header=1, skip_footer=0,usecols=(2,3), unpack=True) #nm, erg/s/cm2/nm
	
	###Read in files programatically
	surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	
	elt_list0=np.array(['0.02bar_250K_0so2_100ppbh2s','0.02bar_250K_0so2_1ppmh2s','0.02bar_250K_0so2_10ppmh2s', '0.02bar_250K_0so2_100ppmh2s', '0.02bar_250K_0so2_1000ppmh2s', '0.02bar_250K_0so2_10000ppmh2s'])
	elt_list1=np.array(['0.2bar_250K_0so2_10ppbh2s','0.2bar_250K_0so2_100ppbh2s','0.2bar_250K_0so2_1ppmh2s', '0.2bar_250K_0so2_10ppmh2s', '0.2bar_250K_0so2_100ppmh2s','0.2bar_250K_0so2_1000ppmh2s'])
	elt_list2=np.array(['2bar_250K_0so2_1ppbh2s', '2bar_250K_0so2_10ppbh2s','2bar_250K_0so2_100ppbh2s','2bar_250K_0so2_1ppmh2s', '2bar_250K_0so2_10ppmh2s', '2bar_250K_0so2_100ppmh2s'])

	for elt in elt_list0:
		wav_dict[elt], surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
	for elt in elt_list1:
		wav_dict[elt], surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
	for elt in elt_list2:
		wav_dict[elt], surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
	########################
	###Plot
	########################
	numelt=len(elt_list0)
	fig, ax=plt.subplots(3, figsize=(16.5*cm2inch,10), sharex=True, sharey=True)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,numelt))
	conclist=np.array([r'pH$_2$S$=2\times10^{-9}$ bar',r'pH$_2$S$=2\times10^{-8}$ bar',r'pH$_2$S$=2\times10^{-7}$ bar',r'pH$_2$S$=2\times10^{-6}$ bar',r'pH$_2$S$=2\times10^{-5}$ bar',r'pH$_2$S$=2\times10^{-4}$ bar'])
	
	ax[0].set_title('pCO2=0.02 bar')
	ax[1].set_title('pCO2=0.2 bar')
	ax[2].set_title('pCO2=2 bar')
	
	ax[0].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[1].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[2].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')

	for ind in range(0, numelt):
		elt0=elt_list0[ind]
		elt1=elt_list1[ind]
		elt2=elt_list2[ind]
		
		ax[0].plot(wav_dict[elt0], surfint_dict[elt0], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[1].plot(wav_dict[elt1], surfint_dict[elt1], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[2].plot(wav_dict[elt2], surfint_dict[elt2], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])

	#ax[0].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	ax[1].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	#ax[2].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')

	ax[2].set_ylim([1.e-5, 1.e4])
	ax[2].set_yscale('log')
	ax[2].set_xlim([100, 500])
	ax[2].set_xlabel('Wavelength (nm)')
	
	ax[0].legend(bbox_to_anchor=[0, 1.13, 1., .152], loc=3, ncol=3, mode='expand', borderaxespad=0., fontsize=10)
	plt.tight_layout(rect=(0,0,1,0.90))
	plt.savefig('./Plots/results_ph2s_pco2.pdf', orientation='portrait',papertype='letter', format='pdf')
	
	plt.show()

########################
###
########################
if plot_pso2_pco2:
	"""
	The purpose of this script is to plot the UV surface fluence for Mars with a CO2-H2O-SO2 atmosphere. The SO2 inventory varies. pCO2=0.02-2 bar is considered (to take into account CO2 multi-scattering effects). 
	pCO2=0.02-2 bar
	pSO2=2e-9 -- 2e-5 bar
	
	Conditions: T_0=250, A=desert, z=0, no TD XCs, no DeltaScaling, no clouds
	"""
	########################
	###Read in Models
	########################
	###Read in TOA input
	toa_wavelengths, toa_intensity=np.genfromtxt('./Solar_Input/general_youngsun_mars_spectral_input.dat', skip_header=1, skip_footer=0,usecols=(2,3), unpack=True) #nm, erg/s/cm2/nm
	
	###Read in files programatically
	surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	
	elt_list0=np.array(['0.02bar_250K_100ppbso2_0h2s','0.02bar_250K_1ppmso2_0h2s','0.02bar_250K_10ppmso2_0h2s', '0.02bar_250K_100ppmso2_0h2s', '0.02bar_250K_1000ppmso2_0h2s'])
	elt_list1=np.array(['0.2bar_250K_10ppbso2_0h2s','0.2bar_250K_100ppbso2_0h2s','0.2bar_250K_1ppmso2_0h2s', '0.2bar_250K_10ppmso2_0h2s', '0.2bar_250K_100ppmso2_0h2s'])
	elt_list2=np.array(['2bar_250K_1ppbso2_0h2s','2bar_250K_10ppbso2_0h2s','2bar_250K_100ppbso2_0h2s', '2bar_250K_1ppmso2_0h2s', '2bar_250K_10ppmso2_0h2s'])

	for elt in elt_list0:
		wav_dict[elt], surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
	for elt in elt_list1:
		wav_dict[elt], surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
	for elt in elt_list2:
		wav_dict[elt], surfint_dict[elt]=np.genfromtxt('./TwoStreamOutput/volcanicmars_'+elt+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
	########################
	###Plot
	########################
	numelt=len(elt_list0)
	fig, ax=plt.subplots(3, figsize=(16.5*cm2inch,10), sharex=True, sharey=True)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,numelt))
	conclist=np.array([r'pSO$_2=2\times10^{-9}$ bar',r'pSO$_2=2\times10^{-8}$ bar',r'pSO$_2=2\times10^{-7}$ bar',r'pSO$_2=2\times10^{-6}$ bar',r'pSO$_2=2\times10^{-5}$ bar'])
	
	ax[0].set_title('pCO2=0.02 bar')
	ax[1].set_title('pCO2=0.2 bar')
	ax[2].set_title('pCO2=2 bar')
	
	ax[0].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[1].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[2].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')

	for ind in range(0, numelt):
		elt0=elt_list0[ind]
		elt1=elt_list1[ind]
		elt2=elt_list2[ind]
		
		ax[0].plot(wav_dict[elt0], surfint_dict[elt0], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[1].plot(wav_dict[elt1], surfint_dict[elt1], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])
		ax[2].plot(wav_dict[elt2], surfint_dict[elt2], marker='s', markersize=markersizeval, linewidth=1, color=colors[ind], label=conclist[ind])

	#ax[0].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	ax[1].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	#ax[2].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')

	ax[2].set_ylim([1.e-5, 1.e4])
	ax[2].set_yscale('log')
	ax[2].set_xlim([100, 500])
	ax[2].set_xlabel('Wavelength (nm)')
	
	ax[0].legend(bbox_to_anchor=[0, 1.13, 1., .152], loc=3, ncol=3, mode='expand', borderaxespad=0., fontsize=10)
	plt.tight_layout(rect=(0,0,1,0.94))
	plt.savefig('./Plots/results_pso2_pco2.pdf', orientation='portrait',papertype='letter', format='pdf')
	
	plt.show()


########################
###
########################
if plot_clouds:
	"""
	The purpose of this script is to plot the UV surface fluence on a cold, dry (CO2+H2O) Martian atmosphere that is clear, except for I) H2O cloud deck of varying optical thickness (at 500 nm) at 3.5 km (3-4 km), and II) CO2 cloud deck of varying optical thickness (at 500 nm) at 20.5 (20-21) km. 
	
	Conditions: pCO2=0.02 bar, T_0=250, A=desert, z=0, no TD XCs, yes DeltaScaling
	"""
	########################
	###Read in Models
	########################
	###Read in TOA input
	toa_wavelengths, toa_intensity=np.genfromtxt('./Solar_Input/general_youngsun_mars_spectral_input.dat', skip_header=1, skip_footer=0,usecols=(2,3), unpack=True) #nm, erg/s/cm2/nm
	
	###Read in computed surface fluences
	cloudtau_list=np.array(['0.1', '1', '10', '100', '1000', '10000']) #list of cloud optical depths (500 nm)
	
	#no temperature dependence
	co2_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm) for co2 cloud decks
	co2_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm) for co2 cloud decks
	
	h2o_surfint_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm) for h2o cloud decks
	h2o_wav_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm) for h2o cloud decks	

	for cloudtau in cloudtau_list:
		h2o_wav_dict[cloudtau], h2o_surfint_dict[cloudtau]=np.genfromtxt('./TwoStreamOutput/colddrymars_0.02bar_250K_z=0_A=desert_noTD_DS_h2ocloudod='+cloudtau+'_z=3.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
		
		co2_wav_dict[cloudtau], co2_surfint_dict[cloudtau]=np.genfromtxt('./TwoStreamOutput/colddrymars_0.02bar_250K_z=0_A=desert_noTD_DS_co2cloudod='+cloudtau+'_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
	
	########################
	###Plot
	########################
	numtau=len(cloudtau_list)
	fig, ax=plt.subplots(2, figsize=(cm2inch*16.5,11.), sharex=True, sharey=False)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,numtau))
	
	ax[0].set_title('CO$_2$ Cloud Deck')
	ax[1].set_title('H$_2$O Cloud Deck')
	ax[0].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	ax[1].plot(toa_wavelengths, toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
	
	for tau_ind in range(0, numtau):
		cloudtau=cloudtau_list[tau_ind]

		ax[0].plot(co2_wav_dict[cloudtau], co2_surfint_dict[cloudtau], marker='s', markersize=markersizeval, linewidth=1, color=colors[tau_ind], label='Cloud OD='+cloudtau)

		ax[1].plot(h2o_wav_dict[cloudtau], h2o_surfint_dict[cloudtau], marker='s', markersize=markersizeval, linewidth=1, color=colors[tau_ind], label='Cloud OD='+cloudtau)
		

	ax[0].set_ylim([1.e-2, 1.e4])
	ax[0].set_yscale('log')
	ax[0].legend(loc=0, fontsize=10)
	ax[0].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	
	ax[1].set_ylim([1.e-2, 1.e4])
	ax[1].set_yscale('log')
	ax[1].legend(loc=0, fontsize=10)
	ax[1].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	ax[1].set_xlabel('Wavelength (nm)')

	plt.tight_layout()
	plt.savefig('./Plots/results_clouds.pdf', orientation='portrait',papertype='letter', format='pdf')
	
	plt.show()

#########################
####
#########################
#if plot_td_dep:
	#"""
	#The purpose of this script is to plot the UV surface fluence on a cold, dry (CO2+H2O) Martian atmosphere for different surface pressures, with and without temperature dependence included.
	
	#Surface Pressures (CO2 dry partial pressure):2e-5, 2e-4, 2e-3 bar 
	#Surface Temperatures: 210K
	#A=desert, z=0. NOTE: Verify this -- might be A=0.2
	#"""
	#########################
	####Read in Models
	#########################
	####Read in TOA input
	#toa_wavelengths, toa_intensity=np.genfromtxt('./Solar_Input/general_youngsun_mars_spectral_input.dat', skip_header=1, skip_footer=0,usecols=(2,3), unpack=True) #nm, erg/s/cm2/nm
	
	####Read in computed surface fluences
	#P_0_list=np.array(['0.00002', '0.0002', '0.002'])
	#T_0='200'
	
	##no temperature dependence
	#surfint_notd_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	#wav_notd_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	##temperature dependence included
	#surfint_td_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	#wav_td_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	
	#for P_0 in P_0_list:
		#label=P_0+'bar_'+T_0+'K'
		
		#wav_notd_dict[label], surfint_notd_dict[label]=np.genfromtxt('./TwoStreamOutput/colddrymars_'+label+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
		
		#wav_td_dict[label], surfint_td_dict[label]=np.genfromtxt('./TwoStreamOutput/colddrymars_'+label+'_z=0_A=desert_TD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)

	
	#########################
	####Plot results
	#########################
	#numP=len(P_0_list)
	
	#fig, ax=plt.subplots(numP, figsize=(16.5*cm2inch,10.), sharex=True, sharey=True)
	#markersizeval=5.
	
	#for P_ind in range(0, numP):
		#P_0=P_0_list[P_ind]
		#ax[P_ind].set_title('pCO$_2$='+P_0+'bar')
		#ax[P_ind].plot(toa_wavelengths,toa_intensity,marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
			
		#label=P_0+'bar_'+T_0+'K'
		
		#ax[P_ind].plot(wav_notd_dict[label], surfint_notd_dict[label], marker='s', markersize=markersizeval, linewidth=1, color='red', label='No Temp Dependence')
		
		#ax[P_ind].plot(wav_td_dict[label], surfint_td_dict[label], marker='s', markersize=markersizeval, linewidth=1, color='blue', label='Temp Dependence')

		
		#ax[P_ind].set_yscale('log')
		#ax[P_ind].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	
	#ax[numP-1].set_ylim([1.e-2, 1.e4])
	#ax[numP-1].set_xlim([100, 500])
	#ax[numP-1].set_xlabel('Wavelength (nm)')
	
	
	####Plot differences
	#ax[0].legend(loc=0, fontsize=12)
	#plt.tight_layout()
	#plt.savefig('./Plots/results_td.pdf', orientation='portrait',papertype='letter', format='pdf')
	
	#plt.show
	
	#fig2, ax=plt.subplots(numP, figsize=(16.5*cm2inch,10.), sharex=True, sharey=True)
	#markersizeval=5.
	
	#for P_ind in range(0, numP):
		#P_0=P_0_list[P_ind]
		#ax[P_ind].set_title(P_0+'bar')
		
		#label=P_0+'bar_'+T_0+'K'
		
		#ax[P_ind].plot(wav_notd_dict[label], (surfint_notd_dict[label]-surfint_td_dict[label])/toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black')

			
		#ax[P_ind].set_yscale('linear')
	

	#ax[numP-1].set_ylim([-1., 1.e-6])
	#ax[numP-1].set_xlim([100, 500])
	#ax[numP-1].set_xlabel('Wavelength (nm)')
	#ax[1].set_ylabel('(TD-noTD)/TOA (Fractional Difference)')

	#plt.tight_layout()
	#plt.show()

########################
###
########################
if plot_td_dep_lowpressure:
	"""
	The purpose of this script is to plot the UV surface fluence on a cold, dry (CO2+H2O) Martian atmosphere for a low surface pressure, with and without temperature dependence included.
	
	Surface Pressures (CO2 dry partial pressure):2e-5
	Surface Temperatures: 200K
	A=desert, z=0. NOTE: Verify this -- might be A=0.2
	"""
	########################
	###Read in Models
	########################
	###Read in TOA input
	toa_wavelengths, toa_intensity=np.genfromtxt('./Solar_Input/general_youngsun_mars_spectral_input.dat', skip_header=1, skip_footer=0,usecols=(2,3), unpack=True) #nm, erg/s/cm2/nm
	
	###Read in computed surface fluences
	P_0_list=np.array(['0.00002'])
	T_0='200'
	
	#no temperature dependence
	surfint_notd_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	wav_notd_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	#temperature dependence included
	surfint_td_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	wav_td_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	
	for P_0 in P_0_list:
		label=P_0+'bar_'+T_0+'K'
		
		wav_notd_dict[label], surfint_notd_dict[label]=np.genfromtxt('./TwoStreamOutput/colddrymars_'+label+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)
		
		wav_td_dict[label], surfint_td_dict[label]=np.genfromtxt('./TwoStreamOutput/colddrymars_'+label+'_z=0_A=desert_TD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)

	
	########################
	###Plot results
	########################
	numP=len(P_0_list)
	
	fig, ax=plt.subplots(numP, figsize=(16.5*cm2inch,4.), sharex=True, sharey=True)
	markersizeval=5.
	
	for P_ind in range(0, numP):
		P_0=P_0_list[P_ind]
		#ax.set_title('pCO$_2$='+P_0+'bar')
		ax.plot(toa_wavelengths,toa_intensity,marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
			
		label=P_0+'bar_'+T_0+'K'
		
		ax.plot(wav_notd_dict[label], surfint_notd_dict[label], marker='s', markersize=markersizeval, linewidth=1, color='red', label='No Temp Dependence')
		
		ax.plot(wav_td_dict[label], surfint_td_dict[label], marker='s', markersize=markersizeval, linewidth=1, color='blue', label='Temp Dependence')

		
		ax.set_yscale('log')
		ax.set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	
	ax.set_ylim([1.e-2, 1.e4])
	ax.set_xlim([100, 500])
	ax.set_xlabel('Wavelength (nm)')
	
	
	ax.legend(loc=0, fontsize=12)
	plt.tight_layout()
	plt.savefig('./Plots/results_td_lowpressure.pdf', orientation='portrait',papertype='letter', format='pdf')
	
	plt.show()
	###Plot differences

	fig2, ax=plt.subplots(numP, figsize=(16.5*cm2inch,10.), sharex=True, sharey=True)
	markersizeval=5.
	
	for P_ind in range(0, numP):
		P_0=P_0_list[P_ind]
		ax.set_title(P_0+'bar')
		
		label=P_0+'bar_'+T_0+'K'
		
		ax.plot(wav_notd_dict[label], (surfint_notd_dict[label]-surfint_td_dict[label])/toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color='black')

			
		ax.set_yscale('linear')
	

	ax.set_ylim([-1., 1.e-6])
	ax.set_xlim([100, 500])
	ax.set_xlabel('Wavelength (nm)')
	ax.set_ylabel('(TD-noTD)/TOA (Fractional Difference)')

	plt.tight_layout()
	plt.show()

########################
###
########################
if plot_tpdep:
	"""
	The purpose of this script is to plot the UV surface fluence on a cold, dry (CO2+H2O) Martian atmosphere for different surface temperatures and pressures. Temperature dependence is NOT included. 
	
	Surface Pressures (CO2 dry partial pressure): 0.02 bar, 0.2 bar, 2 bar
	Surface Temperatures: 210K, 250K, 300K.
	
	No particulates, no TDXC, no delta scaling, SZA=0, A=desert
	"""
	########################
	###Read in Models
	########################
	###Read in TOA input
	toa_wavelengths, toa_intensity=np.genfromtxt('./Solar_Input/general_youngsun_mars_spectral_input.dat', skip_header=1, skip_footer=0,usecols=(2,3), unpack=True) #nm, erg/s/cm2/nm
	
	###Read in computed surface fluences
	P_0_list=np.array(['0.02', '0.2', '2'])
	T_0_list=np.array(['210', '250', '300'])
	surface_intensity_dict={} #initialize dict to hold the surface intensities (ergs/cm2/s/nm)
	wavelengths_dict={}#initialize dict to hold corresponding wavelengths (centers of wavelength bins, in nm)
	
	for P_0 in P_0_list:
		for T_0 in T_0_list:
			label=P_0+'bar_'+T_0+'K'
			
			wavelengths_dict[label], surface_intensity_dict[label]=np.genfromtxt('./TwoStreamOutput/colddrymars_'+label+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(2,6), unpack=True)

	
	########################
	###Plot results
	########################
	numP=len(P_0_list)
	numT=len(T_0_list)
	
	fig, ax=plt.subplots(numP, figsize=(16.5*cm2inch,10.), sharex=True, sharey=True)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,numT))
	
	for P_ind in range(0, numP):
		P_0=P_0_list[P_ind]
		ax[P_ind].set_title('pCO$_2$='+P_0+'bar')
		ax[P_ind].plot(toa_wavelengths,toa_intensity,marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA')
		
		for T_ind in range(0, numT):
			T_0=T_0_list[T_ind]
			
			label=P_0+'bar_'+T_0+'K'
			
			ax[P_ind].plot(wavelengths_dict[label], surface_intensity_dict[label], marker='s', markersize=markersizeval, linewidth=1, color=colors[T_ind], label=T_0+'K')
			
		ax[P_ind].set_yscale('log')
		ax[P_ind].set_ylabel('Integrated Intensity (erg/s/cm2/nm)')
	
	ax[numP-1].set_ylim([1.e-2, 1.e4])
	ax[numP-1].set_xlim([100, 500])
	ax[numP-1].set_xlabel('Wavelength (nm)')
	
	
	ax[0].legend(loc=0, fontsize=12)
	plt.tight_layout()
	plt.savefig('./Plots/results_tpdep.pdf', orientation='portrait',papertype='letter', format='pdf')

	#PLOT DIFFERENCES
	fig2, ax=plt.subplots(numP, figsize=(16.5*cm2inch,10.), sharex=True, sharey=True)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,numT))
	
	for P_ind in range(0, numP):
		P_0=P_0_list[P_ind]
		ax[P_ind].set_title(P_0+'bar')
		
		for T_ind in range(1, numT):
			T_0=T_0_list[T_ind]
			
			label=P_0+'bar_'+T_0+'K'
			label0=P_0+'bar_'+T_0_list[0]+'K'
			
			ax[P_ind].plot(wavelengths_dict[label], np.abs(surface_intensity_dict[label]-surface_intensity_dict[label0])/toa_intensity, marker='s', markersize=markersizeval, linewidth=1, color=colors[T_ind], label=T_0+'K')
			
		ax[P_ind].set_yscale('log')
	

	ax[numP-1].set_ylim([1.e-6, 1.e1])
	ax[numP-1].set_xlim([100, 500])
	ax[numP-1].set_xlabel('Wavelength (nm)')
	ax[1].set_ylabel('abs([Int(T)-Int('+T_0_list[0]+'K)])/TOA (Fractional Difference)')

	ax[0].legend(loc=0, fontsize=12)
	plt.tight_layout()
	plt.show()
