# -*- coding: iso-8859-1 -*-
"""
Purpose of this code is to plot the figures from the Discussion section of our paper. 
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
plot_doses_pco2=True #plot the dependence of dose rate on pCO2 in a CO2-H2O atmosphere (fixed temperature) NOTE: may want to do for low-albedo. Both more physically plausible and avoids weird uptick effect

plot_doses_clouds=True #plot the dependence of dose rate on CO2 cloud optical depth in a CO2-H2O atmosphere (fixed temperature)

plot_doses_dust_pco2=True #plot the dependence of dose rate as a function of dust level for different pCO2
plot_doses_dust_clouds=True #plot the dependence of dose rate as a function of dust level for different cloud levels.

plot_doses_pso2_pco2=True #plot the dependence of dose rate as a function of pSO2 for different pCO2
plot_doses_pso2_clouds=True #plot the dependence of dose rate as a function of dust level for different cloud levels.

plot_doses_ph2s_pco2=True #plot the dependence of dose rate as a function of pH2S for different pCO2
plot_doses_ph2s_clouds=True #plot the dependence of dose rate as a function of pH2S for different pCO2

plot_reldoses_pso2=True #plot the ratio between the "bad" dose rate and the "good" dose rate as function of PSO2. Plot for 1) pCO2=0.02, cloud=1000 and 2) pCO2=2 bar, cloud=0. 
plot_reldoses_ph2s=True #plot the ratio between the "bad" dose rate and the "good" dose rate as function of PH2S. Plot for 1) pCO2=0.02, cloud=1000 and 2) pCO2=2 bar, cloud=0. 
plot_reldoses_dust=True #plot the ratio between the "bad" dose rate and the "good" dose rate as function of dust level. Plot for 1) pCO2=0.02, cloud=1000 and 2) pCO2=2 bar, cloud=0. 


########################
###
########################
if plot_reldoses_dust:
	"""
	The purpose of this script is to plot the relative dose rates of the "stressor" photoprocess compared to the "eustressor" photoprocess as a function of pH2S. We evaluate this for varying levels of atmospheric dust loading. We do this for two cases: 1) pCO2=2 bar, and 2) pCO2=0.02 bar and a tau=1000 cloud deck emplaced at 20.5 km. This is so we can separate out absorption amplification due to cloud decks (relatively flat) and due to Rayleigh scattering (not flat).
	
	#Conditions: T_0=250, A=desert, z=0, no TD XCs, DeltaScaling
	"""
	########################
	###Read in doses
	########################
	###Read in computed doses
	od_list=np.array(['dustod=0.1', 'dustod=1', 'dustod=10']) #list of dust optical depths (500 nm)
	od_axis=np.array([1.e-1, 1., 1.e1])	

	titles_list=np.array([r'pCO$_2$=2 bar, $\tau_{cloud}=0$',r'pCO$_2$=0.02 bar, $\tau_{cloud}=1000$ (unscaled)'])
	num_cases=len(titles_list)
	
	
	num_od=len(od_list)
	
	dose_100_165=np.zeros([num_od,num_cases]) #surface radiance integrated 100-165 nm
	dose_200_300=np.zeros([num_od,num_cases]) #surface radiance integrated 200-300 nm
	dose_ump_193=np.zeros([num_od,num_cases]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=193
	dose_ump_230=np.zeros([num_od,num_cases]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=230
	dose_ump_254=np.zeros([num_od,num_cases]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=254
	dose_cucn3_254=np.zeros([num_od,num_cases]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=254
	dose_cucn3_300=np.zeros([num_od,num_cases]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=300
	
	for ind in range(0, num_od):
		od=od_list[ind]
		dose_100_165[ind,0],dose_200_300[ind,0],dose_ump_193[ind,0],dose_ump_230[ind,0],dose_ump_254[ind,0],dose_cucn3_254[ind,0],dose_cucn3_300[ind,0]=np.genfromtxt('./DoseRates/dose_rates_colddrymars_2bar_250K_z=0_A=desert_noTD_DS_'+od+'.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)
		
		dose_100_165[ind,1],dose_200_300[ind,1],dose_ump_193[ind,1],dose_ump_230[ind,1],dose_ump_254[ind,1],dose_cucn3_254[ind,1],dose_cucn3_300[ind,1]=np.genfromtxt('./DoseRates/dose_rates_colddrymars_0.02bar_250K_z=0_A=desert_noTD_DS_'+od+'_co2cloudod=1000_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)

		
	########################
	###Plot results
	########################
	fig, ax=plt.subplots(num_cases, figsize=(16.5*cm2inch,10.), sharex=True, sharey=False)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,6))
	
	for ind2 in range(0, num_cases):
		ax[ind2].set_title(titles_list[ind2])
		
		ax[ind2].plot(od_axis, dose_ump_193[:,ind2]/dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[0], label=r'UMP-193/CuCN3-254')
		ax[ind2].plot(od_axis, dose_ump_230[:,ind2]/dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[1], label=r'UMP-230/CuCN3-254')
		ax[ind2].plot(od_axis, dose_ump_254[:,ind2]/dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[2], label=r'UMP-254/CuCN3-254')

		ax[ind2].plot(od_axis, dose_ump_193[:,ind2]/dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[3], label=r'UMP-193/CuCN3-300')
		ax[ind2].plot(od_axis, dose_ump_230[:,ind2]/dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[4], label=r'UMP-230/CuCN3-300')
		ax[ind2].plot(od_axis, dose_ump_254[:,ind2]/dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[5], label=r'UMP-254/CuCN3-300')
		
		#ax.set_ylim([1.e-2, 1.e4])
		ax[ind2].set_yscale('log')
		ax[ind2].set_ylabel(r'$\bar{D}_{UMP-X}/\bar{D}_{CuCN3-Y}$')
	#ax.set_xlim([100, 500])
	ax[num_cases-1].set_xscale('log')
	ax[num_cases-1].set_xlabel(r'$\tau_{d}$ (unscaled)', fontsize=12)
	
	
	plt.tight_layout(rect=(0,0,1., 0.9))
	ax[0].legend(bbox_to_anchor=[0, 1.2, 1., 0.7], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10) 
	plt.savefig('./Plots/discussion_doses_reldoses_dust.pdf', orientation='portrait',papertype='letter', format='pdf')

	plt.show()

########################
###
########################
if plot_reldoses_ph2s:
	"""
	The purpose of this script is to plot the relative dose rates of the "stressor" photoprocess compared to the "eustressor" photoprocess as a function of pH2S. We evaluate this for pSO2=2e-9 -- 2e-4 bar. We do this for two cases: 1) pCO2=2 bar, and 2) pCO2=0.02 bar and a tau=1000 cloud deck emplaced at 20.5 km. This is so we can separate out absorption amplification due to cloud decks (relatively flat) and due to Rayleigh scattering (not flat).
	
	#Conditions: T_0=250, A=desert, z=0, no TD XCs, DeltaScaling
	"""
	########################
	###Read in doses
	########################
	###Read in computed doses
	nocloud_list=np.array(['2bar_250K_0so2_1ppbh2s', '2bar_250K_0so2_10ppbh2s','2bar_250K_0so2_100ppbh2s','2bar_250K_0so2_1ppmh2s', '2bar_250K_0so2_10ppmh2s', '2bar_250K_0so2_100ppmh2s']) #list of pH2S for pCO2=2
	cloud_list=np.array(['0.02bar_250K_0so2_100ppbh2s','0.02bar_250K_0so2_1ppmh2s','0.02bar_250K_0so2_10ppmh2s', '0.02bar_250K_0so2_100ppmh2s', '0.02bar_250K_0so2_1000ppmh2s', '0.02bar_250K_0so2_10000ppmh2s']) #list of pH2S for pCO2=0.02
	
	ph2s_axis=np.array([2.e-9, 2.e-8, 2.e-7, 2.e-6, 2.e-5, 2.e-4]) #pH2S in bar
	titles_list=np.array([r'pCO$_2$=2 bar, $\tau_{cloud}=0$',r'pCO$_2$=0.02 bar, $\tau_{cloud}=1000$ (unscaled)'])
	num_cases=len(titles_list)
	
	
	num_h2s=len(nocloud_list)
	
	dose_100_165=np.zeros([num_h2s,num_cases]) #surface radiance integrated 100-165 nm
	dose_200_300=np.zeros([num_h2s,num_cases]) #surface radiance integrated 200-300 nm
	dose_ump_193=np.zeros([num_h2s,num_cases]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=193
	dose_ump_230=np.zeros([num_h2s,num_cases]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=230
	dose_ump_254=np.zeros([num_h2s,num_cases]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=254
	dose_cucn3_254=np.zeros([num_h2s,num_cases]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=254
	dose_cucn3_300=np.zeros([num_h2s,num_cases]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=300
	
	for ind in range(0, num_h2s):
		dose_100_165[ind,0],dose_200_300[ind,0],dose_ump_193[ind,0],dose_ump_230[ind,0],dose_ump_254[ind,0],dose_cucn3_254[ind,0],dose_cucn3_300[ind,0]=np.genfromtxt('./DoseRates/dose_rates_volcanicmars_'+nocloud_list[ind]+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)
		
		dose_100_165[ind,1],dose_200_300[ind,1],dose_ump_193[ind,1],dose_ump_230[ind,1],dose_ump_254[ind,1],dose_cucn3_254[ind,1],dose_cucn3_300[ind,1]=np.genfromtxt('./DoseRates/dose_rates_volcanicmars_'+cloud_list[ind]+'_z=0_A=desert_noTD_DS_co2cloudod=1000_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)

		
	########################
	###Plot results
	########################
	fig, ax=plt.subplots(num_cases, figsize=(16.5*cm2inch,10.), sharex=True, sharey=False)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,6))
	
	for ind2 in range(0, num_cases):
		ax[ind2].set_title(titles_list[ind2])
		
		ax[ind2].plot(ph2s_axis, dose_ump_193[:,ind2]/dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[0], label=r'UMP-193/CuCN3-254')
		ax[ind2].plot(ph2s_axis, dose_ump_230[:,ind2]/dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[1], label=r'UMP-230/CuCN3-254')
		ax[ind2].plot(ph2s_axis, dose_ump_254[:,ind2]/dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[2], label=r'UMP-254/CuCN3-254')

		ax[ind2].plot(ph2s_axis, dose_ump_193[:,ind2]/dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[3], label=r'UMP-193/CuCN3-300')
		ax[ind2].plot(ph2s_axis, dose_ump_230[:,ind2]/dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[4], label=r'UMP-230/CuCN3-300')
		ax[ind2].plot(ph2s_axis, dose_ump_254[:,ind2]/dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[5], label=r'UMP-254/CuCN3-300')
		
		#ax.set_ylim([1.e-2, 1.e4])
		ax[ind2].set_yscale('log')
		ax[ind2].set_ylabel(r'$\bar{D}_{UMP-X}/\bar{D}_{CuCN3-Y}$')
	#ax.set_xlim([100, 500])
	ax[num_cases-1].set_xscale('log')
	ax[num_cases-1].set_xlabel(r'pH$_2$S', fontsize=12)
	
	
	plt.tight_layout(rect=(0,0,1., 0.9))
	ax[0].legend(bbox_to_anchor=[0, 1.2, 1., 0.7], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10) 
	plt.savefig('./Plots/discussion_doses_reldoses_ph2s.pdf', orientation='portrait',papertype='letter', format='pdf')

	plt.show()

########################
###
########################
if plot_reldoses_pso2:
	"""
	The purpose of this script is to plot the relative dose rates of the "stressor" photoprocess compared to the "eustressor" photoprocess as a function of pSO2. We evaluate this for pSO2=2e-9 -- 2e-5 bar. We do this for two cases: 1) pCO2=2 bar, and 2) pCO2=0.02 bar and a tau=1000 cloud deck emplaced at 20.5 km. This is so we can separate out absorption amplification due to cloud decks (relatively flat) and due to Rayleigh scattering (not flat).
	
	#Conditions: T_0=250, A=desert, z=0, no TD XCs, DeltaScaling
	"""
	########################
	###Read in doses
	########################
	###Read in computed doses
	nocloud_list=np.array(['2bar_250K_1ppbso2_0h2s','2bar_250K_10ppbso2_0h2s','2bar_250K_100ppbso2_0h2s', '2bar_250K_1ppmso2_0h2s', '2bar_250K_10ppmso2_0h2s']) #list of pSO2 for pCO2=2
	cloud_list=np.array(['0.02bar_250K_100ppbso2_0h2s','0.02bar_250K_1ppmso2_0h2s','0.02bar_250K_10ppmso2_0h2s', '0.02bar_250K_100ppmso2_0h2s', '0.02bar_250K_1000ppmso2_0h2s']) #list of pSO2 for pCO2=0.02 
	
	pso2_axis=np.array([2.e-9, 2.e-8, 2.e-7, 2.e-6, 2.e-5]) #pSO2 in bar
	titles_list=np.array([r'pCO$_2$=2 bar, $\tau_{cloud}=0$',r'pCO$_2$=0.02 bar, $\tau_{cloud}=1000$ (unscaled)'])
	num_cases=len(titles_list)
	
	
	num_so2=len(nocloud_list)
	
	dose_100_165=np.zeros([num_so2,num_cases]) #surface radiance integrated 100-165 nm
	dose_200_300=np.zeros([num_so2,num_cases]) #surface radiance integrated 200-300 nm
	dose_ump_193=np.zeros([num_so2,num_cases]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=193
	dose_ump_230=np.zeros([num_so2,num_cases]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=230
	dose_ump_254=np.zeros([num_so2,num_cases]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=254
	dose_cucn3_254=np.zeros([num_so2,num_cases]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=254
	dose_cucn3_300=np.zeros([num_so2,num_cases]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=300
	
	for ind in range(0, num_so2):
		dose_100_165[ind,0],dose_200_300[ind,0],dose_ump_193[ind,0],dose_ump_230[ind,0],dose_ump_254[ind,0],dose_cucn3_254[ind,0],dose_cucn3_300[ind,0]=np.genfromtxt('./DoseRates/dose_rates_volcanicmars_'+nocloud_list[ind]+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)
		
		dose_100_165[ind,1],dose_200_300[ind,1],dose_ump_193[ind,1],dose_ump_230[ind,1],dose_ump_254[ind,1],dose_cucn3_254[ind,1],dose_cucn3_300[ind,1]=np.genfromtxt('./DoseRates/dose_rates_volcanicmars_'+cloud_list[ind]+'_z=0_A=desert_noTD_DS_co2cloudod=1000_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)

		
	########################
	###Plot results
	########################
	fig, ax=plt.subplots(num_cases, figsize=(16.5*cm2inch,10.), sharex=True, sharey=False)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,6))
	
	for ind2 in range(0, num_cases):
		ax[ind2].set_title(titles_list[ind2])
		
		ax[ind2].plot(pso2_axis, dose_ump_193[:,ind2]/dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[0], label=r'UMP-193/CuCN3-254')
		ax[ind2].plot(pso2_axis, dose_ump_230[:,ind2]/dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[1], label=r'UMP-230/CuCN3-254')
		ax[ind2].plot(pso2_axis, dose_ump_254[:,ind2]/dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[2], label=r'UMP-254/CuCN3-254')

		ax[ind2].plot(pso2_axis, dose_ump_193[:,ind2]/dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[3], label=r'UMP-193/CuCN3-300')
		ax[ind2].plot(pso2_axis, dose_ump_230[:,ind2]/dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[4], label=r'UMP-230/CuCN3-300')
		ax[ind2].plot(pso2_axis, dose_ump_254[:,ind2]/dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[5], label=r'UMP-254/CuCN3-300')
		
		#ax.set_ylim([1.e-2, 1.e4])
		ax[ind2].set_yscale('log')
		ax[ind2].set_ylabel(r'$\bar{D}_{UMP-X}/\bar{D}_{CuCN3-Y}$')
	#ax.set_xlim([100, 500])
	ax[num_cases-1].set_xscale('log')
	ax[num_cases-1].set_xlabel(r'pSO$_2$', fontsize=12)
	
	
	plt.tight_layout(rect=(0,0,1., 0.9))
	ax[0].legend(bbox_to_anchor=[0, 1.2, 1., 0.7], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10) 
	plt.savefig('./Plots/discussion_doses_reldoses_pso2.pdf', orientation='portrait',papertype='letter', format='pdf')

	plt.show()

########################
###
########################
if plot_doses_ph2s_clouds:
	"""
	The purpose of this script is to plot the doses for a CO2-H2O-SO2 Martian atmosphere with pH2S, with CO2 cloud decks of varying thickness emplaced at 20.5 km (20-21 km). 
	#pCO2=0.02 bar (optically thin in gas scattering)
	#pH2S=2e-9 -- 2e-4 bar
	#CO2 cloud OD varies from 1-1000.
	
	#Conditions: T_0=250, A=desert, z=0, no TD XCs, DeltaScaling
	"""
	########################
	###Read in doses
	########################
	###Read in computed doses
	pso2_list=np.array(['0.02bar_250K_0so2_100ppbh2s','0.02bar_250K_0so2_1ppmh2s','0.02bar_250K_0so2_10ppmh2s', '0.02bar_250K_0so2_100ppmh2s', '0.02bar_250K_0so2_1000ppmh2s', '0.02bar_250K_0so2_10000ppmh2s']) #list of pH2S for pCO2=0.02
	pso2_axis=np.array([2.e-9, 2.e-8, 2.e-7, 2.e-6, 2.e-5, 2.e-4]) #pSO2 in bar
	
	cloud_list=np.array(['co2cloudod=1', 'co2cloudod=10', 'co2cloudod=100','co2cloudod=1000']) 
	cloud_labels=np.array([r'$\tau_{cloud}=1 (unscaled)$',r'$\tau_{cloud}=10$ (unscaled)',r'$\tau_{cloud}=100$ (unscaled)',r'$\tau_{cloud}=1000$ (unscaled)'])
	
	num_so2=len(pso2_list)
	num_cloud=len(cloud_list)
	
	dose_100_165=np.zeros([num_so2,num_cloud]) #surface radiance integrated 100-165 nm
	dose_200_300=np.zeros([num_so2,num_cloud]) #surface radiance integrated 200-300 nm
	dose_ump_193=np.zeros([num_so2,num_cloud]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=193
	dose_ump_230=np.zeros([num_so2,num_cloud]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=230
	dose_ump_254=np.zeros([num_so2,num_cloud]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=254
	dose_cucn3_254=np.zeros([num_so2,num_cloud]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=254
	dose_cucn3_300=np.zeros([num_so2,num_cloud]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=300
	
	for ind in range(0, num_so2):
		for ind2 in range(0, num_cloud):
			pso2=pso2_list[ind]
			cloudod=cloud_list[ind2]
			dose_100_165[ind,ind2],dose_200_300[ind,ind2],dose_ump_193[ind,ind2],dose_ump_230[ind,ind2],dose_ump_254[ind,ind2],dose_cucn3_254[ind,ind2],dose_cucn3_300[ind,ind2]=np.genfromtxt('./DoseRates/dose_rates_volcanicmars_'+pso2+'_z=0_A=desert_noTD_DS_'+cloudod+'_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)
		
	########################
	###Plot results
	########################
	fig, ax=plt.subplots(num_cloud, figsize=(16.5*cm2inch,10.), sharex=True, sharey=False)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,6))
	
	for ind2 in range(0, num_cloud):
		ax[ind2].set_title(cloud_labels[ind2])
		
		#ax[ind2].plot(pso2_axis, dose_200_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[0], label='Radiance 200-300 nm')
		ax[ind2].plot(pso2_axis, dose_ump_193[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[1], label=r'UMP Bond Cleavage ($\lambda_0=193$)')
		ax[ind2].plot(pso2_axis, dose_ump_230[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[2], label=r'UMP Bond Cleavage ($\lambda_0=230$)')
		ax[ind2].plot(pso2_axis, dose_ump_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[3], label=r'UMP Bond Cleavage ($\lambda_0=254$)')
		ax[ind2].plot(pso2_axis, dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[4], label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
		ax[ind2].plot(pso2_axis, dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[5], label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
			
		#ax[ind2].set_ylim([1.e-10, 1.e1])
		ax[ind2].set_yscale('log')
		ax[ind2].set_ylabel(r'Relative Dose Rate $\bar{D}_i$')
	#ax.set_xlim([100, 500])
	ax[num_cloud-1].set_xscale('log')
	ax[num_cloud-1].set_xlabel(r'pH$_2$S (unscaled)', fontsize=12)
	
	
	plt.tight_layout(rect=(0,0,1., 0.9))
	ax[0].legend(bbox_to_anchor=[0, 1.2, 1., 0.7], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10) 
	plt.savefig('./Plots/discussion_doses_ph2s_clouds.pdf', orientation='portrait',papertype='letter', format='pdf')

	plt.show()

########################
###
########################
if plot_doses_ph2s_pco2:
	"""
	The purpose of this script is to plot the doses for a CO2-H2O-H2S Martian atmosphere with varying pH2S for varying background pCO2
	#pCO2=0.02-2 bar
	#pH2S=2e-9 -- 2e-4 bar
	
	#Conditions: T_0=250, A=desert, z=0, no TD XCs, no DeltaScaling, no clouds
	"""
	########################
	###Read in doses
	########################
	###Read in computed doses
	ph2s_list0=np.array(['0.02bar_250K_0so2_100ppbh2s','0.02bar_250K_0so2_1ppmh2s','0.02bar_250K_0so2_10ppmh2s', '0.02bar_250K_0so2_100ppmh2s', '0.02bar_250K_0so2_1000ppmh2s', '0.02bar_250K_0so2_10000ppmh2s']) #list of pH2S for pCO2=0.02
	ph2s_list1=np.array(['0.2bar_250K_0so2_10ppbh2s','0.2bar_250K_0so2_100ppbh2s','0.2bar_250K_0so2_1ppmh2s', '0.2bar_250K_0so2_10ppmh2s', '0.2bar_250K_0so2_100ppmh2s','0.2bar_250K_0so2_1000ppmh2s']) #list of pH2S for pCO2=0.2
	ph2s_list2=np.array(['2bar_250K_0so2_1ppbh2s', '2bar_250K_0so2_10ppbh2s','2bar_250K_0so2_100ppbh2s','2bar_250K_0so2_1ppmh2s', '2bar_250K_0so2_10ppmh2s', '2bar_250K_0so2_100ppmh2s']) #list of pH2S for pCO2=2
	
	ph2s_axis=np.array([2.e-9, 2.e-8, 2.e-7, 2.e-6, 2.e-5, 2e-4]) #pSO2 in bar
	pco2_list=np.array(['0.02bar', '0.2bar', '2bar']) 
	pco2_labels=np.array(['0.02 bar', '0.2 bar', '2 bar'])
	
	num_ph2s=len(ph2s_axis)
	num_pco2=len(pco2_list)
	
	dose_100_165=np.zeros([num_ph2s,num_pco2]) #surface radiance integrated 100-165 nm
	dose_200_300=np.zeros([num_ph2s,num_pco2]) #surface radiance integrated 200-300 nm
	dose_ump_193=np.zeros([num_ph2s,num_pco2]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=193
	dose_ump_230=np.zeros([num_ph2s,num_pco2]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=230
	dose_ump_254=np.zeros([num_ph2s,num_pco2]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=254
	dose_cucn3_254=np.zeros([num_ph2s,num_pco2]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=254
	dose_cucn3_300=np.zeros([num_ph2s,num_pco2]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=300
	
	for ind in range(0, num_ph2s):
		dose_100_165[ind,0],dose_200_300[ind,0],dose_ump_193[ind,0],dose_ump_230[ind,0],dose_ump_254[ind,0],dose_cucn3_254[ind,0],dose_cucn3_300[ind,0]=np.genfromtxt('./DoseRates/dose_rates_volcanicmars_'+ph2s_list0[ind]+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)
		
		dose_100_165[ind,1],dose_200_300[ind,1],dose_ump_193[ind,1],dose_ump_230[ind,1],dose_ump_254[ind,1],dose_cucn3_254[ind,1],dose_cucn3_300[ind,1]=np.genfromtxt('./DoseRates/dose_rates_volcanicmars_'+ph2s_list1[ind]+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)
		
		dose_100_165[ind,2],dose_200_300[ind,2],dose_ump_193[ind,2],dose_ump_230[ind,2],dose_ump_254[ind,2],dose_cucn3_254[ind,2],dose_cucn3_300[ind,2]=np.genfromtxt('./DoseRates/dose_rates_volcanicmars_'+ph2s_list2[ind]+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)

	########################
	###Plot results
	########################
	fig, ax=plt.subplots(num_pco2, figsize=(16.5*cm2inch,10.), sharex=True, sharey=False)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,6))
	
	for ind2 in range(0, num_pco2):
		ax[ind2].set_title('pCO$_2$='+pco2_labels[ind2])
		
		#ax[ind2].plot(ph2s_axis, dose_200_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[0], label='Radiance 200-300 nm')
		ax[ind2].plot(ph2s_axis, dose_ump_193[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[1], label=r'UMP Bond Cleavage ($\lambda_0=193$)')
		ax[ind2].plot(ph2s_axis, dose_ump_230[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[2], label=r'UMP Bond Cleavage ($\lambda_0=230$)')
		ax[ind2].plot(ph2s_axis, dose_ump_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[3], label=r'UMP Bond Cleavage ($\lambda_0=254$)')
		ax[ind2].plot(ph2s_axis, dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[4], label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
		ax[ind2].plot(ph2s_axis, dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[5], label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
			
		#ax[ind2].set_ylim([1.e-, 1.e1])
		ax[ind2].set_yscale('log')
		ax[ind2].set_ylabel(r'Relative Dose Rate $\bar{D}_i$')
	#ax.set_xlim([100, 500])
	ax[num_pco2-1].set_xscale('log')
	ax[num_pco2-1].set_xlabel(r'pH$_2$S (bar)', fontsize=12)
	
	
	plt.tight_layout(rect=(0,0,1., 0.9))
	ax[0].legend(bbox_to_anchor=[0, 1.2, 1., 0.7], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10) 
	plt.savefig('./Plots/discussion_doses_ph2s_pco2.pdf', orientation='portrait',papertype='letter', format='pdf')

	plt.show()


########################
###
########################
if plot_doses_pso2_clouds:
	"""
	The purpose of this script is to plot the doses for a CO2-H2O-SO2 Martian atmosphere with pSO2, with CO2 cloud decks of varying thickness emplaced at 20.5 km (20-21 km). 
	#pCO2=0.02 bar (optically thin in gas scattering)
	#pSO2=2e-9-2e-5 bar
	#CO2 cloud OD varies from 1-1000.
	
	#Conditions: T_0=250, A=desert, z=0, no TD XCs, DeltaScaling
	"""
	########################
	###Read in doses
	########################
	###Read in computed doses
	pso2_list=np.array(['0.02bar_250K_100ppbso2_0h2s','0.02bar_250K_1ppmso2_0h2s','0.02bar_250K_10ppmso2_0h2s', '0.02bar_250K_100ppmso2_0h2s', '0.02bar_250K_1000ppmso2_0h2s']) #list of pSO2 for pCO2=0.02 
	pso2_axis=np.array([2.e-9, 2.e-8, 2.e-7, 2.e-6, 2.e-5]) #pSO2 in bar
	
	cloud_list=np.array(['co2cloudod=1', 'co2cloudod=10', 'co2cloudod=100','co2cloudod=1000']) 
	cloud_labels=np.array([r'$\tau_{cloud}=1 (unscaled)$',r'$\tau_{cloud}=10$ (unscaled)',r'$\tau_{cloud}=100$ (unscaled)',r'$\tau_{cloud}=1000$ (unscaled)'])
	
	num_so2=len(pso2_list)
	num_cloud=len(cloud_list)
	
	dose_100_165=np.zeros([num_so2,num_cloud]) #surface radiance integrated 100-165 nm
	dose_200_300=np.zeros([num_so2,num_cloud]) #surface radiance integrated 200-300 nm
	dose_ump_193=np.zeros([num_so2,num_cloud]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=193
	dose_ump_230=np.zeros([num_so2,num_cloud]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=230
	dose_ump_254=np.zeros([num_so2,num_cloud]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=254
	dose_cucn3_254=np.zeros([num_so2,num_cloud]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=254
	dose_cucn3_300=np.zeros([num_so2,num_cloud]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=300
	
	for ind in range(0, num_so2):
		for ind2 in range(0, num_cloud):
			pso2=pso2_list[ind]
			cloudod=cloud_list[ind2]
			dose_100_165[ind,ind2],dose_200_300[ind,ind2],dose_ump_193[ind,ind2],dose_ump_230[ind,ind2],dose_ump_254[ind,ind2],dose_cucn3_254[ind,ind2],dose_cucn3_300[ind,ind2]=np.genfromtxt('./DoseRates/dose_rates_volcanicmars_'+pso2+'_z=0_A=desert_noTD_DS_'+cloudod+'_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)
		
	########################
	###Plot results
	########################
	fig, ax=plt.subplots(num_cloud, figsize=(16.5*cm2inch,10.), sharex=True, sharey=False)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,6))
	
	for ind2 in range(0, num_cloud):
		ax[ind2].set_title(cloud_labels[ind2])
		
		#ax[ind2].plot(pso2_axis, dose_200_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[0], label='Radiance 200-300 nm')
		ax[ind2].plot(pso2_axis, dose_ump_193[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[1], label=r'UMP Bond Cleavage ($\lambda_0=193$)')
		ax[ind2].plot(pso2_axis, dose_ump_230[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[2], label=r'UMP Bond Cleavage ($\lambda_0=230$)')
		ax[ind2].plot(pso2_axis, dose_ump_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[3], label=r'UMP Bond Cleavage ($\lambda_0=254$)')
		ax[ind2].plot(pso2_axis, dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[4], label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
		ax[ind2].plot(pso2_axis, dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[5], label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
			
		#ax.set_ylim([1.e-2, 1.e4])
		ax[ind2].set_yscale('log')
		ax[ind2].set_ylabel(r'Relative Dose Rate $\bar{D}_i$')
	#ax.set_xlim([100, 500])
	ax[num_cloud-1].set_xscale('log')
	ax[num_cloud-1].set_xlabel(r'pSO$_2$', fontsize=12)
	
	
	plt.tight_layout(rect=(0,0,1., 0.9))
	ax[0].legend(bbox_to_anchor=[0, 1.2, 1., 0.7], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10) 
	plt.savefig('./Plots/discussion_doses_pso2_clouds.pdf', orientation='portrait',papertype='letter', format='pdf')

	plt.show()

########################
###
########################
if plot_doses_pso2_pco2:
	"""
	The purpose of this script is to plot the doses for a CO2-H2O-SO2 Martian atmosphere with varying pSO2 for varying background pCO2
	#pCO2=0.02-2 bar
	#pSO2=2e-9 -- 2e-5 bar
	
	#Conditions: T_0=250, A=desert, z=0, no TD XCs, no DeltaScaling, no clouds
	"""
	########################
	###Read in doses
	########################
	###Read in computed doses
	pso2_list0=np.array(['0.02bar_250K_100ppbso2_0h2s','0.02bar_250K_1ppmso2_0h2s','0.02bar_250K_10ppmso2_0h2s', '0.02bar_250K_100ppmso2_0h2s', '0.02bar_250K_1000ppmso2_0h2s']) #list of pSO2 for pCO2=0.02
	pso2_list1=np.array(['0.2bar_250K_10ppbso2_0h2s','0.2bar_250K_100ppbso2_0h2s','0.2bar_250K_1ppmso2_0h2s', '0.2bar_250K_10ppmso2_0h2s', '0.2bar_250K_100ppmso2_0h2s']) #list of pSO2 for pCO2=0.2
	pso2_list2=np.array(['2bar_250K_1ppbso2_0h2s','2bar_250K_10ppbso2_0h2s','2bar_250K_100ppbso2_0h2s', '2bar_250K_1ppmso2_0h2s', '2bar_250K_10ppmso2_0h2s']) #list of pSO2 for pCO2=2
	
	pso2_axis=np.array([2.e-9, 2.e-8, 2.e-7, 2.e-6, 2.e-5]) #pSO2 in bar
	pco2_list=np.array(['0.02bar', '0.2bar', '2bar']) 
	pco2_labels=np.array(['0.02 bar', '0.2 bar', '2 bar'])
	
	num_pso2=len(pso2_axis)
	num_pco2=len(pco2_list)
	
	dose_100_165=np.zeros([num_pso2,num_pco2]) #surface radiance integrated 100-165 nm
	dose_200_300=np.zeros([num_pso2,num_pco2]) #surface radiance integrated 200-300 nm
	dose_ump_193=np.zeros([num_pso2,num_pco2]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=193
	dose_ump_230=np.zeros([num_pso2,num_pco2]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=230
	dose_ump_254=np.zeros([num_pso2,num_pco2]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=254
	dose_cucn3_254=np.zeros([num_pso2,num_pco2]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=254
	dose_cucn3_300=np.zeros([num_pso2,num_pco2]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=300
	
	for ind in range(0, num_pso2):
		dose_100_165[ind,0],dose_200_300[ind,0],dose_ump_193[ind,0],dose_ump_230[ind,0],dose_ump_254[ind,0],dose_cucn3_254[ind,0],dose_cucn3_300[ind,0]=np.genfromtxt('./DoseRates/dose_rates_volcanicmars_'+pso2_list0[ind]+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)
		
		dose_100_165[ind,1],dose_200_300[ind,1],dose_ump_193[ind,1],dose_ump_230[ind,1],dose_ump_254[ind,1],dose_cucn3_254[ind,1],dose_cucn3_300[ind,1]=np.genfromtxt('./DoseRates/dose_rates_volcanicmars_'+pso2_list1[ind]+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)
		
		dose_100_165[ind,2],dose_200_300[ind,2],dose_ump_193[ind,2],dose_ump_230[ind,2],dose_ump_254[ind,2],dose_cucn3_254[ind,2],dose_cucn3_300[ind,2]=np.genfromtxt('./DoseRates/dose_rates_volcanicmars_'+pso2_list2[ind]+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)

	########################
	###Plot results
	########################
	fig, ax=plt.subplots(num_pco2, figsize=(16.5*cm2inch,10.), sharex=True, sharey=False)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,6))
	
	for ind2 in range(0, num_pco2):
		ax[ind2].set_title('pCO$_2$='+pco2_labels[ind2])
		
		#ax[ind2].plot(pso2_axis, dose_200_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[0], label='Radiance 200-300 nm')
		ax[ind2].plot(pso2_axis, dose_ump_193[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[1], label=r'UMP Bond Cleavage ($\lambda_0=193$)')
		ax[ind2].plot(pso2_axis, dose_ump_230[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[2], label=r'UMP Bond Cleavage ($\lambda_0=230$)')
		ax[ind2].plot(pso2_axis, dose_ump_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[3], label=r'UMP Bond Cleavage ($\lambda_0=254$)')
		ax[ind2].plot(pso2_axis, dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[4], label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
		ax[ind2].plot(pso2_axis, dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[5], label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
			
		#ax.set_ylim([1.e-2, 1.e4])
		ax[ind2].set_yscale('log')
		ax[ind2].set_ylabel(r'Relative Dose Rate $\bar{D}_i$')
	#ax.set_xlim([100, 500])
	ax[num_pco2-1].set_xscale('log')
	ax[num_pco2-1].set_xlabel(r'pSO$_2$ (bar)', fontsize=12)
	
	
	plt.tight_layout(rect=(0,0,1., 0.9))
	ax[0].legend(bbox_to_anchor=[0, 1.2, 1., 0.7], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10) 
	plt.savefig('./Plots/discussion_doses_pso2_pco2.pdf', orientation='portrait',papertype='letter', format='pdf')

	plt.show()

########################
###
########################
if plot_doses_dust_clouds:
	"""
	The purpose of this script is to plot the doses for a cold, dry (CO2+H2O) Martian atmosphere with varying levels of dust for varying background CO2 cloud levels. 
	#pCO2=0.02 bar
	#dust OD=0.1-10
	#cloud OD=1-1000
	
	#Conditions: T_0=250, A=desert, z=0, no TD XCs, DeltaScaling
	"""
	########################
	###Read in doses
	########################
	###Read in computed doses
	od_list=np.array(['dustod=0.1', 'dustod=1', 'dustod=10']) #list of dust optical depths (500 nm)
	od_axis=np.array([1.e-1, 1., 1.e1])
	cloud_list=np.array(['co2cloudod=1', 'co2cloudod=10', 'co2cloudod=100','co2cloudod=1000']) 
	cloud_labels=np.array([r'$\tau_{cloud}=1 (unscaled)$',r'$\tau_{cloud}=10$ (unscaled)',r'$\tau_{cloud}=100$ (unscaled)',r'$\tau_{cloud}=1000$ (unscaled)'])
	
	num_od=len(od_list)
	num_cloud=len(cloud_list)
	
	dose_100_165=np.zeros([num_od,num_cloud]) #surface radiance integrated 100-165 nm
	dose_200_300=np.zeros([num_od,num_cloud]) #surface radiance integrated 200-300 nm
	dose_ump_193=np.zeros([num_od,num_cloud]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=193
	dose_ump_230=np.zeros([num_od,num_cloud]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=230
	dose_ump_254=np.zeros([num_od,num_cloud]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=254
	dose_cucn3_254=np.zeros([num_od,num_cloud]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=254
	dose_cucn3_300=np.zeros([num_od,num_cloud]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=300
	
	for ind in range(0, num_od):
		for ind2 in range(0, num_cloud):
			od=od_list[ind]
			cloudod=cloud_list[ind2]
			dose_100_165[ind,ind2],dose_200_300[ind,ind2],dose_ump_193[ind,ind2],dose_ump_230[ind,ind2],dose_ump_254[ind,ind2],dose_cucn3_254[ind,ind2],dose_cucn3_300[ind,ind2]=np.genfromtxt('./DoseRates/dose_rates_colddrymars_0.02bar_250K_z=0_A=desert_noTD_DS_'+od+'_'+cloudod+'_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)
		
	########################
	###Plot results
	########################
	fig, ax=plt.subplots(num_cloud, figsize=(16.5*cm2inch,10.), sharex=True, sharey=False)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,6))
	
	for ind2 in range(0, num_cloud):
		ax[ind2].set_title(cloud_labels[ind2])
		
		#ax[ind2].plot(od_axis, dose_200_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[0], label='Radiance 200-300 nm')
		ax[ind2].plot(od_axis, dose_ump_193[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[1], label=r'UMP Bond Cleavage ($\lambda_0=193$)')
		ax[ind2].plot(od_axis, dose_ump_230[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[2], label=r'UMP Bond Cleavage ($\lambda_0=230$)')
		ax[ind2].plot(od_axis, dose_ump_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[3], label=r'UMP Bond Cleavage ($\lambda_0=254$)')
		ax[ind2].plot(od_axis, dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[4], label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
		ax[ind2].plot(od_axis, dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[5], label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
			
		#ax.set_ylim([1.e-2, 1.e4])
		ax[ind2].set_yscale('log')
		ax[ind2].set_ylabel(r'Relative Dose Rate $\bar{D}_i$')
	#ax.set_xlim([100, 500])
	ax[num_cloud-1].set_xscale('log')
	ax[num_cloud-1].set_xlabel(r'$\tau_{d}$ (unscaled)', fontsize=12)
	
	
	plt.tight_layout(rect=(0,0,1., 0.9))
	ax[0].legend(bbox_to_anchor=[0, 1.2, 1., 0.7], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10) 
	plt.savefig('./Plots/discussion_doses_dust_clouds.pdf', orientation='portrait',papertype='letter', format='pdf')

	plt.show()

########################
###
########################
if plot_doses_dust_pco2:
	"""
	The purpose of this script is to plot the doses for a cold, dry (CO2+H2O) Martian atmosphere with varying levels of dust for varying background pCO2
	#pCO2=0.02-2 bar
	#dust OD=0.1-10 (100--> absorption optical epth of >20 -->just no flux anywhere)
	
	#Conditions: T_0=250, A=desert, z=0, no TD XCs, DeltaScaling, no clouds
	"""
	########################
	###Read in doses
	########################
	###Read in computed doses
	od_list=np.array(['dustod=0.1', 'dustod=1', 'dustod=10']) #list of dust optical depths (500 nm)
	od_axis=np.array([1.e-1, 1., 1.e1])
	pco2_list=np.array(['0.02bar', '0.2bar', '2bar']) 
	pco2_labels=np.array(['0.02 bar', '0.2 bar', '2 bar'])
	
	num_od=len(od_list)
	num_pco2=len(pco2_list)
	
	dose_100_165=np.zeros([num_od,num_pco2]) #surface radiance integrated 100-165 nm
	dose_200_300=np.zeros([num_od,num_pco2]) #surface radiance integrated 200-300 nm
	dose_ump_193=np.zeros([num_od,num_pco2]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=193
	dose_ump_230=np.zeros([num_od,num_pco2]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=230
	dose_ump_254=np.zeros([num_od,num_pco2]) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=254
	dose_cucn3_254=np.zeros([num_od,num_pco2]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=254
	dose_cucn3_300=np.zeros([num_od,num_pco2]) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=300
	
	for ind in range(0, num_od):
		for ind2 in range(0, num_pco2):
			od=od_list[ind]
			pco2=pco2_list[ind2]
			dose_100_165[ind,ind2],dose_200_300[ind,ind2],dose_ump_193[ind,ind2],dose_ump_230[ind,ind2],dose_ump_254[ind,ind2],dose_cucn3_254[ind,ind2],dose_cucn3_300[ind,ind2]=np.genfromtxt('./DoseRates/dose_rates_colddrymars_'+pco2+'_250K_z=0_A=desert_noTD_DS_'+od+'.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)
		
	########################
	###Plot results
	########################
	fig, ax=plt.subplots(num_pco2, figsize=(16.5*cm2inch,10.), sharex=True, sharey=False)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,6))
	
	for ind2 in range(0, num_pco2):
		ax[ind2].set_title('pCO$_2$='+pco2_labels[ind2])
		
		#ax[ind2].plot(od_axis, dose_200_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[0], label='Radiance 200-300 nm')
		ax[ind2].plot(od_axis, dose_ump_193[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[1], label=r'UMP Bond Cleavage ($\lambda_0=193$)')
		ax[ind2].plot(od_axis, dose_ump_230[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[2], label=r'UMP Bond Cleavage ($\lambda_0=230$)')
		ax[ind2].plot(od_axis, dose_ump_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[3], label=r'UMP Bond Cleavage ($\lambda_0=254$)')
		ax[ind2].plot(od_axis, dose_cucn3_254[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[4], label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
		ax[ind2].plot(od_axis, dose_cucn3_300[:,ind2], marker='s', markersize=markersizeval, linewidth=1, color=colors[5], label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
			
		#ax.set_ylim([1.e-2, 1.e4])
		ax[ind2].set_yscale('log')
		ax[ind2].set_ylabel(r'Relative Dose Rate $\bar{D}_i$')
	#ax.set_xlim([100, 500])
	ax[num_pco2-1].set_xscale('log')
	ax[num_pco2-1].set_xlabel(r'$\tau_{d}$ (unscaled)', fontsize=12)
	
	
	plt.tight_layout(rect=(0,0,1., 0.9))
	ax[0].legend(bbox_to_anchor=[0, 1.2, 1., 0.7], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10) 
	plt.savefig('./Plots/discussion_doses_dust_pco2.pdf', orientation='portrait',papertype='letter', format='pdf')

	plt.show()


########################
###
########################
if plot_doses_clouds:
	"""
	The purpose of this script is to plot the doses for a cold, dry (CO2+H2O) Martian atmosphere with a cloud deck of varying optical thickness at 20-21 km altitude.
	Cloud optical depths: 1-10000 (specified at 500 nm).
	
	#Conditions: pCO2=2 bar, T_0=250, A=desert, SZA=0, no TD XCs, yes DeltaScaling
	NOTE: In final run, may want to use A=tundra to reduce reference to weird effects
	NOTE: In final run, may want to use pCO2=0.02 to deconvolve cloud vs atmosphere.
	"""
	########################
	###Read in doses
	########################
	###Read in computed doses
	elt_list=np.array(['0.1', '1', '10', '100', '1000', '10000']) #list of cloud optical depths (500 nm)
	od_list=np.array([1.e-1, 1., 1.e1, 1.e2, 1.e3, 1.e4])
	
	dose_100_165=np.zeros(np.shape(elt_list)) #surface radiance integrated 100-165 nm
	dose_200_300=np.zeros(np.shape(elt_list)) #surface radiance integrated 200-300 nm
	dose_ump_193=np.zeros(np.shape(elt_list)) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=193
	dose_ump_230=np.zeros(np.shape(elt_list)) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=230
	dose_ump_254=np.zeros(np.shape(elt_list)) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=254
	dose_cucn3_254=np.zeros(np.shape(elt_list)) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=254
	dose_cucn3_300=np.zeros(np.shape(elt_list)) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=300
	
	num_elt=len(elt_list)
	for ind in range(0, num_elt):
		elt=elt_list[ind]
		dose_100_165[ind],dose_200_300[ind],dose_ump_193[ind],dose_ump_230[ind],dose_ump_254[ind],dose_cucn3_254[ind],dose_cucn3_300[ind]=np.genfromtxt('./DoseRates/dose_rates_colddrymars_0.02bar_250K_z=0_A=desert_noTD_DS_co2cloudod='+elt+'_z=20.5_reff=10.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)

	########################
	###Plot results
	########################
	fig, ax=plt.subplots(1, figsize=(16.5*cm2inch,6.), sharex=True, sharey=True)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,6))
	
	#ax.plot(od_list, dose_200_300, marker='s', markersize=markersizeval, linewidth=1, color=colors[0], label='Radiance 200-300 nm')
	ax.plot(od_list, dose_ump_193, marker='s', markersize=markersizeval, linewidth=1, color=colors[1], label=r'UMP Bond Cleavage ($\lambda_0=193$)')
	ax.plot(od_list, dose_ump_230, marker='s', markersize=markersizeval, linewidth=1, color=colors[2], label=r'UMP Bond Cleavage ($\lambda_0=230$)')
	ax.plot(od_list, dose_ump_254, marker='s', markersize=markersizeval, linewidth=1, color=colors[3], label=r'UMP Bond Cleavage ($\lambda_0=254$)')
	ax.plot(od_list, dose_cucn3_254, marker='s', markersize=markersizeval, linewidth=1, color=colors[4], label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax.plot(od_list, dose_cucn3_300, marker='s', markersize=markersizeval, linewidth=1, color=colors[5], label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
		
	#ax.set_ylim([1.e-2, 1.e4])
	ax.set_yscale('log')
	ax.set_ylabel(r'Relative Dose Rate $\bar{D}_i$')
	#ax.set_xlim([100, 500])
	ax.set_xscale('log')
	ax.set_xlabel(r'$\tau_{cloud}$ (unscaled)', fontsize=12)
	
	
	plt.tight_layout(rect=(0,0,1., 0.85))
	ax.legend(bbox_to_anchor=[0, 1.05, 1., 0.7], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10) 
	plt.savefig('./Plots/discussion_doses_co2clouds.pdf', orientation='portrait',papertype='letter', format='pdf')

	plt.show()


########################
###
########################
if plot_doses_pco2:
	"""
	The purpose of this script is to plot the doses for a cold, dry (CO2+H2O) Martian atmosphere for different surface pressures (measure dose as function of pCO2). 
	Surface Pressures (CO2 dry partial pressure): 2, 0.2, 0.02 bar (250K), 2e-3, 2e-4, 2e-5 bar (200 K)
	
	A=desert, SZA=0 (maximum fluence conditions)
	No delta scaling (no clouds), no TD XCs
	NOTE: In final run, may want to use A=tundra to reduce reference to weird effects
	"""
	########################
	###Read in doses
	########################
	###Read in computed doses
	elt_list=np.array(['0.00002bar_200K', '0.0002bar_200K', '0.002bar_200K', '0.02bar_250K', '0.2bar_250K','2bar_250K']) #list of file name identifiers
	P_list=np.array([2.e-5, 2.e-4, 2.e-3, 2.e-2, 2.e-1, 2.]) #list of pressures associated with each file name
	
	dose_100_165=np.zeros(np.shape(elt_list)) #surface radiance integrated 100-165 nm
	dose_200_300=np.zeros(np.shape(elt_list)) #surface radiance integrated 200-300 nm
	dose_ump_193=np.zeros(np.shape(elt_list)) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=193
	dose_ump_230=np.zeros(np.shape(elt_list)) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=230
	dose_ump_254=np.zeros(np.shape(elt_list)) #dose rate for UMP glycosidic bond cleavage, assuming lambda0=254
	dose_cucn3_254=np.zeros(np.shape(elt_list)) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=254
	dose_cucn3_300=np.zeros(np.shape(elt_list)) #dose rate for solvated electron production from tricyanocuprate, assuming lambda0=300
	
	num_elt=len(elt_list)
	for ind in range(0, num_elt):
		elt=elt_list[ind]
		dose_100_165[ind],dose_200_300[ind],dose_ump_193[ind],dose_ump_230[ind],dose_ump_254[ind],dose_cucn3_254[ind],dose_cucn3_300[ind]=np.genfromtxt('./DoseRates/dose_rates_colddrymars_'+elt+'_z=0_A=desert_noTD_noDS_noparticles.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)

	
	########################
	###Plot results
	########################
	fig, ax=plt.subplots(1, figsize=(16.5*cm2inch,6.), sharex=True, sharey=True)
	markersizeval=5.
	colors=cm.rainbow(np.linspace(0,1,6))
	
	#ax.plot(P_list, dose_200_300, marker='s', markersize=markersizeval, linewidth=1, color=colors[0], label='Radiance 200-300 nm')
	ax.plot(P_list, dose_ump_193, marker='s', markersize=markersizeval, linewidth=1, color=colors[1], label=r'UMP Bond Cleavage ($\lambda_0=193$)')
	ax.plot(P_list, dose_ump_230, marker='s', markersize=markersizeval, linewidth=1, color=colors[2], label=r'UMP Bond Cleavage ($\lambda_0=230$)')
	ax.plot(P_list, dose_ump_254, marker='s', markersize=markersizeval, linewidth=1, color=colors[3], label=r'UMP Bond Cleavage ($\lambda_0=254$)')
	ax.plot(P_list, dose_cucn3_254, marker='s', markersize=markersizeval, linewidth=1, color=colors[4], label=r'CuCN$_3$ Photoionization ($\lambda_0=254$)')
	ax.plot(P_list, dose_cucn3_300, marker='s', markersize=markersizeval, linewidth=1, color=colors[5], label=r'CuCN$_3$ Photoionization ($\lambda_0=300$)')
	
	###ax.axvline(6.e-3, color='black', linestyle='--', linewidth='1') #demarcates breakpoint between two temperature regimes. Doesn't matter, but still good to be clear when variable changes.
	
	#ax.set_ylim([1.e-2, 1.e4])
	ax.set_yscale('linear')
	ax.set_ylabel(r'Relative Dose Rate $\bar{D}_i$')
	#ax.set_xlim([100, 500])
	ax.set_xscale('log')
	ax.set_xlabel('pCO$_2$ (bar)')
	
	
	plt.tight_layout(rect=(0,0,1., 0.85))
	ax.legend(bbox_to_anchor=[0, 1.05, 1., 0.7], loc=3, ncol=2, mode='expand', borderaxespad=0., fontsize=10) 
	plt.savefig('./Plots/discussion_doses_pco2.pdf', orientation='portrait',papertype='letter', format='pdf')

	plt.show()