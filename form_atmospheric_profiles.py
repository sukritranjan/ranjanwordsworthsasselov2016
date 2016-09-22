# -*- coding: iso-8859-1 -*-
"""
This file contains function definitions and calls that create the atmospheric profile files (T, P, and gas molar concentrations as a function of altitude), that form inputs into our radiative transfer codes. There are two broad families of such files:

I) VALIDATION CASES: these include form_profiles_primitive_earth_rugheimer and form_profiles_wuttke. These functions, once called, generate atmospheric profiles files as well as files containing the TOA solar input that we can use to reproduce the calculations of Rugheimer et al (2015) and the measurements of Wuttke et al (2006), as validation cases.

II) RESEARCH CASES: these are the calls that create the feedstock files used in our study. They include form_spectral_feedstock_youngmars, which is used to define a uniform TOA solar flux file for the young Mars. They also include calls to generate_profiles_cold_dry_mars and generate_profiles_volcanic_mars, which are defined in the file mars_atmosphere_models.py, to give the atmospheric profiles for atmospheres with user specified P_0(CO2), T_0, and specified SO2 and H2S loading levels (latter file only)

All file generation calls are at the end of the respective section.
"""
import numpy as np
import pdb
import matplotlib.pyplot as plt
import scipy.stats
from scipy import interpolate as interp
import cookbook
import mars_atmosphere_models as mars

bar2Ba=1.0e6 #1 bar in Ba
k=1.3806488e-16 #Boltzmann Constant in erg/K


############################################
###I. VALIDATION CASES
############################################
def form_profiles_primitive_earth_rugheimer():
	"""
	Purpose of this code is to form spectra, mixing ratio files, and T/P profiles for the revised Rugheimer Epoch 0 (3.9 Ga) Earth models (Rugheimer et al 2015).
	"""
	
	#####Zeroth: set value of constants, specify filenames
	filename='./Raw_Data/Rugheimer_Metadata/outchem_Ep0_A0.2_Frac1.0.dat'
	
	
	#####First, form the spectra for comparison.
	importeddata=np.genfromtxt(filename, skip_header=290, skip_footer=1277)
	
	#Remove the first wavelength bin which corresponds to Lyman Alpha and which does not have a bin width that fits with its neighbors.
	rugheimer_wav_centers=importeddata[1:,1]/10. #Convert wavelengths from Angstroms to nm
	rugheimer_s=importeddata[1:,4] #ratio of 4piJ(surf)/I_0
	rugheimer_s[19]=3.16548e-128 #one element of rugheimer_s has value 3.16548e-128. Python has trouble with this and imports as a NaN. Here, we manually set its value.

	###Form wavelength bins from Rugheimer wavelength centers
	rugheimer_wav_bin_leftedges=np.zeros(len(rugheimer_wav_centers))
	rugheimer_wav_bin_rightedges=np.zeros(len(rugheimer_wav_centers))

	#First ten FUV fluxes are 5 nm (50 A) bins (email from srugheimer@gmail.com, 3/12/2015) 
	rugheimer_wav_bin_leftedges[0:9]=rugheimer_wav_centers[0:9]-2.5
	rugheimer_wav_bin_rightedges[0:9]=rugheimer_wav_centers[0:9]+2.5

	#Remainder of FUV fluxes are taken from a file that sarah sent me (srugheimer@gmail.com, 3/12/2015)
	del importeddata
	importeddata=np.genfromtxt('./Raw_Data/Rugheimer_Metadata/Active_M9_Teff2300_photo.pdat', skip_header=1, skip_footer=0)

	rugheimer_wav_bin_leftedges[9:]=importeddata[:,2]*0.1 #convert A to nm
	rugheimer_wav_bin_rightedges[9:]=importeddata[:,3]*0.1 #convert A to nm

	####Check that bins are correct:
	###print np.sum(rugheimer_wav_centers-0.5*(rugheimer_wav_bin_leftedges+rugheimer_wav_bin_rightedges)) #0 to within 1e-12 rounding error.

	###Rebin Claire et al input.
	#Import 0.01-nm resolution Claire et al 3.9 Ga Sun model.
	del importeddata
	importeddata=np.genfromtxt('./Raw_Data/Claire_Model/claire_youngsun_highres.dat', skip_header=1, skip_footer=0)
	claire_wav=importeddata[:,0] #nm, 0.01 nm resolution
	claire_fluxes=importeddata[:,1]#erg/s/cm2/nm

	#Bin Claire et al model to resolution of Rugheimer model
	claire_fluxes_rebinned=np.zeros(len(rugheimer_wav_centers))
	claire_wav_rebinned=np.zeros(len(claire_fluxes_rebinned))#This should be redundant with rugheimer_wav_centers. We include it as a check statistic that the rebinning is proceeding appropriately. 

	for ind in range(0, len(rugheimer_wav_centers)):
		min_wav=rugheimer_wav_bin_leftedges[ind]
		max_wav=rugheimer_wav_bin_rightedges[ind]
		inds=(claire_wav >= min_wav) & (claire_wav <= max_wav)
		claire_fluxes_rebinned[ind]=np.mean(claire_fluxes[inds])
		claire_wav_rebinned[ind]=np.mean(claire_wav[inds]) #check statistic.

	###print np.sum((claire_wav_rebinned-rugheimer_wav_centers)/rugheimer_wav_centers) #check statistic. Good to within 1e-5 in all cases. Any problems caused by slight misalignment from 0.01 due to rounding error. 

	###Compute bottom-of-atmosphere actinic flux, which is what is reported in Rugheimer+2015.
	rugheimer_ground_energies=claire_fluxes_rebinned*rugheimer_s
	
	#Let's print out the results
	spectable=np.zeros([len(rugheimer_wav_bin_leftedges),4])
	spectable[:,0]=rugheimer_wav_bin_leftedges
	spectable[:,1]=rugheimer_wav_bin_rightedges
	spectable[:,2]=rugheimer_wav_centers
	spectable[:,3]=claire_fluxes_rebinned
	
	header='Left Bin Edge (nm)	Right Bin Edge (nm)	Bin Center (nm)		Solar Flux at Earth (erg/s/nm/cm2)\n'

	f=open('./Solar_Input/general_youngsun_earth_spectral_input.dat', 'w')
	f.write(header)
	np.savetxt(f, spectable, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()
	
	###########################################################################################
	###########################################################################################
	###########################################################################################
	
	#####Second, form the mixing ratio files
	importeddata1=np.genfromtxt(filename, skip_header=779, skip_footer=873) #O2, O3, H2O
	importeddata2=np.genfromtxt(filename, skip_header=837, skip_footer=817) #CH4, SO2
	importeddata4=np.genfromtxt(filename, skip_header=958, skip_footer=704) #N2, CO2
	
	#Let's print out the results. We have established that the z values are the same, so can use a common block
	printtable=np.zeros([np.shape(importeddata1)[0],9])
	printtable[:,0]=importeddata1[:,0] #altitude in cm
	#N2 and CO2: We use the values from this block rather than block 1 because rugheimer et al force it to these values in their code, regardless of what the photochemistry code wants to do.
	printtable[:,1]=importeddata4[:,2] #N2. 
	printtable[:,2]=importeddata4[:,1] #CO2
	#The rest are normal
	printtable[:,3]=importeddata1[:,3] #H2O
	printtable[:,4]=importeddata2[:,2] #CH4
	printtable[:,5]=importeddata2[:,9] #SO2
	printtable[:,6]=importeddata1[:,2] #O2
	printtable[:,7]=importeddata1[:,8] #O3
	#printtable[:,8]# H2S; left as zeros since not included in Rugheimer model
	
	#print np.sum(printtable[:,1:],1)
	#pdb.set_trace()
	
	header0='Extracted from Rugheimer outchem_Ep0_A0.2_Frac1.0.dat\n'
	header1='Z (cm)		N2	CO2	H2O	CH4	SO2	O2	O3	H2S \n'

	f=open('./MolarConcentrations/rugheimer_earth_epoch0_molarconcentrations.dat', 'w')
	f.write(header0)
	f.write(header1)
	np.savetxt(f, printtable, delimiter='	', fmt='%1.7e', newline='\n')
	f.close()	

	###########################################################################################
	###########################################################################################
	###########################################################################################
	#####Third, form the T/P profiles
	
	#Extract temperature and pressure profile from climate model output
	#For whatever reason the very last line of the table is doubled. We remove this. 
	importeddata=np.genfromtxt(filename, skip_header=1568, skip_footer=104)
	model_z=importeddata[:-1,0] #altitude in cm
	model_t=importeddata[:-1,1] #temperature in K
	model_n=importeddata[:-1,3] #number density in cm**-3.
	model_p=importeddata[:-1,4] #pressure, in bar (based on text in draft manuscript sent to me by Sarah Rugheimer)

	#Let's print out the results
	printtable=np.zeros([len(model_z)+1,4])
	printtable[1:,0]=model_z
	printtable[1:,1]=model_t
	printtable[1:,2]=model_n
	printtable[1:,3]=model_p
	
	#Rugheimer data file does not explicitly include t, P, n at z=0 (Surface). Our code requires z=0 data. To reconcile, we include these data manually as follows:
	printtable[0,0]=0. #z=0 case
	printtable[0,3]=1. #In the paper, p=1.0 bar at surface is specified
	printtable[0,1]=292.95 #From linear extrapolation from z=0.5 km and z=1.5 km points
	printtable[0,2]= 1.*bar2Ba/(k*292.95)#Compute number density self-consistently from temperature, pressure via Ideal Gas Law as is done elsewhere (n [cm**-3] = p [Barye]/(k*T [K])
	
	header0='Extracted from Rugheimer outchem_Ep0_A0.2_Frac1.0.dat\n'
	header1='Z (cm)	T (K)	DEN (cm**-3)	P (bar) \n'

	f=open('./TPProfiles/rugheimer_earth_epoch0_tpprofile.dat', 'w')
	f.write(header0)
	f.write(header1)
	np.savetxt(f, printtable, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()


def form_profiles_wuttke():
	"""
	Purpose of this code is to form the feedstock files to replicat the Wuttke+2006 Antarctic diffuse radiance measurements.
	
	To replicate these observations, we need a model of the modern atmosphere. We use the atmosphere Rugheimer et al (2013) calculated for the modern Earth. The feedstock files for these measurements were taken from Rugheimer et al (2013) via Ranjan & Sasselov (2016b).
	"""
	#First, form the spectral file.
	#Define spectral bins. 0.25 nm from 280-500 nm, 1 nm from 500-1000 nm. We just go to 900 since that's what our data is good to. Also we start at 292.75 because that's where our graphclicked data starts
	bin_left_edges=np.concatenate((np.arange(292.75,500.,0.25),np.arange(500., 900.,1.)))
	bin_right_edges=np.concatenate((np.arange(293.,500.25,0.25),np.arange(501., 901.,1.)))
	bin_centers=0.5*(bin_left_edges+bin_right_edges)
	
	
	#load BOA diffuse zenith flux from Wuttke+2006 (extracted via GraphClick)	
	importeddata=np.genfromtxt('./Raw_Data/UV_Surface_Measurements/wuttke.csv', skip_header=0, skip_footer=0, delimiter=',')
	dif_wav=importeddata[:,0] #nm
	dif_flux=importeddata[:,1]*2.*np.pi #mW/m2/nm/sr=erg/s/cm2/nm/sr; multiply by 2pi to convert to hemisphere-integrated total surface diffuse radiances
	dif_func=interp.interp1d(dif_wav, dif_flux, kind='linear')
	dif_flux_interp=dif_func(bin_centers)
	
	datatable=np.zeros([len(bin_centers), 2])
	datatable[:,0]=bin_centers
	datatable[:,1]=dif_flux_interp
	
	f=open('./LiteratureSpectra/wuttke2006.dat', 'w')
	np.savetxt(f, datatable, delimiter='		', fmt='%1.7e', newline='\n', header='Bin Center (nm)	Zenith Diffuse Flux (erg/s/nm/cm2)')

	
	#load solar spectrum from Claire et al (2012) models, normalized to 1 au 
	importeddata=np.genfromtxt('./Raw_Data/Claire_Model/claire_modernsun_highres.dat', skip_header=1, skip_footer=0)
	claire_wav=importeddata[:,0] #nm, 0.1 nm resolution, 100-900 nm.
	claire_fluxes=importeddata[:,1]#erg/s/cm2/nm
	
	#rebin claire spectrum
	claire_fluxes_rebinned=cookbook.rebin_uneven(np.arange(99.995,900.005,0.01), np.arange(100.005, 900.015,0.01),claire_fluxes,bin_left_edges, bin_right_edges)   
	
	#Plot to make sure rebinning worked correctly
	fig, ax1=plt.subplots(1, figsize=(6,4))
	ax1.plot(claire_wav, claire_fluxes, marker='s', color='black', label='Claire Fluxes')
	ax1.plot(bin_centers, claire_fluxes_rebinned, marker='s', color='blue', label='Binned Claire Fluxes')	
	ax1.set_yscale('log')
	ax1.set_ylim([1.e-2, 1.e4])
	ax1.set_xlim([280.,900.])
	ax1.set_xlabel('nm')
	ax1.set_ylabel('erg/s/cm2/nm')
	ax1.legend(loc=0)
	plt.show()	
	
	#Let's print out the results
	spectable=np.zeros([len(bin_left_edges),4])
	spectable[:,0]=bin_left_edges
	spectable[:,1]=bin_right_edges
	spectable[:,2]=bin_centers
	spectable[:,3]=claire_fluxes_rebinned
	
	header='Left Bin Edge (nm)	Right Bin Edge (nm)	Bin Center (nm)		Top of Atm Flux (erg/s/nm/cm2)\n'

	f=open('./Solar_Input/modernsun_earth_wuttke2006.dat', 'w')
	f.write(header)
	np.savetxt(f, spectable, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()
	
	###########################################################################################
	###########################################################################################
	###########################################################################################
	#####Second, form the molar concentration files
	#####Form by replicating the Rugheimer modern Earth profile, then scaling down the H2O level and scaling up the O3 level.
	mixingratios=np.genfromtxt('./MolarConcentrations/rugheimer_earth_modern_mixingratios_v2.dat', skip_header=2, skip_footer=0)
	mixingratios[:,3]=mixingratios[:,3]*0.1 #scale down h2o by factor of 10
	mixingratios[:,7]=mixingratios[:,7]*1.25 #scale up ozone by factor of 1.25
	header0='Based on Rugheimer+2013 Modern Earth Model\n'
	header1='Z (cm)		N2	CO2	H2O	CH4	SO2	O2	O3	H2S\n'

	f=open('./MolarConcentrations/wuttke2006_molarconcentrations.dat', 'w')
	f.write(header0)
	f.write(header1)
	np.savetxt(f, mixingratios, delimiter='	', fmt='%1.7e', newline='\n')
	f.close()
	###########################################################################################
	###########################################################################################
	###########################################################################################
	#####Finally, form TP profile
	#####Form by duplicating Rugheimer+2013 modern Earth profile
	tpprofile=np.genfromtxt('./TPProfiles/rugheimer_earth_modern_atmosphereprofile.dat', skip_header=2, skip_footer=0)
	header0='Based on Rugheimer+2013 Modern Earth Model\n'
	header1='Z (cm)	T (K)	DEN (cm**-3)	P (bar) \n'

	f=open('./TPProfiles/wuttke2006_tpprofile.dat', 'w')
	f.write(header0)
	f.write(header1)
	np.savetxt(f, tpprofile, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()

###Call Validation Cases
#form_profiles_primitive_earth_rugheimer()
#form_profiles_wuttke()

############################################
###II. RESEARCH CASES
############################################
def form_spectral_feedstock_youngmars():
	"""
	Purpose of this code is to form the spectral feedstock file (TOA solar flux) to explore formally the dependence of UV surface radiance on various factors.
	"""
	#Define spectral bins.
	bin_left_edges=np.arange(100.,500.,1.)
	bin_right_edges=np.arange(101.,501.,1.)
	bin_centers=0.5*(bin_left_edges+bin_right_edges)
	
	#load solar spectrum at 3.9 Ga from Claire et al (2012) models, normalized to 1 AU. These are really TOA intensities. Multiply by mu_0 to get TOA fluxes. 
	importeddata=np.genfromtxt('./Raw_Data/Claire_Model/claire_youngsun_highres.dat', skip_header=1, skip_footer=0)
	claire_wav=importeddata[:,0] #nm, 0.01 nm resolution, 100-900 nm.
	claire_fluxes=importeddata[:,1]*(1./1.524)**2 #Scale the flux to the Martian semimajor axis, units of erg/s/cm2/nm
	
	#rebin claire spectrum
	claire_fluxes_rebinned=cookbook.rebin_uneven(np.arange(99.995,900.005,0.01), np.arange(100.005, 900.015,0.01),claire_fluxes,bin_left_edges, bin_right_edges)   
	
	#Plot to make sure rebinning worked correctly
	fig, ax1=plt.subplots(1, figsize=(6,4))
	ax1.plot(claire_wav, claire_fluxes, marker='s', color='black', label='Claire Fluxes')
	ax1.plot(bin_centers, claire_fluxes_rebinned, marker='s', color='blue', label='Binned Claire Fluxes')	
	ax1.set_yscale('log')
	ax1.set_ylim([1.e-2, 1.e4])
	ax1.set_xlim([100.,500.])
	ax1.set_xlabel('nm')
	ax1.set_ylabel('erg/s/cm2/nm')
	ax1.legend(loc=0)
	plt.show()	
	
	#Let's print out the results
	spectable=np.zeros([len(bin_left_edges),4])
	spectable[:,0]=bin_left_edges
	spectable[:,1]=bin_right_edges
	spectable[:,2]=bin_centers
	spectable[:,3]=claire_fluxes_rebinned
	
	header='Left Bin Edge (nm)	Right Bin Edge (nm)	Bin Center (nm)		Top of Atm Intensity (erg/s/nm/cm2)\n'

	f=open('./Solar_Input/general_youngsun_mars_spectral_input.dat', 'w')
	f.write(header)
	np.savetxt(f, spectable, delimiter='		', fmt='%1.7e', newline='\n')
	f.close()
#form_spectral_feedstock_youngmars() #form general TOA solar flux file

#Form CO2-H2O atmosphere profile files
mars.generate_profiles_cold_dry_mars(2., 210., 'colddrymars_2bar_210K')
mars.generate_profiles_cold_dry_mars(2., 250., 'colddrymars_2bar_250K')
mars.generate_profiles_cold_dry_mars(2., 300., 'colddrymars_2bar_300K')

mars.generate_profiles_cold_dry_mars(0.2, 210., 'colddrymars_0.2bar_210K')
mars.generate_profiles_cold_dry_mars(0.2, 250., 'colddrymars_0.2bar_250K')
mars.generate_profiles_cold_dry_mars(0.2, 300., 'colddrymars_0.2bar_300K')

mars.generate_profiles_cold_dry_mars(0.02, 210., 'colddrymars_0.02bar_210K')
mars.generate_profiles_cold_dry_mars(0.02, 250., 'colddrymars_0.02bar_250K')
mars.generate_profiles_cold_dry_mars(0.02, 300., 'colddrymars_0.02bar_300K')

mars.generate_profiles_cold_dry_mars(0.002, 200., 'colddrymars_0.002bar_200K')
mars.generate_profiles_cold_dry_mars(0.0002, 200., 'colddrymars_0.0002bar_200K')
mars.generate_profiles_cold_dry_mars(0.00002, 200., 'colddrymars_0.00002bar_200K')

#Form CO2-H2O-SO2/H2S atmosphere profile files
mars.generate_profiles_volcanic_mars(2., 250., 1.e-9, 0., 'volcanicmars_2bar_250K_1ppbso2_0h2s')
mars.generate_profiles_volcanic_mars(2., 250., 1.e-8, 0., 'volcanicmars_2bar_250K_10ppbso2_0h2s')
mars.generate_profiles_volcanic_mars(2., 250., 1.e-7, 0., 'volcanicmars_2bar_250K_100ppbso2_0h2s')
mars.generate_profiles_volcanic_mars(2., 250., 1.e-6, 0., 'volcanicmars_2bar_250K_1ppmso2_0h2s')
mars.generate_profiles_volcanic_mars(2., 250., 1.e-5, 0., 'volcanicmars_2bar_250K_10ppmso2_0h2s')

mars.generate_profiles_volcanic_mars(2., 250., 0., 1.e-9, 'volcanicmars_2bar_250K_0so2_1ppbh2s')
mars.generate_profiles_volcanic_mars(2., 250., 0., 1.e-8, 'volcanicmars_2bar_250K_0so2_10ppbh2s')
mars.generate_profiles_volcanic_mars(2., 250., 0., 1.e-7, 'volcanicmars_2bar_250K_0so2_100ppbh2s')
mars.generate_profiles_volcanic_mars(2., 250., 0., 1.e-6, 'volcanicmars_2bar_250K_0so2_1ppmh2s')
mars.generate_profiles_volcanic_mars(2., 250., 0., 1.e-5, 'volcanicmars_2bar_250K_0so2_10ppmh2s')
mars.generate_profiles_volcanic_mars(2., 250., 0., 1.e-4, 'volcanicmars_2bar_250K_0so2_100ppmh2s')

mars.generate_profiles_volcanic_mars(0.2, 250., 1.e-8, 0., 'volcanicmars_0.2bar_250K_10ppbso2_0h2s')
mars.generate_profiles_volcanic_mars(0.2, 250., 1.e-7, 0., 'volcanicmars_0.2bar_250K_100ppbso2_0h2s')
mars.generate_profiles_volcanic_mars(0.2, 250., 1.e-6, 0., 'volcanicmars_0.2bar_250K_1ppmso2_0h2s')
mars.generate_profiles_volcanic_mars(0.2, 250., 1.e-5, 0., 'volcanicmars_0.2bar_250K_10ppmso2_0h2s')
mars.generate_profiles_volcanic_mars(0.2, 250., 1.e-4, 0., 'volcanicmars_0.2bar_250K_100ppmso2_0h2s')

mars.generate_profiles_volcanic_mars(0.2, 250., 0., 1.e-8, 'volcanicmars_0.2bar_250K_0so2_10ppbh2s')
mars.generate_profiles_volcanic_mars(0.2, 250., 0., 1.e-7, 'volcanicmars_0.2bar_250K_0so2_100ppbh2s')
mars.generate_profiles_volcanic_mars(0.2, 250., 0., 1.e-6, 'volcanicmars_0.2bar_250K_0so2_1ppmh2s')
mars.generate_profiles_volcanic_mars(0.2, 250., 0., 1.e-5, 'volcanicmars_0.2bar_250K_0so2_10ppmh2s')
mars.generate_profiles_volcanic_mars(0.2, 250., 0., 1.e-4, 'volcanicmars_0.2bar_250K_0so2_100ppmh2s')
mars.generate_profiles_volcanic_mars(0.2, 250., 0., 1.e-3, 'volcanicmars_0.2bar_250K_0so2_1000ppmh2s')


mars.generate_profiles_volcanic_mars(0.02, 250., 1.e-8, 0., 'volcanicmars_0.02bar_250K_10ppbso2_0h2s')
mars.generate_profiles_volcanic_mars(0.02, 250., 1.e-7, 0., 'volcanicmars_0.02bar_250K_100ppbso2_0h2s')
mars.generate_profiles_volcanic_mars(0.02, 250., 1.e-6, 0., 'volcanicmars_0.02bar_250K_1ppmso2_0h2s')
mars.generate_profiles_volcanic_mars(0.02, 250., 1.e-5, 0., 'volcanicmars_0.02bar_250K_10ppmso2_0h2s')
mars.generate_profiles_volcanic_mars(0.02, 250., 1.e-4, 0., 'volcanicmars_0.02bar_250K_100ppmso2_0h2s')
mars.generate_profiles_volcanic_mars(0.02, 250., 1.e-3, 0., 'volcanicmars_0.02bar_250K_1000ppmso2_0h2s')

mars.generate_profiles_volcanic_mars(0.02, 250., 0., 1.e-8, 'volcanicmars_0.02bar_250K_0so2_10ppbh2s')
mars.generate_profiles_volcanic_mars(0.02, 250., 0., 1.e-7, 'volcanicmars_0.02bar_250K_0so2_100ppbh2s')
mars.generate_profiles_volcanic_mars(0.02, 250., 0., 1.e-6, 'volcanicmars_0.02bar_250K_0so2_1ppmh2s')
mars.generate_profiles_volcanic_mars(0.02, 250., 0., 1.e-5, 'volcanicmars_0.02bar_250K_0so2_10ppmh2s')
mars.generate_profiles_volcanic_mars(0.02, 250., 0., 1.e-4, 'volcanicmars_0.02bar_250K_0so2_100ppmh2s')
mars.generate_profiles_volcanic_mars(0.02, 250., 0., 1.e-3, 'volcanicmars_0.02bar_250K_0so2_1000ppmh2s')
mars.generate_profiles_volcanic_mars(0.02, 250., 0., 1.e-2, 'volcanicmars_0.02bar_250K_0so2_10000ppmh2s')

#mars.generate_profiles_volcanic_mars(2., 205., 1.e-7, 0., 'volcanicmars_2bar_205K_100ppbso2_0h2s')
#mars.generate_profiles_volcanic_mars(0.02, 205., 1.e-5, 0., 'volcanicmars_0.02bar_205K_10ppmso2_0h2s')
