# -*- coding: iso-8859-1 -*-
"""
This code demonstrates our reproduction of previous model results, including 1) Rugheimer et al (2015, 2013), 2) the Wuttke measurement and 3) the WOUDC measurement
"""
########################
###Import useful libraries
########################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pdb
import cookbook

def cm2inch(cm): #function to convert cm to inches; useful for complying with Astrobiology size guidelines
	return cm/2.54

########################
###Which plots to make:
########################
reproduce_rugheimer_youngEarth=True
reproduce_wuttke=True


########################
###Reproduce Rugheimer+2015 Young Earth Model
########################

if reproduce_rugheimer_youngEarth:
	########################
	###Read in Models
	########################
	###Read in Rugheimer model output
	rugheimer_wav_leftedges, rugheimer_wav_rightedges, rugheimer_wav_centers, rugheimer_int_toa, rugheimer_int_boa=np.genfromtxt('./LiteratureSpectra/rugheimer_epoch0_recomputed_A0.2.dat', skip_header=1, skip_footer=0, usecols=(0,1,2,3,4), unpack=True) # 0:left edges of wavelength bins, nm | 1:right edges of wavelength bins, nm | 2:centers of wavelength bins, nm | 3:top-of-atmosphere intensity, erg/s/cm2/nm | 4:BOA Actinic Flux/Total Intensity, erg/s/cm2/nm. This is what we need to match.


	###Load in our code's reproduction (A=0.20, our cross-sections)
	twostr1_wav_leftedges, twostr1_wav_rightedges, twostr1_wav_centers,twostr1_int_toa, twostr1_flux_gnd, twostr1_int_BOA, twostr1_int_gnd=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_epoch0_reproduce_rugheimer.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True)# 0:left edges of wavelength bins, nm | 1:right edges of wavelength bins, nm | 2:centers of wavelength bins, nm | 3:top-of-atmosphere intensity, erg/s/cm2/nm | 4:#surface flux, erg/s/cm2/nm | 5:total intensity in middle of bottom layer of atmosphere | 6:intensity at planetary surface, =0.5*diffuse total intensity+solar intensity

	
	###Load in our code's reproduction (A=0.20, Rugheimer cross-sections)
	twostr3_wav_leftedges, twostr3_wav_rightedges, twostr3_wav_centers,twostr3_int_toa, twostr3_flux_gnd, twostr3_int_BOA, twostr3_int_gnd=np.genfromtxt('./TwoStreamOutput/rugheimer_earth_epoch0_reproduce_rugheimer_rugheimerXCs.dat', skip_header=1, skip_footer=0,usecols=(0,1,2,3,4,5,6), unpack=True) # 0:left edges of wavelength bins, nm | 1:right edges of wavelength bins, nm | 2:centers of wavelength bins, nm | 3:top-of-atmosphere intensity, erg/s/cm2/nm | 4:surface flux, erg/s/cm2/nm | 5:total intensity in middle of bottom layer of atmosphere | 6:intensity at planetary surface, =0.5*diffuse total intensity+solar intensity

	########################
	###Plot results
	########################
	fig, (ax0, ax1)=plt.subplots(2, figsize=(cm2inch(16.5),10.), sharex=True)
	markersizeval=5.

	ax0.plot(rugheimer_wav_centers, rugheimer_int_toa,  marker='s', markersize=markersizeval, linewidth=1, color='black', label='TOA Int.')
	ax0.plot(rugheimer_wav_centers, rugheimer_int_boa,  marker='s', markersize=markersizeval, linewidth=1, color='red', label='BOA Int. \n(Rugheimer+2015)')
	ax0.plot(twostr1_wav_centers, twostr1_int_BOA,  marker='s', markersize=markersizeval, linewidth=1, color='blue', label='BOA Int. (Our Model)')
	ax0.plot(twostr3_wav_centers, twostr3_int_BOA,  marker='s', markersize=markersizeval, linewidth=1, color='orange', label='BOA Int. (Our Model, \nRugheimer+2015 XCs)')
	ax0.set_yscale('log')
	ax0.set_ylim([1.e-2, 1e4])
	ax0.set_ylabel('Total Intensity \n (erg s$^{-1}$cm$^{-2}$nm$^{-1}$)')
	ax0.legend(loc='lower center', fontsize=12)
	ax0.set_xlim([130, 860])


	#if both our codes report zero, ratio should be 1.
	rugheimer_int_boa_plot_3=np.copy(rugheimer_int_boa)
	rugheimer_int_boa_plot_1=np.copy(rugheimer_int_boa)
	twostr1_int_BOA_plot=np.copy(twostr1_int_BOA)
	twostr3_int_BOA_plot=np.copy(twostr3_int_BOA)

	foo=np.where((rugheimer_int_boa==0.) & (twostr3_int_BOA==0.))
	twostr3_int_BOA_plot[foo]=1.
	rugheimer_int_boa_plot_3[foo]=1.

	bar=np.where((rugheimer_int_boa==0.) & (twostr1_int_BOA==0.))
	twostr1_int_BOA_plot[bar]=1.
	rugheimer_int_boa_plot_1[bar]=1.


	ax1.plot(twostr1_wav_centers, (twostr1_int_BOA-rugheimer_int_boa)/rugheimer_int_toa, marker='s', linewidth=1, color='blue', label='Our XCs')
	ax1.plot(twostr3_wav_centers, (twostr3_int_BOA-rugheimer_int_boa)/rugheimer_int_toa, marker='s', linewidth=1, color='orange', label='Rugheimer+2015 XCs')
	ax1.set_yscale('linear')
	ax1.set_ylim([-0.3, 0.1])
	ax1.set_xlabel('Wavelength (nm)')
	ax1.set_ylabel('(Our Model-R+2015)/TOA Int.')
	ax1.legend(loc='lower center', fontsize=12)
	ax1.set_xlim([130, 860])

	ax0.set_title('Reproduction of Rugheimer+2015 Model')
	bunk=(twostr1_int_BOA_plot-rugheimer_int_boa_plot_1)/rugheimer_int_toa
	thunk=(twostr3_int_BOA_plot-rugheimer_int_boa_plot_3)/rugheimer_int_toa
	print np.max(np.abs(bunk))
	print np.max(np.abs(thunk))
	plt.tight_layout()
	plt.savefig('./Plots/reproduce_rugheimer_epoch0.eps', orientation='portrait',papertype='letter', format='eps')
	plt.show()
	
########################
###Reproduce Wuttke+2006 modern Earth measurements
########################

if reproduce_wuttke:
	########################
	###Read in Models
	########################
	###Read in TOA input
	input_wav, input_toa=np.genfromtxt('./Solar_Input/modernsun_earth_wuttke2006.dat', skip_header=1, skip_footer=0, usecols=(2,3), unpack=True) # 0:left edges of wavelength bins, nm | centers of wavelength bins, nm | top-of-atmosphere intensity, erg/s/cm2/nm 
	
	wuttke_wav, wuttke_surf_measured=np.genfromtxt('./LiteratureSpectra/wuttke2006.dat', skip_header=1, skip_footer=0, usecols=(0,1), unpack=True)


	###Load in our code's reproduction
	twostr_wav_centers, twostr_toa, twostr_surf_diffint=np.genfromtxt('./TwoStreamOutput/wuttke2006_reproduce_wuttke2006.dat', skip_header=1, skip_footer=0,usecols=(2,3,7), unpack=True)# 0:left edges of wavelength bins, nm | 1:right edges of wavelength bins, nm | 2:centers of wavelength bins, nm | 3:top-of-atmosphere intensity, erg/s/cm2/nm | 4:#surface flux, erg/s/cm2/nm | 5:total intensity in middle of bottom layer of atmosphere | 6:intensity at planetary surface, =0.5*diffuse total intensity+solar intensity

	twostr_surf_diffint_smoothed=cookbook.movingaverage(twostr_surf_diffint, 10)
	
	########################
	###Some quick checks to make sure that things that should be the same are in fact the same
	########################
	print np.max(np.abs((input_wav-twostr_wav_centers)/input_wav))
	print np.max(np.abs((input_wav-wuttke_wav)/input_wav))
	#are the wavelength scales the same?
	
	print np.max(np.abs((input_toa-twostr_toa)/input_toa))
	#are the TOA intensities the same?
	
	########################
	###Plot results
	########################
	fig, (ax1, ax2)=plt.subplots(2, figsize=(cm2inch(16.5),10.), sharex=True)
	markersizeval=5.

	ax1.plot(input_wav, input_toa, marker='s', color='black', label='TOA Intensity')
	ax1.plot(twostr_wav_centers, twostr_surf_diffint_smoothed, marker='s', color='red', label='Calculated Diffuse Radiance')
	ax1.plot(wuttke_wav, wuttke_surf_measured, marker='s', color='blue', label='Measured Diffuse Radiance')
	ax1.set_yscale('log')
	ax1.set_ylim([1.e-3, 1.e4])
	ax1.set_xlim([290.,900.])
	ax1.set_xlabel('nm')
	ax1.set_ylabel('erg/s/cm2/nm')
	ax1.legend(loc=0)
	ax2.plot(wuttke_wav, (wuttke_surf_measured-twostr_surf_diffint_smoothed)/(input_toa), marker='s', color='black', label='(Measured BOA-Modelled BOA)/TOA')
	ax2.set_xlabel('nm')
	ax2.set_yscale('linear')
	ax2.set_ylim([-0.3,.3])
	ax2.legend(loc=0)

	plt.savefig('./Plots/reproduce_wuttke_paper.eps', orientation='portrait',papertype='letter', format='eps')
	plt.show()