# -*- coding: iso-8859-1 -*-
"""
Purpose of this file is to run the uv_radtrans function from radiativetransfer.py to generate the surface radiance calculations used to derive the results in our paper.

This script: run for the Mars-approximate case (A=desert)
"""

import radiativetransfer as rt
import numpy as np
import pdb
##################################
###First, validation cases
##################################
####Replicate Rugheimer+2015 3.9 Ga young Earth calculation.
#rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel='rugheimer_earth_epoch0', outputfilelabel='reproduce_rugheimer', inputspectrafile='general_youngsun_earth_spectral_input.dat',TDXC=False, DeltaScaling=False, SZA_deg=60., albedoflag='uniformalbedo',uniformalbedo=0.2, includedust=False, includeco2cloud=False,includeh2ocloud=False)
##Replicate Rugheimer+2015 3.9 Ga young Earth calculation using Rugheimer cross-sections
rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel='rugheimer_earth_epoch0', outputfilelabel='reproduce_rugheimer_rugheimerXCs', inputspectrafile='general_youngsun_earth_spectral_input.dat',TDXC=False, DeltaScaling=False, SZA_deg=60., albedoflag='uniformalbedo',uniformalbedo=0.2, includedust=False, includeco2cloud=False,includeh2ocloud=False) ###NEED TO UNCOMMENT RUGHEIMER XCs BLOCK IN CODE -- BE SURE TO RECOMMEND AFTERWARDS

#####Replicate Wuttke+2006 surface radiance measurement calculation
#rt.uv_radtrans(z_upper_limit=60.e5, z_step=60.e3, inputatmofilelabel='wuttke2006', outputfilelabel='reproduce_wuttke2006', inputspectrafile='modernsun_earth_wuttke2006.dat',TDXC=False, DeltaScaling=False, SZA_deg=51.2, albedoflag='nonuniformalbedo',nonuniformalbedo=np.array([0.,0.,0.,0.,1.]), includedust=False, includeco2cloud=False,includeh2ocloud=False)

#######################################################################################################################
#CO2-H2O Atmosphere (no clouds/dust)
#######################################################################################################################

#################
####Particulate-free atmospheres, no TDXC, no delta-scaling, variable P_0 and T_0, SZA=0, A=desert
#################
inputatmofilelabel_elt_list=np.array(['colddrymars_0.02bar_210K', 'colddrymars_0.02bar_250K', 'colddrymars_0.02bar_300K', 'colddrymars_0.2bar_210K', 'colddrymars_0.2bar_250K', 'colddrymars_0.2bar_300K','colddrymars_2bar_210K', 'colddrymars_2bar_250K', 'colddrymars_2bar_300K'])

for inputatmofilelabel_elt in inputatmofilelabel_elt_list:
	rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel=inputatmofilelabel_elt, outputfilelabel='z=0_A=desert_noTD_noDS_noparticles', inputspectrafile='general_youngsun_mars_spectral_input.dat',TDXC=False, DeltaScaling=False, SZA_deg=0., albedoflag='nonuniformalbedo',nonuniformalbedo=np.array([0.,0.,1.,0.,0.]), includedust=False, includeco2cloud=False,includeh2ocloud=False)


#################
####Particulate-free atmospheres, T_0=200, P_0=2e-5-2e-3, SZA=0, A=desert, with and without TD XCs
#################
inputatmofilelabel_elt_list=np.array(['colddrymars_0.00002bar_200K', 'colddrymars_0.0002bar_200K', 'colddrymars_0.002bar_200K'])

for inputatmofilelabel_elt in inputatmofilelabel_elt_list:
	#Without TDXCs:
	rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel=inputatmofilelabel_elt, outputfilelabel='z=0_A=desert_noTD_noDS_noparticles', inputspectrafile='general_youngsun_mars_spectral_input.dat',TDXC=False, DeltaScaling=False, SZA_deg=0., albedoflag='nonuniformalbedo',nonuniformalbedo=np.array([0.,0.,1.,0.,0.]), includedust=False, includeco2cloud=False,includeh2ocloud=False)

	#with TDXCs
	rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel=inputatmofilelabel_elt, outputfilelabel='z=0_A=desert_TD_noDS_noparticles', inputspectrafile='general_youngsun_mars_spectral_input.dat',TDXC=True, DeltaScaling=False, SZA_deg=0., albedoflag='nonuniformalbedo',nonuniformalbedo=np.array([0.,0.,1.,0.,0.]), includedust=False, includeco2cloud=False,includeh2ocloud=False)


#########################################################################################################################
###CO2-H2O Atmosphere (clouds)
#########################################################################################################################
##################
#####Atmospheres with H2O and CO2 clouds, T_0=250, P_0=0.02, SZA=0, A=desert, no TDXCs, yes delta scaling
##################

cloudtaulabel_list=np.array(['0.1', '1', '10', '100', '1000', '10000']) #list of cloud optical depths (500 nm)
cloudtauvalue_list=np.array([0.1, 1., 10., 100., 1000., 10000.]) #list of cloud optical depths (500 nm)

for ind in range(0, len(cloudtaulabel_list)):
	cloudtaulabel=cloudtaulabel_list[ind]
	cloudtauvalue=cloudtauvalue_list[ind]
	
	#H2O Clouds
	rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel='colddrymars_0.02bar_250K', outputfilelabel='z=0_A=desert_noTD_DS_h2ocloudod='+cloudtaulabel+'_z=3.5_reff=10', inputspectrafile='general_youngsun_mars_spectral_input.dat',TDXC=False, DeltaScaling=True, SZA_deg=0., albedoflag='nonuniformalbedo',nonuniformalbedo=np.array([0.,0.,1.,0.,0.]), includedust=False, includeco2cloud=False,includeh2ocloud=True,h2ocloudlayerinds=np.array([60]), h2ocloudlayerods=np.array([cloudtauvalue]),h2oiceparamsfile='cloud_h2o_reff10_vareff0p1_lognormal.pickle')


	##CO2 Clouds
	rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel='colddrymars_0.02bar_250K', outputfilelabel='z=0_A=desert_noTD_DS_co2cloudod='+cloudtaulabel+'_z=20.5_reff=10', inputspectrafile='general_youngsun_mars_spectral_input.dat',TDXC=False, DeltaScaling=True, SZA_deg=0., albedoflag='nonuniformalbedo',nonuniformalbedo=np.array([0.,0.,1.,0.,0.]), includedust=False, includeh2ocloud=False,includeco2cloud=True,co2cloudlayerinds=np.array([43]), co2cloudlayerods=np.array([cloudtauvalue]),co2iceparamsfile='cloud_co2_reff10_vareff0p1_lognormal.pickle')

#########################################################################################################################
###CO2-H2O-SO2/H2S atmospheres. No particulate, no TDXC, no delta-scaling, T_0=250, pCO2=0.02-2bar, SZA=0, A=fresh snow
#########################################################################################################################

#################
###SO2 (pSO2=2e-9 -- 2e-5 bar for pCO2=2e-2-2 bar)
#################
inputatmofilelabel_elt_list=np.array(['volcanicmars_0.02bar_250K_100ppbso2_0h2s', 'volcanicmars_0.02bar_250K_1ppmso2_0h2s', 'volcanicmars_0.02bar_250K_10ppmso2_0h2s', 'volcanicmars_0.02bar_250K_100ppmso2_0h2s', 'volcanicmars_0.02bar_250K_1000ppmso2_0h2s','volcanicmars_0.2bar_250K_10ppbso2_0h2s', 'volcanicmars_0.2bar_250K_100ppbso2_0h2s', 'volcanicmars_0.2bar_250K_1ppmso2_0h2s', 'volcanicmars_0.2bar_250K_10ppmso2_0h2s', 'volcanicmars_0.2bar_250K_100ppmso2_0h2s','volcanicmars_2bar_250K_1ppbso2_0h2s', 'volcanicmars_2bar_250K_10ppbso2_0h2s', 'volcanicmars_2bar_250K_100ppbso2_0h2s', 'volcanicmars_2bar_250K_1ppmso2_0h2s', 'volcanicmars_2bar_250K_10ppmso2_0h2s'])

for inputatmofilelabel_elt in inputatmofilelabel_elt_list:
	rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel=inputatmofilelabel_elt, outputfilelabel='z=0_A=desert_noTD_noDS_noparticles', inputspectrafile='general_youngsun_mars_spectral_input.dat',TDXC=False, DeltaScaling=False, SZA_deg=0., albedoflag='nonuniformalbedo',nonuniformalbedo=np.array([0.,0.,1.,0.,0.]), includedust=False, includeh2ocloud=False,includeco2cloud=False)


#################
###H2S (pH2S=2e-9 -- 2e-4 bar for pCO2=2e-2-2 bar)
#################
inputatmofilelabel_elt_list=np.array(['volcanicmars_0.02bar_250K_0so2_100ppbh2s', 'volcanicmars_0.02bar_250K_0so2_1ppmh2s', 'volcanicmars_0.02bar_250K_0so2_10ppmh2s', 'volcanicmars_0.02bar_250K_0so2_100ppmh2s', 'volcanicmars_0.02bar_250K_0so2_1000ppmh2s', 'volcanicmars_0.02bar_250K_0so2_10000ppmh2s', 'volcanicmars_0.2bar_250K_0so2_10ppbh2s', 'volcanicmars_0.2bar_250K_0so2_100ppbh2s', 'volcanicmars_0.2bar_250K_0so2_1ppmh2s', 'volcanicmars_0.2bar_250K_0so2_10ppmh2s', 'volcanicmars_0.2bar_250K_0so2_100ppmh2s', 'volcanicmars_0.2bar_250K_0so2_1000ppmh2s', 'volcanicmars_2bar_250K_0so2_1ppbh2s', 'volcanicmars_2bar_250K_0so2_10ppbh2s', 'volcanicmars_2bar_250K_0so2_100ppbh2s', 'volcanicmars_2bar_250K_0so2_1ppmh2s', 'volcanicmars_2bar_250K_0so2_10ppmh2s', 'volcanicmars_2bar_250K_0so2_100ppmh2s'])

for inputatmofilelabel_elt in inputatmofilelabel_elt_list:
	rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel=inputatmofilelabel_elt, outputfilelabel='z=0_A=desert_noTD_noDS_noparticles', inputspectrafile='general_youngsun_mars_spectral_input.dat',TDXC=False, DeltaScaling=False, SZA_deg=0., albedoflag='nonuniformalbedo',nonuniformalbedo=np.array([0.,0.,1.,0.,0.]), includedust=False, includeh2ocloud=False,includeco2cloud=False)



#########################################################################################################################
###CO2-H2O-SO2/H2S atmospheres. CO2 clouds @ 20-21 km (variable thickness), no TDXC, no delta-scaling, T_0=250, P_0=0.02, SZA=0, A=fresh snow
#########################################################################################################################

##################
####SO2 (pSO2=2e-9 -- 2e-5 bar for pCO2=2e-2 bar, CO2 cloud OD=1-1000)
##################
inputatmofilelabel_elt_list=np.array(['volcanicmars_0.02bar_250K_100ppbso2_0h2s', 'volcanicmars_0.02bar_250K_1ppmso2_0h2s', 'volcanicmars_0.02bar_250K_10ppmso2_0h2s', 'volcanicmars_0.02bar_250K_100ppmso2_0h2s', 'volcanicmars_0.02bar_250K_1000ppmso2_0h2s'])
cloudtaulabel_list=np.array(['1', '10', '100', '1000']) #list of cloud optical depths (500 nm)
cloudtauvalue_list=np.array([1., 10., 100., 1000.]) #list of cloud optical depths (500 nm)

for inputatmofilelabel_elt in inputatmofilelabel_elt_list:
	for ind in range(0, len(cloudtaulabel_list)):
		cloudtaulabel=cloudtaulabel_list[ind]
		cloudtauvalue=cloudtauvalue_list[ind]
		
		rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel=inputatmofilelabel_elt, outputfilelabel='z=0_A=desert_noTD_DS_co2cloudod='+cloudtaulabel+'_z=20.5_reff=10', inputspectrafile='general_youngsun_mars_spectral_input.dat',TDXC=False, DeltaScaling=True, SZA_deg=0., albedoflag='nonuniformalbedo',nonuniformalbedo=np.array([0.,0.,1.,0.,0.]), includedust=False, includeh2ocloud=False,includeco2cloud=True,co2cloudlayerinds=np.array([43]), co2cloudlayerods=np.array([cloudtauvalue]),co2iceparamsfile='cloud_co2_reff10_vareff0p1_lognormal.pickle')

##################
####H2S (pSO2=2e-9 -- 2e-4 bar for pCO2=2e-2 bar, CO2 cloud OD=1-1000)
##################
inputatmofilelabel_elt_list=np.array(['volcanicmars_0.02bar_250K_0so2_100ppbh2s', 'volcanicmars_0.02bar_250K_0so2_1ppmh2s', 'volcanicmars_0.02bar_250K_0so2_10ppmh2s', 'volcanicmars_0.02bar_250K_0so2_100ppmh2s', 'volcanicmars_0.02bar_250K_0so2_1000ppmh2s', 'volcanicmars_0.02bar_250K_0so2_10000ppmh2s'])
cloudtaulabel_list=np.array(['1', '10', '100', '1000']) #list of cloud optical depths (500 nm)
cloudtauvalue_list=np.array([1., 10., 100., 1000.]) #list of cloud optical depths (500 nm)

for inputatmofilelabel_elt in inputatmofilelabel_elt_list:
	for ind in range(0, len(cloudtaulabel_list)):
		cloudtaulabel=cloudtaulabel_list[ind]
		cloudtauvalue=cloudtauvalue_list[ind]
		
		rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel=inputatmofilelabel_elt, outputfilelabel='z=0_A=desert_noTD_DS_co2cloudod='+cloudtaulabel+'_z=20.5_reff=10', inputspectrafile='general_youngsun_mars_spectral_input.dat',TDXC=False, DeltaScaling=True, SZA_deg=0., albedoflag='nonuniformalbedo',nonuniformalbedo=np.array([0.,0.,1.,0.,0.]), includedust=False, includeh2ocloud=False,includeco2cloud=True,co2cloudlayerinds=np.array([43]), co2cloudlayerods=np.array([cloudtauvalue]),co2iceparamsfile='cloud_co2_reff10_vareff0p1_lognormal.pickle')

#########################################################################################################################
###CO2-H2O atmospheres. Varying levels of exponentially distributed dust, no TDXC, no delta-scaling, T_0=250, P_0=0.02, SZA=0, A=fresh snow
#########################################################################################################################

##################
####Just dust (tau_d=0.01-10, pCO2=2e-2--2 bar)
##################
inputatmofilelabel_elt_list=np.array(['colddrymars_0.02bar_250K','colddrymars_0.2bar_250K','colddrymars_2bar_250K'])
dusttaulabel_list=np.array(['0.1', '1', '10']) #list of total dust optical depths (500 nm)
dusttauvalue_list=np.array([0.1, 1., 10.]) #list of total dust optical depths (500 nm)

for inputatmofilelabel_elt in inputatmofilelabel_elt_list:
	for ind in range(0, len(dusttaulabel_list)):
		dusttaulabel=dusttaulabel_list[ind]
		dusttauvalue=dusttauvalue_list[ind]
		
		rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel=inputatmofilelabel_elt, outputfilelabel='z=0_A=desert_noTD_DS_dustod='+dusttaulabel, inputspectrafile='general_youngsun_mars_spectral_input.dat',TDXC=False, DeltaScaling=True, SZA_deg=0., albedoflag='nonuniformalbedo',nonuniformalbedo=np.array([0.,0.,1.,0.,0.]), includeh2ocloud=False,includeco2cloud=False, includedust=True,tau_d=dusttauvalue,dustparamsfile='dust_wolff_pangajello_reff1p5_vareff0p5_lognormal.pickle')

##################
####Dust and clouds (Varying tau_d and tau_cloud)
##################
dusttaulabel_list=np.array(['0.1', '1', '10']) #list of total dust optical depths (500 nm)
dusttauvalue_list=np.array([0.1, 1., 10.]) #list of total dust optical depths (500 nm)

cloudtaulabel_list=np.array(['1', '10', '100', '1000']) #list of cloud optical depths (500 nm)
cloudtauvalue_list=np.array([1., 10., 100., 1000.]) #list of cloud optical depths (500 nm)

for ind_dust in range(0, len(dusttaulabel_list)):
	for ind_cloud in range(0, len(cloudtaulabel_list)):
		dusttaulabel=dusttaulabel_list[ind_dust]
		dusttauvalue=dusttauvalue_list[ind_dust]
		
		cloudtaulabel=cloudtaulabel_list[ind_cloud]
		cloudtauvalue=cloudtauvalue_list[ind_cloud]
		
		rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel='colddrymars_0.02bar_250K', outputfilelabel='z=0_A=desert_noTD_DS_dustod='+dusttaulabel+'_co2cloudod='+cloudtaulabel+'_z=20.5_reff=10', inputspectrafile='general_youngsun_mars_spectral_input.dat',TDXC=False, DeltaScaling=True, SZA_deg=0., albedoflag='nonuniformalbedo',nonuniformalbedo=np.array([0.,0.,1.,0.,0.]), includeh2ocloud=False,includeco2cloud=True,co2cloudlayerinds=np.array([43]), co2cloudlayerods=np.array([cloudtauvalue]),co2iceparamsfile='cloud_co2_reff10_vareff0p1_lognormal.pickle', includedust=True,tau_d=dusttauvalue,dustparamsfile='dust_wolff_pangajello_reff1p5_vareff0p5_lognormal.pickle')

##################
####Dust and clouds (varying tau_d and cloud position)
##################
dusttaulabel_list=np.array(['0.1', '1', '10']) #list of total dust optical depths (500 nm)
dusttauvalue_list=np.array([0.1, 1., 10.]) #list of total dust optical depths (500 nm)

cloudpositionslabel_list=np.array(['0.5', '20.5', '40.5']) #list of cloud altitudes (km)
cloudindexvalue_list=np.array([63, 43, 23]) #indices corresponding to those cloud altitudes

for ind_dust in range(0, len(dusttaulabel_list)):
	for ind_cloud in range(0, len(cloudpositionslabel_list)):
		dusttaulabel=dusttaulabel_list[ind_dust]
		dusttauvalue=dusttauvalue_list[ind_dust]
		
		cloudpositionslabel=cloudpositionslabel_list[ind_cloud]
		cloudindexvalue=cloudindexvalue_list[ind_cloud]
		
		rt.uv_radtrans(z_upper_limit=64.e5, z_step=1.e5, inputatmofilelabel='colddrymars_0.02bar_250K', outputfilelabel='z=0_A=desert_noTD_DS_dustod='+dusttaulabel+'_co2cloudod=100_z='+cloudpositionslabel+'_reff=10', inputspectrafile='general_youngsun_mars_spectral_input.dat',TDXC=False, DeltaScaling=True, SZA_deg=0., albedoflag='nonuniformalbedo',nonuniformalbedo=np.array([0.,0.,1.,0.,0.]), includeh2ocloud=False,includeco2cloud=True,co2cloudlayerinds=np.array([cloudindexvalue]), co2cloudlayerods=np.array([100.]),co2iceparamsfile='cloud_co2_reff10_vareff0p1_lognormal.pickle', includedust=True,tau_d=dusttauvalue,dustparamsfile='dust_wolff_pangajello_reff1p5_vareff0p5_lognormal.pickle')
