#! /usr/local/bin python
# -*- coding: iso-8859-1 -*-

""" 
Inputs:
- The list of objects name, coordinates, and galex NUV image path: galaxy_list.txt
- A galex NUV image

Process:
- Gets R images from sdss corresponding to the object
- (branch for galex auto-import currently abandonned)
- Convolves sdss images to Galex PSF
- Projects/Register Galex image on SDSS field
- and then plays with the data obtained

Outputs:

"""




# This code requires Montage installed on your system:
# http://montage.ipac.caltech.edu

# Only if trying to get data from Galex 
# (does not work yet, and maybe never will)
# it requires casjobs: 
# install it from http://galex.stsci.edu/casjobs/casjobscl.aspx
# casjobs = "java -jar ~/sandbox/CasJobsCL/casjobs.jar" 


#################
### Packages ####
#################
#import matplotlib
#matplotlib.use('Agg')
import numpy as np
from astropy.io import fits
import pyfits
import sys
#import glob
import os
from astropy.convolution import Gaussian2DKernel, convolve_fft
#import aplpy
from PIL import Image #If you want to make RGB image
import urllib
from astropy.table import Table
from astropy.io import ascii
#import bz2
import subprocess
from scipy.spatial.distance import cdist


#######################
# Galaxy initial data #
#######################

# Download galex data from:
# http://galex.stsci.edu/GR6/?page=mastform

galaxies = Table.read("galaxy_list.txt", format="ascii")

for galaxy in galaxies:

	print galaxy

	name = galaxy["Name"]
	RA = galaxy["RA"]
	Dec = galaxy["Dec"]
	GalexFile = galaxy["GalexFilePath"] 

	# Output files directory
	out = "out/"
	if not os.path.exists(out):
		os.mkdir(out)
	out = out+name+"/"
	if not os.path.exists(out):
		os.mkdir(out)



	###################################################
	# Downloads directly from GALEX using coordinates:#
	###################################################
	""" Does not work (yet?). Instead, just download the image yourself.

	casjobs = "java -jar ~/sandbox/CasJobsCL/casjobs.jar" 

	galex_table_file = out+"galextable.csv"

	casjobs_cmd = casjobs+" execute "
	casjobs_cmd = casjobs_cmd + "\"SELECT TOP 100 p.objid, dbo.fHasSpectrum(p.objid) as specObjID, n.distance as distance_arcmin, dbo.fIAUFromEq(p.ra,p.dec) as IAUName,p.ra,p.dec,p.fuv_mag, p.nuv_mag,p.fuv_flux, p.nuv_flux,p.e_bv,p.isThereSpectrum,p.fuv_fwhm_world,p.nuv_fwhm_world,p.vsn,p.prod,p.tilenum,p.try,p.img,p.band,p.id,p.subvisit,p.leg,p.ow,p.type,p.htmID FROM PhotoObjAll as p LEFT OUTER JOIN SpecObjAll as s on p.objid = s.objid, dbo.fGetNearbyObjEq("
	casjobs_cmd = casjobs_cmd + str(RA)+" , "+str(Dec)
	casjobs_cmd = casjobs_cmd +" , 0.5) as n WHERE p.band=1 AND p.objID = n.objID ORDER BY n.distance ASC,p.nuv_flux ASC ,p.band ASC ,p.e_bv ASC\""
	casjobs_cmd = casjobs_cmd + ">"+galex_table_file
	os.system(casjobs_cmd)

	galex_data = ascii.read(galex_table_file, format='csv', header_start=1)

	survey = "AIS"
	galex_filter = "nd" # nd for NUV, fd for FUV
	vsn = galex_data["[vsn]:Integer"][0]
	vsn2 = ("%2i" % int(vsn)).replace(' ','0') + "-vsn"
	tile = str(galex_data["[tilenum]:Integer"][0])
	tile2 = tile +"-"+survey+"_"+tile[-3:]
	tryy = galex_data["[try]:Integer"][0]
	tryy2 = ("%2i" % int(tryy)).replace(' ','0') + "-try"
	subvisit = str(galex_data["[subvisit]:Integer"][0])

	imagename = survey+"_"+tile[-3:]+"_0001_sg"+subvisit+"-"+galex_filter+"-int.fits.gz"

	galeximageurl = "http://galex.stsci.edu/data/GR6/pipe/"+vsn2+"/"+tile2+"/d/00-visits/0001-img/"+tryy2+"/"+imagename

	urllib.urlretrieve(galeximageurl, out+"GalexNUV.fits.gz")
	urllib.urlcleanup()

	#urllib.urlretrieve("http://galex.stsci.edu/data/GR6/pipe/02-vsn/50112-AIS_112/d/00-visits/0001-img/07-try/AIS_112_0001_sg68-nd-int.fits.gz", out+"GalexNUV.fits.gz")


	sys.exit()
	"""



	###################################################
	# Downloads directly from SDSS using coordinates: #
	###################################################

	# TODO:
	# Make a mosaic !
	# crap, automatic mosaic does work:
	# urllib.urlretrieve("http://data.sdss3.org/mosaic-server/mosaic?onlyprimary=False&pixelscale=0.396&ra=155.86248&filters=r&dec=19.89849&size=0.3", "test.tar")
	# gets a 403

	sdss_table_file = out+"sdsstable.csv"
	urllib.urlretrieve("http://skyserver.sdss.org/dr10/en/tools/search/x_radial.aspx?ra="+str(RA)+"&dec="+str(Dec)+"&radius=0.2&format=csv&limit=20", sdss_table_file)
	urllib.urlcleanup()
	sdss_data = ascii.read(sdss_table_file, format='csv', header_start=1)

	run = sdss_data['run'][0] #4862
	rerun = sdss_data['rerun'][0] #301
	camcol = sdss_data['camcol'][0] #4
	field = sdss_data['field'][0] #92

	run2 = ("%6i" % int(run)).replace(' ', '0')
	field2 = ("%4i" % int(field)).replace(' ', '0')
	
	SDSSFile = out+"sdssR.fits"
	urllib.urlretrieve("http://dr10.sdss3.org/sas/dr10/boss/photoObj/frames/"+rerun+"/"+run+"/"+camcol+"/frame-r-"+run2+"-"+camcol+"-"+field2+".fits.bz2", SDSSFile+".bz2")
	urllib.urlcleanup()
	try:
		subprocess.call(["bunzip2",SDSSFile+".bz2"])
	except OSError:
		print "Error in unzipping the sdss file: "+SDSSFile+".bz2"
		raise


	################################
	###### Convolution part ########
	################################

	# PSF to convolve SDSS to Galex resolution
	galex_psf_arcsec = 5. #arcsec FWHM
	galex_psf_deg = galex_psf_arcsec / 360.

	# Infos from Galex file
	galex_hdulist = fits.open(GalexFile)
	galex_cdelt1 = galex_hdulist[0].header["CDELT1"]
	galex_cdelt2 = galex_hdulist[0].header["CDELT2"]
	galex_hdulist.close()

	# Sanity check
	if np.abs(galex_cdelt1) != np.abs(galex_cdelt2):
		print "CDELT 1 and 2 in galex are different. That could mess the psf (or not, but I stop anyway)"
		sys.exit()

	# Create psf
	psf = Gaussian2DKernel(stddev=galex_psf_deg/np.abs(galex_cdelt1)/2.355)


	# Opens SDSS and convolves it to psf
	sdss_hdulist = fits.open(SDSSFile)
	sdss_lowres = convolve_fft(sdss_hdulist[0].data, psf)

	# Name for the low resolution sdss file
	sdss_lowres_file = out+'sdss_lowres.fits'

	# Saves the low resolution sdss
	hdu = fits.PrimaryHDU(sdss_lowres, header = sdss_hdulist[0].header)
	hdulist_sdss_lowres = fits.HDUList([hdu])
	hdulist_sdss_lowres.writeto(sdss_lowres_file, clobber=True)
	sdss_hdulist.close()
	hdulist_sdss_lowres.close()


	################################
	###### Projection part  ########
	################################

	# File in which we will create the header template requested by Montage 
	#(http://montage.ipac.caltech.edu/docs/mProject.html)
	sdssheader_tmp = out+"sdsstmp.hdr"

	# Projected Galex image on the header system of the sdss image
	galex_projected = out+"galex_projected.fits"

	# Creates the header template
	subprocess.call(["mGetHdr", SDSSFile, sdssheader_tmp])

	# Does the projection
	subprocess.call(["mProject", GalexFile, galex_projected, sdssheader_tmp])



	################################
	######    RGB Image     ########
	################################

	# Make a cutout (thumbnails) of the images
	tRfile = out+"thumbR.fits"
	tNUVfile = out+"thumbNUV.fits"

	subprocess.call(["mSubimage", sdss_lowres_file, tRfile, str(RA), str(Dec), "0.1"])
	subprocess.call(["mSubimage", galex_projected, tNUVfile, str(RA), str(Dec), "0.1"])


	thumbR = (fits.open(tRfile))[0].data
	thumbNUV = (fits.open(tNUVfile))[0].data
	
	#sdss_lr_image = sdss_lr_image / np.nansum(sdss_lr_image)
	#galex_proj_image = galex_proj_image / np.nansum(galex_proj_image)

	rgbArray = np.zeros((thumbR.shape[0],thumbR.shape[1],3), 'uint8')
	rgbArray[..., 0] = thumbR / np.nanmax(thumbR)*255.
	#rgbArray[..., 1] = thumbR * 0. # Green is already 0
	rgbArray[..., 2] = thumbNUV / np.nanmax(thumbR)*255.*15.

	img = Image.fromarray(rgbArray)

	img.save(out+name+'.png')
	
	#sys.exit()

