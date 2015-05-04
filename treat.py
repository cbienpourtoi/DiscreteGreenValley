#! /usr/local/bin python
# -*- coding: iso-8859-1 -*-

""" 
- Gets R images from sdss 
- Imports the image you chose from Galex 
(branch for galex auto-import currently abandonned)
- Convolves sdss images to Galex images
- Projects/Register sdss on Galex field
- and then plays with the data obtained"""


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
import bz2

############
### Init ###
############

# Temporary files directory
tmp = "tmp/"
if not os.path.exists(tmp):
	os.mkdir(tmp)

# Name for the llow resolution sdss file
sdss_lowres_file = tmp+'sdss_lowres.fits'


# Initial files (sdss for R-band, galex for NUV band)
#GalexFile = "galex/MAST_2015-04-29T2046/GALEX/6371021344036880384/AIS_3_sg33-nd-int.fits"
#SDSSFile = "sdss/frame-r-004682-4-0092.fits"




######################
# Galaxy coordinates #
######################
# galex data from:
# http://galex.stsci.edu/GR6/?page=mastform

galaxies = Table.read("galaxy_list.txt", format="ascii")

for galaxy in galaxies:

	name = galaxy["Name"]
	RA = galaxy["RA"]
	Dec = galaxy["Dec"]
	GalexFile = galaxy["GalexFilePath"] 
	
	"""
	# NGC6125 E
	RA = 244.7958
	Dec = 57.9842

	# NGC4185 S
	#RA = 183.3417
	#Dec = 28.5103
	"""


	###################################################
	# Downloads directly from GALEX using coordinates:#
	###################################################
	""" Does not work (yet?). Instead, just download the image yourself.

	casjobs = "java -jar ~/sandbox/CasJobsCL/casjobs.jar" 

	galex_table_file = tmp+"galextable.csv"

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

	urllib.urlretrieve(galeximageurl, tmp+"GalexNUV.fits.gz")
	urllib.urlcleanup()

	#urllib.urlretrieve("http://galex.stsci.edu/data/GR6/pipe/02-vsn/50112-AIS_112/d/00-visits/0001-img/07-try/AIS_112_0001_sg68-nd-int.fits.gz", tmp+"GalexNUV.fits.gz")


	sys.exit()
	"""



	###################################################
	# Downloads directly from SDSS using coordinates: #
	###################################################

	sdss_table_file = tmp+"sdsstable.csv"
	urllib.urlretrieve("http://skyserver.sdss.org/dr10/en/tools/search/x_radial.aspx?ra="+str(RA)+"&dec="+str(Dec)+"&radius=0.2&format=csv&limit=20", sdss_table_file)
	urllib.urlcleanup()
	sdss_data = ascii.read(sdss_table_file, format='csv', header_start=1)

	run = sdss_data['run'][0] #4862
	rerun = sdss_data['rerun'][0] #301
	camcol = sdss_data['camcol'][0] #4
	field = sdss_data['field'][0] #92

	run2 = ("%6i" % int(run)).replace(' ', '0')
	field2 = ("%4i" % int(field)).replace(' ', '0')
	
	SDSSFile = tmp+"sdssR.fits.bz2"
	urllib.urlretrieve("http://dr10.sdss3.org/sas/dr10/boss/photoObj/frames/"+rerun+"/"+run+"/"+camcol+"/frame-r-"+run2+"-"+camcol+"-"+field2+".fits.bz2", SDSSFile)
	urllib.urlcleanup()



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
	sdss_hdulist = fits.open(bz2.BZ2File(SDSSFile))
	sdss_lowres = convolve_fft(sdss_hdulist[0].data, psf)

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
	sdssheader_tmp = tmp+"sdsstmp.hdr"

	# Projected Galex image on the header system of the sdss image
	galex_projected = tmp+"galex_projected.fits"

	# Creates the header template
	os.system("mGetHdr "+SDSSFile+" "+ sdssheader_tmp)

	# Does the projection
	os.system("mProject "+GalexFile+" "+ galex_projected+" "+sdssheader_tmp)



	################################
	######    RGB Image     ########
	################################

	sdss_lr_image = (fits.open(sdss_lowres_file))[0].data
	galex_proj_image = (fits.open(galex_projected))[0].data
	last_color = sdss_lr_image * 0.

	sdss_lr_image = sdss_lr_image / np.nansum(sdss_lr_image)
	galex_proj_image = galex_proj_image / np.nansum(galex_proj_image)

	rgbArray = np.zeros((sdss_lr_image.shape[0],sdss_lr_image.shape[1],3), 'uint8')
	rgbArray[..., 0] = sdss_lr_image * 25
	rgbArray[..., 1] = galex_proj_image * 250
	rgbArray[..., 2] = last_color

	img = Image.fromarray(rgbArray)

	img.save('myimg.jpeg')












################
#### Drafts ####
################


"""

aplpy.make_rgb_cube([sdss_lowres_file, sdss_lowres_file, galex_projected], 'rgb.fits')

aplpy.make_rgb_image('rgb.fits', 'rgb.png', vmin_r=20, vmax_r=400, vmin_g=0, vmax_g=150, vmin_b=-2,vmax_b=50)

fig = aplpy.FITSFigure('rgb.fits')

fig.show_rgb('rgb.png')

# Add contours
fig.show_contour('sc.fits', cmap='gist_heat', levels=[0.2,0.4,0.6,0.8,1.0])

# Overlay a grid
fig.add_grid()
fig.grid.set_alpha(0.5)

# Save image
fig.save('plot.png')
                     
"""
"""
img = np.zeros((sdss_lr_image.shape[0], sdss_lr_image.shape[1], 3), dtype=float)
img[:,:,0] = img_scale.sqrt(sdss_lr_image, scale_min=0, scale_max=10000)
img[:,:,1] = img_scale.sqrt(galex_proj_image, scale_min=0, scale_max=10000)
img[:,:,2] = img_scale.sqrt(last_color, scale_min=0, scale_max=10000)
"""



"""
RA = 244.7958
Dec = 57.9842

GalexFile = "galex/MAST_2015-04-29T2046/GALEX/6371021344036880384/AIS_3_sg33-nd-int.fits"
outgalex = "galex_stamp_l.fits"

SDSSFile = "sdss/frame-r-004682-4-0092.fits"
outsdss = "sdss_stamp_l.fits"

cutout.cutout(GalexFile, RA, Dec, 0.05, 0.05, units='wcs', outfile=outgalex, coordsys='celestial')
cutout.cutout(SDSSFile, RA, Dec, 0.04, 0.04, units='wcs', outfile=outsdss, coordsys='celestial')


images_to_align = sorted(glob.glob(outgalex))
ref_image = outsdss

images_to_align = sorted(glob.glob(outsdss))
ref_image = outgalex

#images_to_align = sorted(glob.glob("sdss/frame-r-004682-4-0092.fits"))
#ref_image = "galex/MAST_2015-04-29T2046/GALEX/6371021344036880384/AIS_3_sg33-nd-int.fits"

identifications = alipy.ident.run(ref_image, images_to_align, visu=False, hdu=1)
# That's it !
# Put visu=True to get visualizations in form of png files (nice but much slower)
# On multi-extension data, you will want to specify the hdu (see API doc).

# The output is a list of Identification objects, which contain the transforms :
for id in identifications: # list of the same length as images_to_align.
        if id.ok == True: # i.e., if it worked

                print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
                # id.trans is a alipy.star.SimpleTransform object. Instead of printing it out as a string,
                # you can directly access its parameters :
                #print id.trans.v # the raw data, [r*cos(theta)  r*sin(theta)  r*shift_x  r*shift_y]
                #print id.trans.matrixform()
                #print id.trans.inverse() # this returns a new SimpleTransform object

        else:
                print "%20s : no transformation found !" % (id.ukn.name)


# Minimal example of how to align images :

outputshape = alipy.align.shape(ref_image)
# This is simply a tuple (width, height)... you could specify any other shape.

for id in identifications:
        if id.ok == True:

                # Variant 1, using only scipy and the simple affine transorm :
                alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, makepng=True)

                # Variant 2, using geomap/gregister, correcting also for distortions :
                #alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars, shape=outputshape, makepng=False)
                # id.uknmatchstars and id.refmatchstars are simply lists of corresponding Star objects.

                # By default, the aligned images are written into a directory "alipy_out".






"""
"""
import aplpy

fileR = 'sdss/frame-r-004682-4-0092.fits'
#fileNUV = "galex/MAST_2015-04-29T2046/GALEX/6371021344036880384/AIS_3_sg33-fd-int.fits.gz"




f = aplpy.FITSFigure(fileR)
#f.show_grayscale()


"""

"""

hdulist = fits.open(fileR)
R = hdulist[1].data
if hdulist[0].header["FILTER"] != 'r': print "Error: wrong filter for SDSS"
hdulist.close()

hdulist = fits.open(fileR)

imgplot = plt.imshow(hdulist[0].data)
#imgplot.set_cmap('binary_r')

#plt.show()
#hdulist.close()
"""
