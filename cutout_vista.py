import numpy as np
import os, sys
import astropy.io.fits as fits
import astropy.io.fits.compression
import matplotlib.pyplot as plt
import cStringIO as cS

def makestamp(ra, dec, filename, chipno, size=60.0, cmap='gray', annotate=None, ctype='png', db = "vista", onlytiles = True, conf=False):
	"""
	Return a postage stamp as png binary data

	ra,dec in decimal degrees
	size in arcsec
	cmap is colormap: http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps?action=AttachFile&do=get&target=colormaps3.png
	"""
	import xmlrpclib
	import ConfigParser

	config_file = "srpylib.cfg"
	server = config.get("vista cutouts", server)
	
	xms = xmlrpclib.ServerProxy(server)
	dd = {'filename': filename, 'hdu': chipno, 'ra': ra, 'dec': dec, 'size': size, 'cmap': cmap, 'ctype': ctype, "db": db}
	# By default catalogue sources are overplotted. To turn that off use instead of the above line this
	# # dd = {'filename': filename, 'hdu': chipno, 'ra': ra, 'dec': dec, 'size': size, 'cmap': cmap, 'nocat': True}
	if annotate:
		dd['annotate']=annotate
	if conf:
		dd['conf']='1'
	if "vst" in dd["filename"] and onlytiles and dd["filename"][-6:-4] == "st":
		cutout = xms.cutout(dd)
		return cutout.data, dd["filename"], dd["hdu"]
	elif "vista" in dd["filename"]:
		cutout = xms.cutout(dd)
		return cutout.data, dd["filename"], dd["hdu"]
	elif not onlytiles:
		cutout = xms.cutout(dd)
		return cutout.data, dd["filename"], dd["hdu"]
	else:
		try:
			cutout = xms.cutout(dd)
		except xmlrpclib.Fault as e:
			cutout = np.zeros(shape = (10,10), dtype = "int8")
			return cutout, dd["filename"], dd["hdu"]
		return cutout.data, dd["filename"], dd["hdu"]


def querydb(ra,dec, size = 60.0, prefix='img',db="vista", onlytiles=True, ctype='png', allfiles=False, conf = False):
	"""
	Query database for objects around position and save postage stamp
	"""
	import xmlrpclib
	import ConfigParser

	config_file = "srpylib.cfg"
	server = config.get("vista cutouts", server)

	xms = xmlrpclib.ServerProxy(server)
	dd={'ra': ra, 'dec': dec, 'db': db, 'onlytiles': onlytiles, 'allfiles': allfiles}
	filenames = []
	apm_filenames = []
	chipnos = []
	n = 0
	for item in xms.queryradec(dd):
		# The line below produces the postage stamp and returns png data which needs to be saved into a file
		# # The annotate part is optional
		png, apm_file, chipno = makestamp(ra, dec, os.path.join(item['filepath'], item['filename']), item['chipno'], annotate='%s -- %s\n%s' % (db.upper(), item['filtername'], item['filename']), ctype=ctype, size = size, db = db, onlytiles = onlytiles, conf = conf)
		# Wihtout annotating the plot it would be just
	    # png = makestamp(ra, dec, os.path.join(item['filepath'], item['filename']), item['chipno'])
		output='/tmp/%s_%s_%s_%s%s.%s' % (prefix, item['filename'].replace('.fit',''), item['chipno'], item['filtername'], {True: '_conf', False: ''}[conf], ctype)
		#output = "/tmp/" + str(ra) + "_" + str(dec) + ".fits"
		if png is not None:
			open(output,'wb').write(png)
			filenames.append(output)
			chipnos.append(chipno)
			apm_filenames.append(apm_file)
	return filenames, apm_filenames, chipnos

