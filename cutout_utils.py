"""
Functions to make various cutouts and plots
"""

def MAD(l, med):
	import numpy as np
	return np.median(abs( l - med ))

def SED_plot(mags, magErrs, id):

	"""
	Makes plots of the magnitudes in different bands.
	Uses id as title
	Assumes VHS and WISE are in Vega and converts them
	Uses conversions from WISE explanatory Supplement and CASU website
	Sets 99 DES objects to 25.0 with an error of 1.0
	Also sets WISE limits to 25.0 and nan magnitudes to 25.0
	Returns rather than draws the plot
	"""

	import matplotlib.pyplot as plt
	import numpy as np
	import pylab

	mids = np.array([4800.0, 6400.0, 7800.0, 9200.0, 9900.0, 13140.0, 16890.0, 22130.0, 34000.0, 46000.0])
	filters = np.array(["g", "r", "i", "z", "Y", "J", "H", "K", "W1", "W2"])
	convs = np.array([0.0,0.0,0.0,0.0,0.0,0.937,1.384,1.839,2.699,3.339, 5.174, 6.620])

	fig = plt.figure()
	ax = fig.add_subplot(111)

	ids = np.where( (mags <> mags) | (mags == 99.0) | (magErrs <> magErrs) )[0]

	magErrs[ids] = 1.0
	mags[ids] = 25.0
	ax.errorbar(mids, mags+convs, yerr = [magErrs, magErrs], xerr = [[0.0]*len(mags), [0.0]*len(mags)], ls = ".", ms = 0, color = "black")
	ax.plot(mids, mags+convs, "k.", label = "Actual values")
	ax.invert_yaxis()
	ax.set_title(id)
	pylab.xticks(mids, filters)

	return fig

def make_cutout(filename, RA, DEC, width, nhdu = 0, w_units = "arcsecs", verbose = False):

	"""
	Makes cutouts from a file passed as either a filename or hdulist
	Returns a new hdulist with updated WCS and the cutout as the data
	The specified width should be in arcsecs or pixels.
	It is the full width.
	"""

	from astropy import wcs
	import wcsutil
	import astropy.io.fits as fits
	import astropy.io.fits.compression
	import numpy as np
	import gc
	import matplotlib.pyplot as plt
	import warnings
	from astropy.utils.exceptions import AstropyWarning
	import copy

	im_status = "good"

	if not verbose:
		warnings.filterwarnings('ignore', category=UserWarning, append=True)
		warnings.simplefilter('ignore', category=AstropyWarning)

	if isinstance(filename, basestring):
		with fits.open(filename) as h:
			hdr = h[nhdu].header
			data = h[nhdu].data
	else:
		try:
			hdr = filename[nhdu].header
			#It is this line data = filename[nhdu].data
			data = filename[nhdu].data
		except RuntimeError:
			print "File Corrupt"
			return None, "bad"
	try:
		test = hdr["NAXIS1"]
	except KeyError:
		return None, "bad"

	#Figure out how big the square would be at the centre
	test_pix = [[hdr["NAXIS1"]/2, hdr["NAXIS2"]/2]]

	try:
		from astropy import wcs
		w = wcs.WCS(hdr, naxis = 2)

		if w_units == "arcsecs":
			width = width/3600.0/2.0
			w_coord = np.array([[RA - width/np.cos(np.deg2rad(DEC)), DEC - width], [RA + width/np.cos(np.deg2rad(DEC)), DEC + width]], np.float_)
			pix_coord = w.wcs_world2pix(w_coord, 1)
			coords = [int(pix_coord[0][0]), int(pix_coord[0][1]), int(pix_coord[1][0]), int(pix_coord[1][1])]
			test_coord = w.wcs_pix2world(test_pix, 1)
			test_edges = np.array([[test_coord[0][0], test_coord[0][1] - width], \
								[test_coord[0][0], test_coord[0][1] + width]], np.float_)
			test_pix = w.wcs_world2pix(test_edges, 1)
			[[tx1, ty1], [tx2, ty2]] = test_pix

		elif w_units == "pixels":
			w_coord = np.array([[RA, DEC]], np.float_)
			pix_coord = w.wcs_world2pix(w_coord, 1)
			coords = [int(pix_coord[0][0]+width/2.0), int(pix_coord[0][1]-width/2.0), int(pix_coord[0][0]-width/2.0), int(pix_coord[0][1]+width/2.0)]
			#Create a background of zeros
			blank = np.zeros((width+1, width+1))
	except (UnboundLocalError, wcs.wcs.InconsistentAxisTypesError) as e:
		try:
			wcs = wcsutil.WCS(hdr)
		except KeyError:
			hdr = h[nhdu-1].header
			wcs = wcsutil.WCS(hdr)

		width = width/3600.0/2.0
		test_coord = wcs.image2sky(test_pix[0][0], test_pix[0][1])
		tx1, ty1 = wcs.sky2image(test_coord[0], test_coord[1] - width)
		tx2, ty2 = wcs.sky2image(test_coord[0], test_coord[1] + width)

		x1, y1 = wcs.sky2image(RA - width/np.cos(np.deg2rad(DEC)), DEC - width)
		x2, y2 = wcs.sky2image(RA + width/np.cos(np.deg2rad(DEC)), DEC + width)
		coords = [int(x1), int(y1), int(x2), int(y2)]

	if tx1 > tx2:
		x_width = tx1 - tx2
	else: 
		x_width = tx1 - tx2
	if ty1 > ty2:
		y_width = ty1 - ty2
	else:
		y_width = ty2 - ty1
	if x_width > y_width:
		b_width = x_width
	else:
		b_width = y_width
	blank = np.zeros((int(np.ceil(b_width))+1, int(np.ceil(b_width))+1))

	#coords = [x1, y1, x2, y2]
	coords_clean = coords + []
	for (n,p) in enumerate(coords):
		if p < 0.0: 
			coords_clean[n] = 0.0
			im_status = "bad"

		elif (p > hdr["NAXIS1"] and n % 2 == 0):
			coords_clean[n] = hdr["NAXIS1"]
			im_status = "bad"
		elif (p > hdr["NAXIS2"] and n % 2 <> 0):
			coords_clean[n] = hdr["NAXIS2"]
			im_status = "bad"

	if len(data.shape) > 2:
		data = data[0,0,:,:]
	if w_units == "pixels":
		im = data[coords_clean[1]:coords_clean[3], coords_clean[2]:coords_clean[0]]
	if w_units == "arcsecs":
		im = data[int(coords_clean[1]):int(coords_clean[3]), int(coords_clean[2]):int(coords_clean[0])]

	if coords[0] > hdr["NAXIS1"]:
		by2 = hdr["NAXIS1"]-coords[2]
		im_status = "bad"
	else:
		by2 = im.shape[1]

	if coords[2] < 0:
		by1 = np.fabs(coords[2])
		im_status = "bad"
		by2 += by1
	else:
		by1 = 0

	if coords[3] > hdr["NAXIS2"]:
		bx2 = hdr["NAXIS2"]-coords[1]
		im_status = "bad"
	else:
		bx2 = im.shape[0]

	if coords[1] < 0:
		bx1 = np.fabs(coords[1])
		im_status = "bad"
		bx2 += bx1
	else:
		bx1 = 0

	if (coords_clean[1] == coords_clean[3]) or (coords_clean[0] == coords_clean[2]):
		im = blank
		im_status = "bad"

	else:
		try:
			blank[int(bx1):int(bx2), int(by1):int(by2)] = im
		except ValueError:
			#blank = np.zeros((max(im.shape), max(im.shape)))
			blank = np.zeros((max([bx2, by2]), max([bx2, by2])))
			blank[int(bx1):int(bx2), int(by1):int(by2)] = im
		im = blank

	phdu = fits.PrimaryHDU()
	h2 = copy.deepcopy(hdr)
	h2["CRPIX1"] = hdr["CRPIX1"]-coords[2]
	h2["CRPIX2"] = hdr["CRPIX2"]-coords[1]
	imhdu = fits.ImageHDU(header = h2, data = im)
	hdulist = fits.HDUList([phdu, imhdu])

	return hdulist, im_status

def draw_crosshairs(ax, w, ra, dec):

	"""
	Draws crosshairs on an image axis
	Requires the axis to draw on, ax the wcs information w as an astropy wcs thing
	eg w = wcs.WCS(hdr, naxis = 2) and the ra and dec of object
	Also puts on N and E arrows.
	Not extensivly tested may require tweaking
	"""

	import wcsutil
	from astropy import wcs
	import numpy as np
	import matplotlib.pyplot as plt

	cen_world = [[ra, dec]]
	try:
		cen_coord = w.wcs_world2pix(cen_world, 1)
	except AttributeError:
		cen_coord = [w.sky2image(cen_world[0][0], cen_world[0][1])]
	cen_x = cen_coord[0][0]
	cen_y = np.floor(cen_coord[0][1])

	ax.autoscale(False)
	ax.plot([cen_x, cen_x], [cen_y + 7, cen_y + 27], "k", lw = 3)
	ax.plot([cen_x, cen_x], [cen_y - 7, cen_y - 27], "k", lw = 3)
	ax.plot([cen_x + 7, cen_x + 27], [cen_y, cen_y], "k", lw = 3)
	ax.plot([cen_x - 7, cen_x - 27], [cen_y, cen_y], "k", lw = 3)

	eline = [[ra + 5.0/3600.0, dec - 10.0/3600.0], [ra + 10.0/3600.0, dec - 10.0/3600.0], [ra + 10.0/3600.0, dec - 12.0/3600.0]]
	try:
		eline = w.wcs_world2pix(eline, 1)
	except AttributeError:
		eline1 = w.sky2image(eline[0][0], eline[0][1])
		eline2 = w.sky2image(eline[1][0], eline[1][1])
		eline3 = w.sky2image(eline[2][0], eline[2][1])
		eline = [eline1, eline2, eline3]
	ax.arrow(eline[0][0], eline[0][1], eline[1][0]-eline[0][0], eline[1][1]-eline[0][1], fc = "k", ec = "k", head_width = 5.0, head_length = 5.0, lw = 2)
	ax.annotate("E", xy = (eline[2][0], eline[2][1]), size = 12, weight = "bold")

	nline = [[ra + 5.0/3600.0, dec - 10.0/3600.0], [ra + 5.0/3600.0, dec - 5.0/3600.0], [ra + 9.0/3600.0, dec - 5.0/3600.0]]
	try:
		nline = w.wcs_world2pix(nline, 1)
	except AttributeError:
		nline1 = w.sky2image(nline[0][0], nline[0][1])
		nline2 = w.sky2image(nline[1][0], nline[1][1])
		nline3 = w.sky2image(nline[2][0], nline[2][1])
		nline = [nline1, nline2, nline3]
	ax.arrow(nline[0][0], nline[0][1], nline[1][0]-nline[0][0], nline[1][1]-nline[0][1], fc = "k", ec = "k", head_width = 5.0, head_length = 5.0, lw = 2)

	ax.annotate("N", xy = (nline[2][0]-2.0, nline[2][1]), size = 12, weight = "bold")

	return ax

def find_north(RA, DEC, w):

	import wcsutil
	from astropy import wcs

	nline = [[RA + 5.0/3600.0, DEC - 10.0/3600.0], [RA + 5.0/3600.0, DEC - 5.0/3600.0], [RA + 9.0/3600.0, DEC - 5.0/3600.0]]
	try:
		nline = w.wcs_world2pix(nline, 1)
	except AttributeError:
		nline1 = w.sky2image(nline[0][0], nline[0][1])
		nline2 = w.sky2image(nline[1][0], nline[1][1])
		nline3 = w.sky2image(nline[2][0], nline[2][1])
		nline = [nline1, nline2, nline3]
	if int(nline[1][0])-int(nline[0][0]) > 0:
		north = "right"
	elif int(nline[1][1])-int(nline[0][1]) < 0:
		north = "up"
	elif int(nline[1][1])-int(nline[0][1]) > 0:
		north = "down"
	elif int(nline[1][0])-int(nline[0][0]) < 0:
		north = "left"

	return north

def PAngle(ra1, dec1, ra2, dec2, frame = "icrs", verbose = True, debug = False, test=False):

	"""
	Calculates the posistion angle between two objects
	"""
	from astropy import units as u
	from astropy.coordinates import SkyCoord
	from astropy.coordinates import ICRS, Galactic, FK4, FK5
	from astropy.coordinates import Angle, Latitude, Longitude
	import numpy
	import astropy.coordinates as coord

	radec1 = SkyCoord(ra1, dec1, unit = (u.degree, u.degree), frame = frame)
	radec2 = SkyCoord(ra2, dec2, unit = (u.degree, u.degree), frame = frame)

	sep = radec1.separation(radec2).arcsecond

	PA = radec1.position_angle(radec2).degree

	return sep, PA

def make_DES_finding_chart(RA, DEC, mag, id, tile, run, outdir, filename, release = "Y1A1", band = "z", bright_star = False, width = 4.0*60.0, invert = False, invert_east = False):

	"""
	If Y2 then put tile and run as ""
	"""

	import DES_utils as Du
	import numpy as np
	import matplotlib.pyplot as plt
	import wcsutil
	from astropy import units as u
	from astropy.coordinates import SkyCoord
	import sqlutil
	from matplotlib.patches import Ellipse, Circle
	import match_lists
	import matplotlib.patheffects as path_effects
	import plotid

	if invert:
		cmap = "gray_r"
		font_col = "white"
		line_col = "black"

	else:
		cmap = "gray"
		font_col = "lawngreen"
		line_col = "white"

	fig = plt.figure(figsize = (8.27, 11.69))
	ax = fig.add_axes([0.0, 0.3, 1.0, 0.6])
	ax1 = fig.add_axes([0.05, 0.05, 0.90, 0.20])

	ax.axes.get_xaxis().set_visible(False)
	ax.axes.get_yaxis().set_visible(False)
	ax1.axes.get_xaxis().set_visible(False)
	ax1.axes.get_yaxis().set_visible(False)
	for i in ax1.spines.itervalues():
		i.set_linewidth(0.0)
	for i in ax.spines.itervalues():
		i.set_linewidth(2.0)

	if "Y2" in release:
		files, bands = Du.DES_Y2_image_file(RA, DEC, bands = [band], download = True)
		nhdu = 0
		filename = files[0]

	else:
		filename = "/data/desardata/" + release + "/" + tile + "/" + tile + "_" + band + ".fits.fz"
		nhdu = 1

	hdulist, im_status = make_cutout(filename, RA, DEC, width, nhdu = nhdu, w_units = "arcsecs")

	if (hdulist is None or im_status == "bad") and "Y2" in release:
		n = 0
		while n < len(files[1:]):
			hdulist, im_status = make_cutout(files[n], RA, DEC, width, nhdu = nhdu, w_units = "arcsecs")
			if hdulist is not None and im_status == "good":
				n = len(files)
			n += 1

	hdr = hdulist[1].header
	image = hdulist[1].data

	data = image.flatten()
	med = np.median(data)
	sigma_MAD = 1.4826 * MAD(data, med)
	vmax = med + 5.0*sigma_MAD
	vmin = med - 2.0*sigma_MAD
	try:
		from astropy import wcs
		w = wcs.WCS(hdr, naxis = 2)
	except (UnboundLocalError, wcs.wcs.InconsistentAxisTypesError) as e:
		try:
			w = wcsutil.WCS(hdr)
		except KeyError:
			hdr = h[nhdu-1].header
			w = wcsutil.WCS(hdr)

	north = find_north(RA, DEC, w)

	if north == "up":
		rot_mat = np.matrix( ((1,1), (1,1)) )
		ax.imshow(image, vmin = vmin, vmax = vmax, picker = True, interpolation = "nearest", cmap = cmap)
	elif north == "down":
		image = np.flipud(image)
		#PLaceholder as just want to subtract(add) centre y from each y.
		rot_mat = np.matrix( ((1,1), (1,1)) )
		ax.imshow(image, vmin = vmin, vmax = vmax, picker = True, interpolation = "nearest", cmap = cmap)
	elif north == "right":
		image = np.rot90(image, 1)
		rot_mat = np.matrix( ((0,-1), (1,0)) )
		ax.imshow(zip(*image)[::-1], vmin = vmin, vmax = vmax, picker = True, interpolation = "nearest", cmap = cmap)
	elif north == "left":
		image = np.rot90(image, 3)
		rot_mat = np.matrix( ((0,1), (-1,0)) )
		ax.imshow(image, vmin = vmin, vmax = vmax, picker = True, interpolation = "nearest", cmap = cmap)

	cen_world = [[RA, DEC]]
	try:
		cen_coord = w.wcs_world2pix(cen_world, 1)
	except AttributeError:
		cen_coord = [w.sky2image(cen_world[0][0], cen_world[0][1])]
	cen_x = cen_coord[0][0]
	cen_y = np.floor(cen_coord[0][1])

	ax.autoscale(False)
	ax.plot([cen_x, cen_x], [cen_y + 10, cen_y + 30], color = line_col, lw = 2)
	ax.plot([cen_x, cen_x], [cen_y - 10, cen_y - 30], color = line_col, lw = 2)
	ax.plot([cen_x + 10, cen_x + 30], [cen_y, cen_y], color = line_col, lw = 2)
	ax.plot([cen_x - 10, cen_x - 30], [cen_y, cen_y], color = line_col, lw = 2)

	#Add arrows
	if invert_east:
		ax.annotate("", xy=(0.15, 0.05), xytext=(0.05, 0.05), textcoords = "axes fraction", \
				color = line_col, fontsize = 15, xycoords = "axes fraction", \
				arrowprops=dict(fc=line_col, width = 2.0, frac = 0.3, ec = line_col))
		text = ax.annotate("E", xy = (0.17,0.04), xycoords = "axes fraction", \
				color = font_col, fontsize = 15)
		text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
		ax.annotate("", xy=(0.05, 0.15), xytext=(0.05, 0.05), textcoords = "axes fraction", \
				color = line_col, fontsize = 15, xycoords = "axes fraction", \
				arrowprops=dict(fc= line_col, width = 2.0, frac = 0.3, ec = line_col))
		text = ax.annotate("N", xy = (0.04,0.17), xycoords = "axes fraction", \
				color = font_col, fontsize = 15)
		text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

	if not invert_east:
		ax.annotate("", xy=(0.05, 0.15), xytext=(0.15, 0.15), textcoords = "axes fraction", \
				color = line_col, fontsize = 15, xycoords = "axes fraction", \
				arrowprops=dict(fc=line_col, width = 2.0, frac = 0.3, ec = line_col))
		text = ax.annotate("E", xy = (0.02,0.14), xycoords = "axes fraction", \
				color = font_col, fontsize = 15)
		text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
		ax.annotate("", xy=(0.15, 0.25), xytext=(0.15, 0.15), textcoords = "axes fraction", \
				color = line_col, fontsize = 15, xycoords = "axes fraction", \
				arrowprops=dict(fc= line_col, width = 2.0, frac = 0.3, ec = line_col))
		text = ax.annotate("N", xy = (0.14,0.27), xycoords = "axes fraction", \
				color = font_col, fontsize = 15)
		text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

	#Add box
	try:
		s_coords = [[RA + 30.0/3600.0/np.cos(np.deg2rad(DEC)), DEC + 30.0/3600.0], [RA - 30.0/3600.0/np.cos(np.deg2rad(DEC)), DEC - 30.0/3600.0]]
		[[raplus, decplus], [raminus, decminus]] = w.wcs_world2pix(s_coords, 1) 
	except AttributeError:
		raplus, decplus = w.sky2image(RA + 30.0/3600.0/np.cos(np.deg2rad(DEC)), DEC + 30.0/3600.0)
		raminus, decminus = w.sky2image(RA - 30.0/3600.0/np.cos(np.deg2rad(DEC)), DEC - 30.0/3600.0)

	ax.plot([raplus, raplus], [decminus, decplus], color = line_col, lw = 2)
	ax.plot([raminus, raminus], [decminus, decplus], color = line_col, lw = 2)
	ax.plot([raminus, raplus], [decminus, decminus], color = line_col, lw = 2)
	ax.plot([raminus, raplus], [decplus, decplus], color = line_col, lw = 2)

	if invert_east:
		text = ax.annotate('60"', xy = (raplus-20, cen_y), color = font_col, fontsize = 15)

	else:
		text = ax.annotate('60"', xy = (raplus-50, cen_y), color = font_col, fontsize = 15)
	text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

	#Add circles
	c_width = abs(raplus - raminus)
	ax.add_patch(Ellipse(xy = (cen_x, cen_y), width = c_width*2, height = c_width*2, angle = 0.0, fill = False, lw = 1, color = line_col)) 
	ax.add_patch(Ellipse(xy = (cen_x, cen_y), width = c_width*4, height = c_width*4, angle = 0.0, fill = False, lw = 1, color = line_col)) 

	#Add object details
	c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='icrs')
	coord_string = "%02.0f:%02.0f:%06.3f  %02.0f:%02.0f:%05.2f" % (c.ra.hms[0], c.ra.hms[1], c.ra.hms[2], c.dec.dms[0], np.fabs(c.dec.dms[1]), np.fabs(c.dec.dms[2]))

	ax1.text(0.01, 0.88, "Name        RA             Dec           " + band + " [AB]   PA       ROA", fontweight = "bold", family = "monospace")
	ax1.text(0.01, 0.73, "Object      " + coord_string + "   %0.2f " % mag, family = "monospace")

	radius = "%0.5f" % (60.0/3600.0/2.0)
	if release == "Y1A1":
		query = "select RA, DEC, MAG_PSF_" + band.upper() + " from des_y1a1.coadd_objects WHERE \
				q3c_radial_query(ra, dec, " + str(RA) + ", " + str(DEC) + \
				", " + radius + ")"

		ras_info, decs_info, mags_info = sqlutil.get(query, db = "wsdb", host = "cappc127", user = "sophiereed", password = "5ef23edc6")

	if release == "Y2Q1":
		query = "select RA, DEC, MAG_PSF_" + band.upper() + " from des_y2q1.objects WHERE \
				q3c_radial_query(ra, dec, " + str(RA) + ", " + str(DEC) + \
				", " + radius + ")"

		ras_info, decs_ino, mags_info = sqlutil.get(query, db = "wsdb", host = "cappc127", user = "sophiereed", password = "5ef23edc6")

	for (ra_info, dec_info, mag_info) in zip(ras_info, decs_info, mags_info):
		try:
			xy_info = w.wcs_world2pix([[ra_info, dec_info]], 1)
		except AttributeError:
			xy_info = [w.sky2image(ra_info, dec_info)]
		if north <> "down":
			xy_info = np.array((xy_info[0][0], xy_info[0][1]))*rot_mat
			xy_info = (xy_info.getA()[0][0], xy_info.getA()[0][1])
		else:
			y_info = xy_info[0][1]
			y_info = cen_y + (cen_y - y_info)
			xy_info = (xy_info[0][0], y_info)
		if invert_east:
			text1 = ax.annotate("%0.2f" % mag_info, xy = (xy_info[0] - 10.0, xy_info[1] - 0.0) , color = "k", fontweight = "light")
		else:
			text1 = ax.annotate("%0.2f" % mag_info, xy = (xy_info[0] + 10.0, xy_info[1] + 0.0) , color = "k", fontweight = "light")

	if bright_star:
		radius = str(width/2.0*np.sqrt(2.0)/3600.0)

		if release == "Y1A1":
			query = "select RA, DEC, MAG_PSF_" + band.upper() + " from des_y1a1.coadd_objects WHERE \
					q3c_radial_query(ra, dec, " + str(RA) + ", " + str(DEC) + \
					", " + radius + ") and MAG_PSF_" + band.upper() + "< 20 and MAG_PSF_" + band.upper() + "> 16"

			ras_bs, decs_bs, mags_bs = sqlutil.get(query, db = "wsdb", host = "cappc127", user = "sophiereed", password = "5ef23edc6")

		if release == "Y2Q1":
			print band.upper(), "band"
			query = "select RA, DEC, MAG_PSF_" + band.upper() + " from des_y2q1.objects WHERE \
					q3c_radial_query(ra, dec, " + str(RA) + ", " + str(DEC) + \
					", " + radius + ") and MAG_PSF_" + band.upper() + " < 19 and MAG_PSF_" + band.upper() + " > 16"

			ras_bs, decs_bs, mags_bs = sqlutil.get(query, db = "wsdb", host = "cappc127", user = "sophiereed", password = "5ef23edc6")

		"""
		pm_ra, pm_dec, ra_proMot, dec_proMot = sqlutil.get(pm_query, db = "wsdb", host = "cappc127", user = "sophiereed", password = "5ef23edc6")

		dists, inds = match_lists.match_lists(ras_bs, decs_bs, pm_ra, pm_dec, 2.0/3600.0)
		ids = np.where( (inds <> len(pm_ra)) )[0]
		print ids
		ras_bs = ras_bs[ids]
		decs_bs = decs_bs[ids]
		print ra_proMot[inds[ids]]
		"""

		dists = []
		PAs = []
		for (ra_bs, dec_bs) in zip(ras_bs, decs_bs):
			dist, PA = PAngle(RA, DEC, ra_bs, dec_bs)
			PAs.append(PA)
			dists.append(dist)

		info_bs = zip(ras_bs, decs_bs, mags_bs, dists, PAs)
		info_bs.sort(key = lambda tup: tup[3])
		info_bs = info_bs[1:]

		letters = ["A", "B", "C", "D", "E", "F", "G", "H"]

		m = 0
		y_pos = 0.58
		for (ra_bs, dec_bs, mag_bs, dist, PA) in info_bs:
			try:
				xy_bs = w.wcs_world2pix([[ra_bs, dec_bs]], 1)
			except AttributeError:
				xy_bs = [w.sky2image(ra_bs, dec_bs)]
			if north <> "down":
				xy_bs = np.array((xy_bs[0][0], xy_bs[0][1]))*rot_mat
				xy_bs = (xy_bs.getA()[0][0], xy_bs.getA()[0][1])
			else:
				y_bs = xy_bs[0][1]
				y_bs = cen_y + (cen_y - y_bs)
				xy_bs = (xy_bs[0][0], y_bs)
			ax.add_patch(Ellipse(xy = xy_bs, width = 20.0, height = 20.0, angle = 0.0, fill = False, lw = 1, color = font_col))
			if m < len(letters):
				print letters[m], ra_bs, dec_bs
				if invert_east:
					text = ax.annotate(letters[m], xy = (xy_bs[0] - 10.0, xy_bs[1] - 10.0) , color = font_col, fontweight = "bold")
					#text1 = ax.annotate("%0.2f" % mag_bs, xy = (xy_bs[0] - 20.0, xy_bs[1] - 10.0) , color = "k", fontweight = "light")
				else:
					text = ax.annotate(letters[m], xy = (xy_bs[0] + 15.0, xy_bs[1] + 15.0) , color = font_col, fontweight = "bold")
					#text1 = ax.annotate("%0.2f" % mag_bs, xy = (xy_bs[0] + 30.0, xy_bs[1] + 15.0) , color = "k", fontweight = "light")

				text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
				c_bs = SkyCoord(ra=ra_bs*u.degree, dec=dec_bs*u.degree, frame='icrs')
				coord_string_bs = "%02.0f:%02.0f:%06.3f  %02.0f:%02.0f:%05.2f" % \
				(c_bs.ra.hms[0], c_bs.ra.hms[1], c_bs.ra.hms[2], c_bs.dec.dms[0], np.fabs(c_bs.dec.dms[1]), np.fabs(c_bs.dec.dms[2]))

				ROA = PA - 90.0

				letter_info = letters[m] + "           " + coord_string_bs + "   %04.2f    %06.2f   %06.2f" % (mag_bs, PA, ROA)
				ax1.text(0.01, y_pos, letter_info, family = "monospace")
				m += 1
				y_pos -= 0.08


	if invert_east:
		ax.invert_xaxis()
	ax.set_title(filename + "    " + tile)
	plotid.plotid()
	name = " - VDES%02.0f%02.0f%02.0f%02.0f" % (c.ra.hms[0], c.ra.hms[1], c.dec.dms[0], np.fabs(c.dec.dms[1]))
	fig.suptitle(id + name, fontweight = "bold")
	#plt.show()
	if invert_east:
		out_name = outdir + "Finding_Chart_" + str(id) + "_" + band + "_%02.0f_WE.pdf" % width
	else:
		out_name = outdir + "Finding_Chart_" + str(id) + "_" + band + "_%02.0f_EW.pdf" % width
	plt.savefig(out_name, dpi = 300)
	plt.close()

def cutout_scale(im, num_min = 2.0, num_max = 5.0):

	"""
	Takes an image array and returns the vmin and vmax required to scale the image 
	between median + 5 * sigma MAD and median - 2 * sigma MAD
	"""

	import numpy as np

	data = im.flatten()
	try:
		med = np.median(data[data <> 0.0])
		sigma_MAD = 1.4826 *MAD(data[data <> 0.0], med)
	except IndexError:
		med = 0.0
		sigma_MAD = 0.0
	vmax = med + num_max * sigma_MAD
	vmin = med - num_min * sigma_MAD

	return vmin, vmax

def cutout_WISE(outdir, RA, DEC, id1, save = False, width = 30.0, small = False, release = "WISE_Images", verbose = True, return_val = False):

	"""
	Makes cutouts from the WISE images stored here.
	If the image doesn't exist downloads it.
	Width is the full width of the cutout and is in arcsecs
	Small removes the title and reduces the white space around the plot
	Missing data is given as a blank cutout
	"""

	import atpy
	import numpy as np
	import matplotlib.gridspec as gridspec
	import matplotlib.pyplot as plt
	import os
	import subprocess
	import astropy.io.fits as fits

	t = atpy.Table("/data/rgm/wise/wise_allwise_metadata_thin.fits")

	id = np.where( (RA < t["ra1"]) & (RA < t["ra4"]) & (RA > t["ra2"]) & (RA > t["ra3"]) & (DEC < t["dec3"]) & (DEC < t["dec4"]) & (DEC > t["dec2"]) & (DEC > t["dec1"]) )[0]
	datas = []

	n = 0
	while n < len(t):
		ras = sorted([t["ra1"][n], t["ra2"][n], t["ra3"][n], t["ra4"][n]])
		decs = sorted([t["dec1"][n], t["dec2"][n], t["dec3"][n], t["dec4"][n]])
		min_ra = ras[1]
		max_ra = ras[2]
		min_dec = decs[1]
		max_dec = decs[2]
		if max_ra > min_ra + 10.0:
			min_ra, max_ra = max_ra, min_ra
		if DEC > min_dec and DEC < max_dec and RA > min_ra and RA < max_ra:
			tile_id_unWISE = t["coadd_id"][n][0:8]
			f2 = tile_id_unWISE[0:3]
			wise_dir = "/data/wiseardata/" + release + "/p3am_cdd/" + f2 + "/" + tile_id_unWISE + "/"
			if verbose:
				print wise_dir
		n += 1

	if len(id) > 0:
		tile_id = t["coadd_id"][id[0]][0:13]
		f1 = tile_id[0:2]
		f2 = tile_id[0:4]
		bands = ["w1", "w2", "w3", "w4"]
		fig = plt.figure()
		gs = gridspec.GridSpec(1,4)
		s = 0
		for band in bands:
			filename = "/data/wiseardata/" + release + "/p3am_cdd/" + f1 + "/" + f2 + "/" + tile_id + "/" + tile_id + "-" + band + "-int-3.fits.gz"
			nhdu = 0
			if release == "unWISE":
				filename = wise_dir + "unwise-" + tile_id_unWISE + "-" + band + "-img-m.fits"
				nhdu = 1
			ax = fig.add_subplot(gs[0,s])
			for i in ax.spines.itervalues():
				i.set_linewidth(3.0)
			if not os.path.exists(filename) and os.path.exists(filename + ".fz"):
				filename = filename + ".fz"
			if not os.path.exists(filename):
				print "Downloading....", filename
				subprocess.call(["bash", "/home/sr525/bash_scripts/wise_tile_get.bash", filename])

			try:
				with fits.open(filename) as h:
					if verbose:
						print filename
					hdulist, im_status = make_cutout(h, RA, DEC, width, nhdu = nhdu)
					im = hdulist[1].data
					vmin, vmax = cutout_scale(im)
					ax.axes.get_xaxis().set_visible(False)
					ax.axes.get_yaxis().set_visible(False)
					ax.imshow(im, vmax = vmax, vmin = vmin, interpolation = "none")
					datas.append((im, vmin, vmax))
					ax.set_title(band)

			except IOError as e:
				#print e
				null_image = np.zeros(shape = (10,10), dtype = "int8")
				ax.imshow(null_image)
				ax.annotate("No data", xy = (0.01,0.9), size = "small")
				datas.append((np.zeros(29, 29), 0, 1))

			s += 1

		if not small: plt.suptitle(id)

	else:
		fig = plt.figure()
		null_image = np.zeros(shape = (10,10), dtype = "int8")
		plt.imshow(null_image)
		plt.annotate("No data", xy = (0.01,0.9), size = "small")

	if small:
		fig.subplots_adjust(wspace = 0.00, left = 0.05, right = 0.95, top = 1.0, bottom = 0.0, hspace = 0.0)
		fig = plt.gcf()
		fig.set_size_inches(6.7,2.0)
		fig.patch.set_facecolor("none")

	if return_val:
		plt.close(fig)
		return datas

	if save:
		plt.savefig(outdir + "/Cutouts_WISE_" + str(id1) + ".png", transparent = True)
		plt.close()

	else: return fig

def cutout_HST(outdir, RA, DEC, id, save = True, sigma = 0, width = 30.0, type = "png", return_data = False):

	"""
	Makes cutouts from the HST image of the COSMOS field
	Width is the full width and needs to be given in arcsecs
	Can return a fits file or a png
	Can also return just the array used to make the image and the scaling parameters
	sigma applies a gaussian filter
	"""

	import astropy.io.fits as fits
	from astropy import wcs
	import numpy as np
	import matplotlib.pyplot as plt
	from scipy.ndimage import gaussian_filter

	datas = []
	with open("/data/desardata/COSMOS/acs_mosaic_2.0/tiles/corners.txt", "r") as f:
		for line in f:
			l = line.split(",")
			ra_min = float(l[0])
			ra_max = float(l[2])
			dec_min = float(l[1])
			dec_max = float(l[3])
			filename = l[4]
			if RA < ra_min and RA > ra_max and DEC > dec_min and DEC < dec_max:
				print filename
				with fits.open(filename[:-1]) as h:
					hdulist, im_status = make_cutout(h, RA, DEC, width, nhdu = 0)
					im = hdulist[1].data
					vmin, vmax = cutout_scale(im)
					if sigma > 0:
						im = gaussian_filter(im, sigma)
					if type == "fits":
						hdulist.writeto(outdir + "/" + id + "_HST_" + band + ".fits")
					fig = plt.figure()
					ax = fig.add_subplot(111)
					ax.axes.get_xaxis().set_visible(False)
					ax.axes.get_yaxis().set_visible(False)
					ax.imshow(im, vmax = vmax, vmin = vmin)
					ax.set_title(id)
					if return_data:
						datas.append((im, vmin, vmax))
					fig = plt.gcf()
					fig.set_size_inches(2.5,2.5)

					if save:
						plt.savefig(outdir + "/Cutouts_HST_" + str(id) + ".png", transparent = True)
					if return_data:
						return datas
					else:
						return fig

def DES_Y2_cutout(RA, DEC, id, outdir, width = 30.0, save = False, cmap = False, crosshairs = False, bands = ["g", "r", "i", "z", "Y"], vmin_scale = False, vmax_scale = False, verbose = True, fits_out = False):

	"""
	Makes cutouts from the DES Y2 single epoch images
	Can specify the cmap with cmap = whatever you want
	If vmin_scale and vmax_scale are not set it scales image between median + 
	5*sigma MAD and median -2 * sigma MAD else uses median - vmin_scale * 
	sigma MAD and median + 5 * sigma MAD
	crosshairs = True draws cross hairs and a N and E arrow
	"""

	import matplotlib.pyplot as plt
	import matplotlib.gridspec as gridspec
	import os
	import numpy as np
	import DES_utils
	import astropy.io.fits as fits
	import gc
	import wcsutil

	files, out_bands = DES_utils.DES_Y2_image_file(RA, DEC, bands = bands)

	if len(files) >= 20:
		files = files[:20]
		out_bands = out_bands[:20]

	fig = plt.figure()
	print len(files)
	if len(files) < 5:
		gs = gridspec.GridSpec(1,len(files))
	else:
		x = int(np.ceil(len(files)/5.0))
		gs = gridspec.GridSpec(x,5)

	n = 0
	m = 0
	c = 0
	for file in files:
		ax = fig.add_subplot(gs[m, n])
		ax.axes.get_xaxis().set_visible(False)
		ax.axes.get_yaxis().set_visible(False)
		for i in ax.spines.itervalues():
			i.set_linewidth(3.0)
		if os.path.exists(file):
			if verbose: print file
			try:
				with fits.open(file) as h:
					hdulist, im_status = make_cutout(h, RA, DEC, width, nhdu = 0, w_units = "arcsecs")
					try:
						hdr = hdulist[1].header
						image = hdulist[1].data
					except (ValueError, TypeError):
						image = np.zeros((60,60))
						
					#image = hdulist[1].data
					#hdr = hdulist[1].header
			except IOError:
				image = np.zeros((10,10))

			data = image.flatten()
			med = np.median(data)
			sigma_MAD = 1.4826 * MAD(data, med)
			if not vmax_scale:
				vmax = med + 5.0*sigma_MAD
			else:
				vmax = med + vmax_scale*sigma_MAD
			if not vmin_scale:
				vmin = med - 2.0*sigma_MAD
			else:
				vmin = med - vmin_scale*sigma_MAD
			if vmax == vmin:
				vmax = np.max(data)
				vmin = np.min(data)
			if fits_out and hdulist is not None:
				hdulist.writeto(outdir + "/" + str(id) + "_" + out_bands[c] + ".fits", clobber = True)

			if cmap:
				ax.imshow(image, vmin = vmin, vmax = vmax, picker = True, interpolation = "nearest", cmap = cmap)
			else:
				ax.imshow(image, vmin = vmin, vmax = vmax, picker = True, interpolation = "nearest")
			if crosshairs and hdulist is not None:
				try:
					from astropy import wcs
					w = wcs.WCS(hdr, naxis = 2)
				except (UnboundLocalError, wcs.wcs.InconsistentAxisTypesError) as e:
					try:
						w = wcsutil.WCS(hdr)
					except KeyError:
						hdr = h[nhdu-1].header
						w = wcsutil.WCS(hdr)

				ax = draw_crosshairs(ax, w, RA, DEC)
			ax.set_title(out_bands[c])

			gc.collect()
			gc.collect()
			gc.collect()

		if len(files) > 5 and (n+1) % 5 == 0:
			n = 0
			m += 1
		else:
			n += 1
		c += 1

	if save:
		plt.savefig(outdir + "DESY2N_cutout_" + str(id) + ".png")
	else:
		return fig

def HST_DES_i(outdir, RA, DEC, tile, run, id, width = 30.0, save = False, release = "SVA1"):

	"""
	Makes an i band cutout from DES and HST scaled to not saturate in the middle
	Written for Gourav's webpage things
	"""

	import matplotlib.pyplot as plt

	HST_im = cutout_HST(outdir, RA, DEC, id, save = False, sigma = 0, width = width, return_val = "data")
	DES_im = cutout_image(outdir, RA, DEC, tile, run, id, save = False, width = width, field = release, bands = ["i"], return_val = "data")
	DES_im1 = cutout_image(outdir, RA, DEC, tile, run, id, save = False, width = width, field = release, bands = ["i"], return_val = "data", vmin_num = -2.0, vmax_num = 50.0)

	plt.close()
	plt.close()
	plt.close()
	fig = plt.figure()
	ax1 = fig.add_subplot(131)
	ax2 = fig.add_subplot(132)
	ax3 = fig.add_subplot(133)
	ax1.imshow(HST_im[0][0], vmin = HST_im[0][1], vmax = HST_im[0][2])
	ax2.imshow(DES_im[0][0], vmin = DES_im[0][1], vmax = DES_im[0][2])
	ax3.imshow(DES_im1[0][0], vmin = DES_im1[0][1], vmax = DES_im1[0][2])

	for axis in [ax1, ax2, ax3]:
		axis.axes.get_xaxis().set_visible(False)
		axis.axes.get_yaxis().set_visible(False)
		for i in axis.spines.itervalues():
			i.set_linewidth(3.0)
	fig.subplots_adjust(wspace = 0.00, left = 0.02, right = 0.98, top = 0.98, bottom = 0.02, hspace = 0.0)

	if save:
		plt.savefig(outdir + "HST_DES_i_" + id + ".png")
		plt.close()

	else: plt.show()

def composite(outdir, RA, DEC, tile, run, id, width = 30.0, save = False, field = "SVA1"):

	"""
	Creates two compostites from the DES images, gri and izY
	"""

	import os
	import numpy as np
	import matplotlib.pyplot as plt
	import img_scale
	
	ims = []
	sfs = [1.0, 1.0, 1.0, 1.0, 0.9]
	for (i,band) in enumerate(["g", "r", "i", "z", "Y"]):

		filename = "/data/desardata/" + field + "/" + tile + "/" + tile + "_" + band + ".fits.fz"
		if not os.path.exists(filename):
			subprocess.call(["bash", "/home/sr525/bash_scripts/wget_ims.bash", tile, str(run), field])
		hdulist, im_status = make_cutout(filename, RA, DEC, width, nhdu = 1)
		im = hdulist[1].data
		ims.append(im*sfs[i])

	smax, smin = img_scale.MAD(ims[0:3], 30, -10)
	img_blue = np.zeros((ims[2].shape[0], ims[2].shape[1], 3), dtype=float)
	j = 0
	for i in [2,1,0]:
		img_blue[:,:,j] = img_scale.linear(ims[i], scale_min = smin, scale_max = smax)
		gc.collect()
		gc.collect()
		j += 1	

	smax, smin = img_scale.MAD(ims[2:], 5, -10)
	img_red = np.zeros((ims[2].shape[0], ims[2].shape[1], 3), dtype=float)
	j = 0
	for i in [4,3,2]:
		img_red[:,:,j] = img_scale.linear(ims[i], scale_min = smin, scale_max = smax)
		j += 1
	
	figs = []
	j = 0
	for im in [img_blue, img_red]:
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.axes.get_xaxis().set_visible(False)
		ax.axes.get_yaxis().set_visible(False)
		for i in ax.spines.itervalues():
			i.set_linewidth(3.0)
		ax.imshow(im, aspect = "equal")
		figs.append(fig)
		if not save:
			plt.show()

		if save:
			if j == 0:
				plt.title("gri")
				fig.savefig(outdir + "Composite_gri_" + id + ".png")
			elif j == 1:
				plt.title("izY")
				fig.savefig(outdir + "Composite_izY_" + id + ".png")

def cutout_Spitzer(outdir, RA, DEC, id, save = False, width = 30.0):

	"""
	Makes a cutout from the COSMOS spitzer image
	Either saves the figure or returns it
	"""

	import matplotlib.pyplot as plt

	dir = "/data/desardata/COSMOS/Spitzer/irac/"
	bands = ["1", "2", "3", "4"]
	fig = plt.figure()
	gs = gridspec.GridSpec(1,4)
	s = 0
	for band in bands:
		filename = dir + "irac_ch" + band + "_go2_sci_10.fits"
		ax = fig.add_subpl8ot(gs[0,s])
		hdulist, im_status = make_cutout(filename, RA, DEC, width)
		im = hdulist[1].data
		vmin, vmax = cutout_scale(im)
		ax.axes.get_xaxis().set_visible(False)
		ax.axes.get_yaxis().set_visible(False)
		ax.imshow(im, vmin = vmin, vmax = vmax)
		ax.set_title("ch" + band)
		s += 1
	plt.suptitle(id)

	if save:
		plt.savefig(outdir + "/Cutouts_Spitzer_" + id + ".png", transparent = True)
	else:
		return fig

def cutout_UVISTA(outdir, RA, DEC, id, save = False, width = 30.0, type = "png"):

	"""
	Makes cutouts from the UVISTA image of the COSMOS field
	Width is in arcsecs
	If save = False then returns the figure
	"""

	import astropy.fits.io as fits
	import matplotlib.pyplot as plt

	dir = "/data/desardata/COSMOS/VISTA/"
	bands = ["H_26_04_11", "J_03_04_11", "Ks_15_12_10", "NB118_20_08_11", "Y_27_04_11"]
	fig = plt.figure()
	gs = gridspec.GridSpec(1,5)
	s = 0
	for band in bands:
		if band == "NB118_20_08_11":
			l = "a"
		else:
			l = ""
		filename = dir + "UVISTA_" + band + "_skysub_015_v1" + l + ".fits"
		ax = fig.add_subplot(gs[0,s])
		hdulist, im_status = make_cutout(filename, RA, DEC, width)
		im = hdulist[1].data
		vmin, vmax = cutout_scale(im)
		ax.axes.get_xaxis().set_visible(False)
		ax.axes.get_yaxis().set_visible(False)
		ax.imshow(im, vmin = vmin, vmax = vmax)
		ax.set_title(band)
		s += 1
		if type == "fits":
			hdulist.writeto(outdir + "/" + str(id) + "_" + band + ".fits", clobber = True)

	plt.suptitle(id)

	if save:
		plt.savefig(outdir + "/Cutouts_UVISTA_" + id + ".png", transparent = True)
	else:
		return fig

def cutout_GALEX(outdir, RA, DEC, id, save = False, width = 30.0):

	"""
	Makes cutouts from the GALEX image of the COSMOS fields
	If save = False then returns the figure
	Width is in arcsecs
	"""

	import matplotlib.pyplot as plt
	import matplotlib.gridspec as gridspec
	import astropy.io.fits as fits

	dir = "/data/desardata/COSMOS/GALEX/"
	bands = ["nd", "fd"]
	fig = plt.figure()
	gs = gridspec.GridSpec(1,8)
	s = 0
	for band in bands:
		with open(dir + band + "_corners.txt", "r") as f_corners:
			for line in f_corners:
				l = line.split(",")
				ra_min = float(l[0])
				ra_max = float(l[2])
				dec_min = float(l[1])
				dec_max = float(l[3])
				filename = l[4]
				if RA < ra_min and RA > ra_max and DEC > dec_min and DEC < dec_max and band in filename:
					ax = fig.add_subplot(gs[0,s])
					hdulist, im_status = make_cutout(filename, RA, DEC, width)
					im = hdulist[1].data
					vmin, vmax = cutout_scale(im)
					ax.axes.get_xaxis().set_visible(False)
					ax.axes.get_yaxis().set_visible(False)
					ax.imshow(im, vmin = vmin, vmax = vmax)
					ax.set_title(band)
					s += 1
	plt.suptitle(id)

	if save:
		plt.savefig(outdir + "/Cutouts_GALEX_" + id + ".png", transparent = True)
	else:
		return fig

def cutout_VLA(outdir, RA, DEC, id, save = False, width = 30.0):

	"""
	Makes cutouts from the VLA image of the COSMOS field
	width is in arcsecs
	If save = false returns the figure
	"""

	import matplotlib.pyplot as plt
	import matplotlib.gridspec as gridspec
	import astropy.io.fits as fits

	dir = "/data/desardata/COSMOS/VLA/"
	bands = ["lg_sin", "dp_sin", "pilot_sin"]
	fig = plt.figure()
	gs = gridspec.GridSpec(1,3)
	s = 0
	for band in bands:
		filename = dir + "vla_20cm_" + band + "_10.fits"
		ax = fig.add_subplot(gs[0,s])
		hdulist, im_status = make_cutout(filename, RA, DEC, width)
		im = hdulist[1].data
		vmin, vmax = cutout_scale(im)
		ax.axes.get_xaxis().set_visible(False)
		ax.axes.get_yaxis().set_visible(False)
		ax.imshow(im, vmin = vmin, vmax = vmax)
		ax.set_title(band)
		s += 1
	plt.suptitle(id)

	if save:
		fig = plt.gcf()
		fig.set_size_inches(7,3.5)
		plt.savefig(outdir + "/Cutouts_GALEX_" + id + ".png", transparent = True)
	else:
		return fig

def cutout_XMM(outdir, RA, DEC, id, save = False, width = 30.0, interpolation = "nearest", sens = False, exp = False, info = True):

	"""
	Makes cutouts from the XMM image of the COSMOS field
	I can't remember what all the options do as I've never used this code.
	"""

	import matplotlib.pyplot as plt
	import matplotlib.gridspec as gridspec
	import astropy.io.fits as fits

	dir = "/data/desardata/COSMOS/XMM/"
	bands = ["C0520", "C2080", "C4510"]
	fig = plt.figure()
	m = 3
	l = 1
	if sens: l += 1
	if exp: l += 1
	gs = gridspec.GridSpec(l,m)
	s = 0
	for band in bands:
		filename = dir + "xmm_" + band + ".img.fits"
		ax = fig.add_subplot(gs[0,s])
		hdulist, im_status = make_cutout(filename, RA, DEC, width)
		im = hdulist[1].data
		vmin, vmax = cutout_scale(im)
		ax.axes.get_xaxis().set_visible(False)
		ax.axes.get_yaxis().set_visible(False)
		ax.imshow(im, vmin = vmin, vmax = vmax, interpolation = interpolation)
		if info:
			med = np.median(im.flatten())
			sigma_mad = stat.MAD(im.flatten(), med)*1.4826
			ax.set_title(band + "\n" + med + "\n" + sigma_mad)
		else:
			ax.set_title(band)
		s += 1

	s = 0
	if sens:
		for band in bands:
			filename = dir + "xmm_C_" + band[1:] + "-10.sens.fits"
			ax = fig.add_subplot(gs[1,s])
			hdulist, im_status = make_cutout(filename, RA, DEC, width)
			im = hdulist[1].data
			vmin, vmax = cutout_scale(im)
			ax.axes.get_xaxis().set_visible(False)
			ax.axes.get_yaxis().set_visible(False)
			ax.imshow(im, vmin = vmin, vmax = vmax, interpolation = interpolation)
			ax.set_title("sens " + band)
			s += 1

	s = 0
	if exp:
		for band in bands:
			filename = dir + "xmm_" + band + "-vigyes-w.exp.fits"
			ax = fig.add_subplot(gs[2,s])
			hdulist, im_status = make_cutout(filename, RA, DEC, width)
			im = hdulist[1].data
			vmin, vmax = cutout_scale(im)
			ax.axes.get_xaxis().set_visible(False)
			ax.axes.get_yaxis().set_visible(False)
			ax.imshow(im, vmin = vmin, vmax = vmax, interpolation = interpolation)
			ax.set_title("exp " + band)
			s += 1

	plt.suptitle(id)

	if save:
		plt.savefig(outdir + "/Cutouts_XMM_" + id + ".png", transparent = True)
	else:
		return fig

def cutout_Chandra(outdir, RA, DEC, id, save = False, width = 30.0, interpolation = "nearest"):

	"""
	Makes cutouts from the Chandra image of the COSMOS field
	Width is in arcsecs
	If save = false then it returns the figure
	"""

	import matplotlib.pyplot as plt
	import matplotlib.gridspec as gridspec
	import astropy.io.fits as fits

	dir = "/data/desardata/COSMOS/Chandra/"
	bands = ["CC0520", "CC0570", "CC2070"]
	fig = plt.figure()
	gs = gridspec.GridSpec(1,3)
	s = 0
	for band in bands:
		filename = dir + band + "_img.fits"
		ax = fig.add_subplot(gs[0,s])
		hdulist, im_status = make_cutout(filename, RA, DEC, width)
		im = hdulist[1].data
		ax.axes.get_xaxis().set_visible(False)
		ax.axes.get_yaxis().set_visible(False)
		ax.imshow(im, interpolation = interpolation)
		ax.set_title(band)
		s += 1
	plt.suptitle(id)

	if save:
		plt.savefig(outdir + "/Cutouts_Chandra_" + id + ".png", transparent = True)
	else:
		return fig

def cutout_image(outdir, RA, DEC, tile, run, id, save = False, sigma = 0, cat_info = False, width = 30.0, type = "png", get_image = True, field = "SVA1", circle = False, circle_coords = [], small = False, bands = ["g", "r", "i", "z", "Y", "det"], return_val = False, vmax_num = 5, vmin_num =2, cat_band = False, crosshairs = False, cmap = False, verbose = True, depth = False):

	"""
	A function that has far more options than is sensible
	Makes multiband cutouts from the DES images
	outdir - where things are saved in save = True
	tile - the DES tile, needed to find the stored images at the IoA
	run - needed to download the image if we don't already have it, if we definitley
	have the images this can be set to ""
	id - used for the title
	sigma - if not 0 smooths the image
	cat_info - circles everything in the catalogue in the image
	width - in arcsecs
	type - png or fits
	get_image - governs whether or not files are downloaded. If false and the file 
	does not exist here the code returns blue boxes.
	circle - draws ellipses around the central object scaled by the kron radius and 
	the A, B and theta for the object
	circle coords - allows the user to specify where to draw some circles
	small - makes the resulting png smaller in size, eg. no title
	return_val - if = "data" code reutrns the image area and the scale factors
	vmin_num, vmax_num - scale factors
	Yband_cat - Uses the Y band based catalogues for information
	crosshairs - draws crosshairs
	cmap - allows the user to specify the colour map used
	"""

	import matplotlib.pyplot as plt
	import numpy as np
	import matplotlib.gridspec as gridspec
	import os
	import astropy.io.fits as fits
	import astropy.io.fits.compression
	import match_lists
	import subprocess
	from scipy.ndimage import gaussian_filter
	import matplotlib.cm as cm
	from astropy.table import Table, vstack, hstack
	from astropy import wcs
	import gc
	from matplotlib.patches import Ellipse, PathPatch
	import img_scale

	if not verbose:
		import warnings
		from astropy.utils.exceptions import AstropyWarning
		warnings.filterwarnings('ignore', category= AstropyWarning, append=True)
		warnings.filterwarnings('ignore', category= RuntimeWarning, append=True)

	fig = plt.figure()
	datas = []
	if len(bands) > 3:
		ls = len(bands) + 2
	else:
		ls = len(bands)
	gs = gridspec.GridSpec(1,ls)
	path = "/data/desardata/" + field + "/"
	if depth:
		path = "/data/desardata/" + field + "/" + depth + "/"
	path_list = [path + tile + "/"]
	a = 0
	im_scale = [0.40, 0.60, 0.85, 0.65, 0.15, 1.0]  
	s_bands = [1,1,1,1,1]
	im_list = []
	scale_list = []

	for path_use in path_list:
		outfile = None
		s = 0

		for band in bands:
			path_cat = "/data/desardata/" + field + "/" + tile + "/" + tile + "_" + band + "_cat.fits"
			if cat_info or circle:
				circle_details = []
				if cat_band:
					path_cat = "/data/desardata/" + field + "/" + tile + "/" + cat_band + "Detection/" + tile + "_" + band + "_cat.fits"

				if os.path.exists(path_use + tile + "_" + band + "_cat.fits") == False and band <> "det":
					print "Getting cat"
					subprocess.call(["bash", "/home/sr525/bash_scripts/wget_cat.bash", tile, str(run), field])

				if band == "det" or os.path.exists(path_use + tile + "_" + band + "_cat.fits") == False:
					a = 0.0
					b = 0.0
					theta = 0.0

				else:
					with fits.open(path_cat) as h:
						data = h[1].data
						dists, inds = match_lists.match_lists([RA], [DEC], data["ALPHAPEAK_J2000"], data["DELTAPEAK_J2000"], width*np.sqrt(2)*2, numNei = 50)
						ids = np.where((inds <> len(data)))
						circle_coords = zip(data[inds[ids]]["ALPHAPEAK_J2000"],data[inds[ids]]["DELTAPEAK_J2000"])
						if circle or (cat_info and circle_coords == []):
							circle_coords = [(RA, DEC)]
						for (ra, dec) in circle_coords:
							dist, ind = match_lists.match_lists([ra], [dec], data["ALPHAPEAK_J2000"], data["DELTAPEAK_J2000"], 0.0006)
							if ind < len(data):
								if cat_band:
									k = 10.0
								else:
									k = data["KRON_RADIUS"][ind]
								a = data["A_IMAGE"][ind]*k
								b = data["B_IMAGE"][ind]*k
								theta = data["THETA_IMAGE"][ind]
							else:
								a = 0
								b = 0
								theta = 0
							circle_details.append((ra, dec, a, b, theta))

			else:
				a = 0
				b = 0
				theta = 0

			p = bands.index(band)
			ax = fig.add_subplot(gs[0, s])
			ax.axes.get_xaxis().set_visible(False)
			ax.axes.get_yaxis().set_visible(False)
			for i in ax.spines.itervalues():
				i.set_linewidth(3.0)
			im_name = tile + "_"
			filename = path_use + im_name + band + ".fits.fz"
			if get_image and not os.path.exists(filename):
				subprocess.call(["bash", "/home/sr525/bash_scripts/wget_ims.bash", tile, str(run), field])
			if verbose:
				print filename
			if os.path.exists(filename):
				with fits.open(filename) as h:
					hdulist, im_status = make_cutout(h, RA, DEC, width, nhdu = 1, w_units = "arcsecs", verbose = verbose)
					image = hdulist[1].data
					hdr = hdulist[1].header
					w = wcs.WCS(hdr, naxis = 2)
					if circle or cat_info:
						ells = []
						for (ra, dec, a, b, theta) in circle_details:
							xy_c = w.wcs_world2pix(np.array([[ra, dec]], np.float_), 1)
							xy_c = (xy_c[0][0], xy_c[0][1])
							ells.append(Ellipse(xy = xy_c, width = a, height = b, angle = theta, fill = False))

					if type == "fits":
						hdulist.writeto(outdir + "/" + str(id) + "_" + band + ".fits", clobber = True)
					image = image * (1.0 / im_scale[p])
					vmin, vmax = cutout_scale(image, num_min = vmin_num, num_max = vmax_num)
					if sigma > 0:
						image = gaussian_filter(image, sigma)
					scale_list.append((vmax, vmin))
					if return_val == "data":
						datas.append((image, vmin, vmax))
					if cmap:
						ax.imshow(image, vmin = vmin, vmax = vmax, picker = True, interpolation = "nearest", cmap = cmap)
					else:
						ax.imshow(image, vmin = vmin, vmax = vmax, picker = True, interpolation = "nearest")
					if crosshairs:
						ax = draw_crosshairs(ax, w, RA, DEC)
					ax.set_title(band)
					if return_val == "ax":
						datas.append(ax)
					im_list.append(image)
					if circle or cat_info:
						for ell in ells:
							ax.add_patch(ell)

			else:
				null_image = np.zeros(shape = (10,10), dtype = "int8")
				ax.imshow(null_image)
				ax.annotate("No " + band + " band data", xy = (0.01,0.9), size = "small")
				if s == 0:
					image = np.zeros(shape = (10,10), dtype = "int8")
				else:
					image = np.zeros((im_list[s-1].shape[0], im_list[s-1].shape[1]), dtype = float)
				im_list.append(image)
				ax.set_title(band)

			s += 1

	gc.collect()
	gc.collect()
	gc.collect()

	if len(bands) > 3:
		img_blue = np.zeros((im_list[0].shape[0], im_list[0].shape[1], 3), dtype=float)
		j = 0
		smin, smax = cutout_scale(np.array(im_list[0:3]).flatten(), num_min = 1, num_max = 30)
		for i in [2,1,0]:
			img_blue[:,:,j] = img_scale.linear(im_list[i], scale_min = smin, scale_max = smax)
			gc.collect()
			gc.collect()
			j += 1	

		smin, smax = cutout_scale(np.array(im_list[2:5]).flatten(), num_min = 1, num_max = 30)
		img_red = np.zeros((im_list[2].shape[0], im_list[2].shape[1], 3), dtype=float)
		j = 0
		for i in [4,3,2]:
			img_red[:,:,j] = img_scale.linear(im_list[i], scale_min = smin, scale_max = smax)
			j += 1
	
		for im in [img_blue, img_red]:
			ax = fig.add_subplot(gs[0, s])
			ax.axes.get_xaxis().set_visible(False)
			ax.axes.get_yaxis().set_visible(False)
			for i in ax.spines.itervalues():
					i.set_linewidth(3.0)
			ax.imshow(im, aspect = "equal")
			if s == 6:
				ax.set_title("gri")
			if s == 7:
				ax.set_title("izY")
			s += 1

	if not small:
		fig.suptitle(id)
		fig.subplots_adjust(wspace = 0.001, left = 0.01, right = 0.99, top = 0.85, bottom = 0.01)
		fig = plt.gcf()
		fig.set_size_inches(7,3.5)
		fig.patch.set_facecolor("none")

	if small:
		fig.subplots_adjust(wspace = 0.00, left = 0.05, right = 0.95, top = 1.0, bottom = 0.0, hspace = 0.0)
		fig = plt.gcf()
		fig.set_size_inches(8,1.5)
		fig.patch.set_facecolor("none")
	
	if save:
		if sigma > 0:
			fig.savefig(outdir + "/Cutouts_Smoothed_" + str(id) + ".png", transparent = True)
		elif circle:
			fig.savefig(outdir + "/Cutouts_Circle_" + str(id) + ".png", transparent = True)
		elif cat_info or (cat_info and circle):
			fig.savefig(outdir + "/Cutouts_Cat_Info_" + str(id) + ".png", transparent = True)
		else:
			fig.savefig(outdir + "/Cutouts_" + str(id) + ".png", transparent = True)
		plt.close(fig)
	if return_val:
		plt.close(fig)
		return datas
	else:
		return fig

def des_vhs_cutout(outdir, RA, DEC, tile, run, id, save = False, sigma = 0, cat_info = False, width = 30.0, get_image = False, field = "SVA1", circle = False, circle_coords = [], small = False, Yband_cat = False, crosshairs = False, cmap = False, return_val = False):

	"""
	Code that is very similar to cutout_image but also puts VHS cutouts on the end
	Probably full of unnecessary repetition
	"""

	import matplotlib.pyplot as plt
	import numpy as np
	import matplotlib.gridspec as gridspec
	import os
	import astropy.io.fits as fits
	import astropy.io.fits.compression
	import match_lists
	import subprocess
	from scipy.ndimage import gaussian_filter
	import matplotlib.cm as cm
	from astropy.table import Table, vstack, hstack
	from astropy import wcs
	import gc
	from matplotlib.patches import Ellipse, PathPatch
	import cutout_vista
	import xmlrpclib

	fig = plt.figure()
	datas = []
	bands = ["g", "r", "i", "z", "Y"]
	gs = gridspec.GridSpec(1,len(bands)+2)
	path = "/data/desardata/" + field + "/" + tile + "/"
	s = 0
	j = 0
	k = 0

	for band in bands:
		path_cat = "/data/desardata/" + field + "/" + tile + "/" + tile + "_" + band + "_cat.fits"
		if cat_info or circle:
			circle_details = []
			if Yband_cat:
				path_cat = "/data/desardata/" + field + "/" + tile + "/YDetection/" + tile + "_" + band + "_cat.fits"

			if os.path.exists(path_use + tile + "_" + band + "_cat.fits") == False and band <> "det":
				print "Getting cat"
				subprocess.call(["bash", "/home/sr525/bash_scripts/wget_cat.bash", tile, str(run), field])

			if band == "det" or os.path.exists(path_use + tile + "_" + band + "_cat.fits") == False:
				a = 0.0
				b = 0.0
				theta = 0.0

			else:
				with fits.open(path_cat) as h:
					data = h[1].data
					dists, inds = match_lists.match_lists([RA], [DEC], data["ALPHAWIN_J2000"], data["DELTAWIN_J2000"], width*np.sqrt(2)*2, numNei = 50)
					ids = np.where((inds <> len(data)))
					circle_coords = zip(data[inds[ids]]["ALPHAWIN_J2000"],data[inds[ids]]["DELTAWIN_J2000"])
					if circle or (cat_info and circle_coords == []):
						circle_coords = [(RA, DEC)]
					for (ra, dec) in circle_coords:
						dist, ind = match_lists.match_lists([ra], [dec], data["ALPHAWIN_J2000"], data["DELTAWIN_J2000"], 0.0006)
						if ind < len(data):
							if Yband_cat:
								k = 10.0
							else:
								k = data["KRON_RADIUS"][ind]
							a = data["A_IMAGE"][ind]*k
							b = data["B_IMAGE"][ind]*k
							theta = data["THETA_IMAGE"][ind]
						else:
							a = 0
							b = 0
							theta = 0
						circle_details.append((ra, dec, a, b, theta))

		else:
			a = 0
			b = 0
			theta = 0

		p = bands.index(band)
		ax = fig.add_subplot(gs[0, s])
		ax.axes.get_xaxis().set_visible(False)
		ax.axes.get_yaxis().set_visible(False)
		for i in ax.spines.itervalues():
			i.set_linewidth(3.0)
		im_name = tile + "_"
		filename = path + im_name + band + ".fits.fz"
		if get_image and not os.path.exists(filename):
			subprocess.call(["bash", "/home/sr525/bash_scripts/wget_ims.bash", tile, str(run), field])
		print filename
		if os.path.exists(filename):
			with fits.open(filename) as h:
				hdulist, im_status = make_cutout(h, RA, DEC, width, nhdu = 1, w_units = "arcsecs")
				image = hdulist[1].data
				hdr = hdulist[1].header
				w = wcs.WCS(hdr, naxis = 2)
				if circle or cat_info:
					ells = []
					for (ra, dec, a, b, theta) in circle_details:
						xy_c = w.wcs_world2pix(np.array([[ra, dec]], np.float_), 1)
						xy_c = (xy_c[0][0], xy_c[0][1])
						ells.append(Ellipse(xy = xy_c, width = a, height = b, angle = theta, fill = False))

				if type == "fits":
					hdulist.writeto(outdir + "/" + str(id) + "_" + band + ".fits", clobber = True)
				data = image.flatten()
				med = np.median(data)
				sigma_MAD = 1.4826 * MAD(data, med)
				vmax = med + 5.0*sigma_MAD
				vmin = med - 2.0*sigma_MAD
				if sigma > 0:
					image = gaussian_filter(image, sigma)
				if cmap:
					ax.imshow(image, vmin = vmin, vmax = vmax, picker = True, interpolation = "nearest", cmap = cmap)
				else:
					ax.imshow(image, vmin = vmin, vmax = vmax, picker = True, interpolation = "nearest")
				datas.append((image, vmin, vmax))
				if crosshairs:
					ax = draw_crosshairs(ax, w, RA, DEC)
				#ax.imshow(image, interpolation = "nearest")
				ax.set_title(band, fontsize = 25)
				if circle or cat_info:
					for ell in ells:
						ax.add_patch(ell)

		else:
			null_image = np.zeros(shape = (10,10), dtype = "int8")
			ax.imshow(null_image)
			ax.annotate("No " + band + " band data", xy = (0.01,0.9), size = "small")
			if s == 0:
				image = np.zeros(shape = (10,10), dtype = "int8")
			else:
				image = np.zeros((image.shape[0], image.shape[1]), dtype = float)
			datas.append((image, 0, 1))
			ax.set_title(band)

		s += 1

		gc.collect()
		gc.collect()
		gc.collect()

	files = []
	try:
		files, apm_files, chipnos = cutout_vista.querydb(RA, DEC, size = 30, db="vista", onlytiles=True, ctype='fit')
		l = len(files)
	except xmlrpclib.Fault as e:
		if e.faultCode == 1:
			#faultCode 1 is IOError
			print "Data not found"
			null_image = np.zeros(shape = (10,10), dtype = "int8")
			ax = fig.add_subplot(111)
			ax.imshow(null_image)
			ax.annotate("No data", xy = (0.01,0.9), size = "small")
			ax.set_title("Data not found")
		else:
			raise e

	for file in files:
		if ("J" in file and j < 1) or ("K" in file and k < 1):
			ax = fig.add_subplot(gs[0,s])
			plt.tick_params(axis='x',which='both', bottom='off', top='off', labelbottom='off')
			plt.tick_params(axis='y',which='both', left='off', right='off', labelleft='off')
			for i in ax.spines.itervalues():
				i.set_linewidth(3.0)

			im_file = fits.open(file)
			hdr = im_file[0].header
			w = wcs.WCS(hdr, naxis = 2)
			im = im_file[0].data
			data = im.flatten()
			med = np.median(data)
			sigma_MAD = 1.4826 * MAD(data, med)
			vmax = med + 5.0*sigma_MAD
			vmin = med - 2.0*sigma_MAD
			if cmap:
				ax.imshow(im, vmax = vmax, vmin = vmin, cmap = cmap, interpolation = "nearest")
			else:
				ax.imshow(im, vmax = vmax, vmin = vmin, interpolation = "nearest")
			datas.append((im, vmin, vmax))
			if crosshairs:
				ax = draw_crosshairs(ax, w, RA, DEC)

			med = str(med)
			sigma_MAD = str(sigma_MAD)
			if "J" in file:
				ax.set_title("J", fontsize = 25)
				j += 1
			elif "K" in file:
				ax.set_title("K", fontsize = 25)
				k += 1

			s += 1

	fig.subplots_adjust(wspace = 0.001, left = 0.01, right = 0.99, top = 0.85, bottom = 0.01)
	fig = plt.gcf()
	fig.set_size_inches(7,3.5)
	fig.patch.set_facecolor("none")

	if return_val:
		plt.close(fig)
		return datas

	if save:
		if sigma > 0:
			fig.savefig(outdir + "/Cutouts_Smoothed_" + str(id) + ".png", transparent = True)
		else:
			fig.savefig(outdir + "/Cutouts_" + str(id) + ".png", transparent = True)
		plt.close(fig)
	else:
		return fig

def cutout_seg(outdir, RA, DEC, tile, run, id, save = True, sigma = 0, small = False, release = "Y1A1", width_arcsecs = 30.0, out_type = "fits", bands = ["g", "r", "i", "z", "Y"]):

	"""
	Makes cutouts of the DES segmentation maps
	"""

	import matplotlib.pyplot as plt
	import matplotlib.gridspec as gridspec
	import os
	import astropy.io.fits as fits
	from scipy.ndimage import gaussian_filter
	import subprocess
	import numpy as np
	import matplotlib.cm as cm

	fig = plt.figure()
	gs = gridspec.GridSpec(1,5)
	path = "/data/desardata/" + release + "/" + tile + "/"
	path_list = [path]
	s_bands = [1,1,1,1,1]
	im_list = []
	scale_list = []
	filename_use = False
	datas = []

	for path_use in path_list:
		outfile = None
		s = 0

		for band in bands:
			p = bands.index(band)
			ax = fig.add_subplot(gs[0, s])
			ax.axes.get_xaxis().set_visible(False)
			ax.axes.get_yaxis().set_visible(False)
			for i in ax.spines.itervalues():
				i.set_linewidth(3.0)
			im_name = tile + "_"
			filenames = [path_use + im_name + band + "_seg.fits.fz", path_use + im_name + band + "_seg.fits"]
			if not os.path.exists(filenames[0]) and not os.path.exists(filenames[1]):
				subprocess.call(["bash", "/home/sr525/bash_scripts/wget_segmap.bash", tile, run, release])
			for filename in filenames:
				if os.path.exists(filename):
					filename_use = filename
			if filename_use is not False:
				with fits.open(filename_use) as h:
					hdulist, im_status = make_cutout(h, RA, DEC, width_arcsecs, nhdu = 1, w_units = "arcsecs")
					image = hdulist[1].data
					hdr = hdulist[1].header
					if out_type == "fits":
						hdulist.writeto(outdir + "/" + str(id) + "_" + band + "_seg.fits", clobber = True)
					if out_type == "data":
						datas.append(hdulist)
				if sigma > 0:
					image = gaussian_filter(image, sigma)
				ax.imshow(image, picker = True, cmap = cm.binary)
				ax.set_title(band)

			else:
				null_image = np.zeros(shape = (10,10), dtype = "int8")
				ax.imshow(null_image)
				ax.annotate("No " + band + " band data", xy = (0.01,0.9), size = "small")
				image = np.zeros(shape = (10,10), dtype = "int8")
				ax.set_title(band)

			s += 1

	if not small:
		fig.suptitle(id)
		fig.subplots_adjust(wspace = 0.001, left = 0.01, right = 0.99, top = 0.85, bottom = 0.01)
		fig = plt.gcf()
		fig.set_size_inches(7,3.5)
		fig.patch.set_facecolor("none")

	if small:
		fig.subplots_adjust(wspace = 0.00, left = 0.05, right = 0.95, top = 1.0, bottom = 0.0, hspace = 0.0)
		fig = plt.gcf()
		fig.set_size_inches(8,1.5)
		fig.patch.set_facecolor("none")

	if save:
		if sigma > 0:
			fig.savefig(outdir + "/Cutouts_seg_smoothed_" + str(id) + ".png", transparent = True)
		else:
			fig.savefig(outdir + "/Cutouts_seg_" + str(id) + ".png", transparent = True)
		plt.close(fig)
		if out_type == "data":
			return datas
	else:
		if out_type == "data":
			return datas, fig
		else:
			return fig

def cutout_weight(outdir, RA, DEC, tile, id, save = True, sigma = 0, small = True, release = "SVA1", type = "png", bands = ["g", "r", "i", "z", "Y", "det"], width = 30.0):

	"""
	Makes cutouts from the DES weight maps
	They are 30" wide
	"""

	import matplotlib.pyplot as plt
	from scipy.ndimage import gaussian_filter
	import matplotlib.gridspec as gridspec
	import os
	import numpy as np
	import astropy.io.fits as fits

	fig = plt.figure()
	gs = gridspec.GridSpec(1,6)
	path = "/data/desardata/" + release + "/"
	path_list = [path + tile + "/"]
	a = 0
	s_bands = [1,1,1,1,1]
	im_list = []
	scale_list = []

	for path_use in path_list:
		outfile = None
		s = 0

		for band in bands:
			p = bands.index(band)
			ax = fig.add_subplot(gs[0, s])
			ax.axes.get_xaxis().set_visible(False)
			ax.axes.get_yaxis().set_visible(False)
			for i in ax.spines.itervalues():
				i.set_linewidth(3.0)
			im_name = tile + "_"
			filename = path_use + im_name + band + ".fits.fz"
			if os.path.exists(filename):
				with fits.open(filename) as hlist:
					hdulist, im_status = make_cutout(hlist, RA, DEC, width, nhdu = 2)
					image = hdulist[1].data
				if type == "fits":
					hdulist.writeto(outdir + "/Weight_" + str(id) + "_" + band + ".fits", clobber = True)
				im_list.append(image)
				vmin, vmax = cutout_scale(image)
				if sigma > 0:
					image = gaussian_filter(image, sigma)
				if band == "det":
					ax.imshow(image)
				else:
					ax.imshow(image, vmin = vmin, vmax = vmax, picker = True)
				scale_list.append((vmax, vmin))
				ax.set_title(band)

			else:
				null_image = np.zeros(shape = (10,10), dtype = "int8")
				ax.imshow(null_image)
				ax.annotate("No " + band + " band data", xy = (0.01,0.9), size = "small")
				if s == 0:
					image = np.zeros(shape = (10,10), dtype = "int8")
					im_list.append(image)
				else:
					image = np.zeros((im_list[s-1].shape[0], im_list[s-1].shape[1]), dtype = float)
					im_list.append(image)
				ax.set_title(band)

			s += 1

	if not small:
		fig.suptitle(id)
		fig.subplots_adjust(wspace = 0.001, left = 0.01, right = 0.99, top = 0.85, bottom = 0.01)
		fig = plt.gcf()
		fig.set_size_inches(7,3.5)
		fig.patch.set_facecolor("none")

	if small:
		fig.subplots_adjust(wspace = 0.00, left = 0.05, right = 0.95, top = 1.0, bottom = 0.0, hspace = 0.0)
		fig = plt.gcf()
		fig.set_size_inches(8,2.0)
		fig.patch.set_facecolor("none")

	if save:
		if sigma > 0:
			fig.savefig(outdir + "/Cutouts_weight_smoothed_" + str(id) + ".png", transparent = True)
		else:
			fig.savefig(outdir + "/Cutouts_weight_" + str(id) + ".png", transparent = True)
		plt.close(fig)
	else:
		return fig

def cutout_zoom(outdir, RA, DEC, tile, id, save = False, small = False, release = "SVA1", depth = False):

	"""
	Makes cutouts of the same object at a range of zooms in each band
	"""

	import matplotlib.pyplot as plt
	import os
	import matplotlib.gridspec as gridspec
	import astropy.io.fits as fits
	import astropy.io.fits.compression
	import gc
	import numpy as np

	fig = plt.figure()
	gs = gridspec.GridSpec(5,5)
	path = "/data/desardata/" + release + "/"
	if depth:
		path = "/data/desardata/" + release + "/" + depth + "/"
	path_list = [path + tile + "/"]

	for path_use in path_list:
		outfile = None
		r = 0
		bands = ["g", "r", "i", "z", "Y"]
		for band in bands:
			s1 = 0
			im_name = tile + "_"
			filename = path_use + im_name + band + ".fits.fz"

			if os.path.exists(filename):
				with fits.open(filename) as h:

					widths = [10, 20, 40, 80, 160]
					for width in widths:
						ax = fig.add_subplot(gs[r,s1])
						plt.tick_params(axis='x',which='both', bottom='off', top='off', labelbottom='off')
						plt.tick_params(axis='y',which='both', left='off', right='off', labelleft='off')
						for i in ax.spines.itervalues():
							i.set_linewidth(3.0)
						hdulist, im_status = make_cutout(h, RA, DEC, width, nhdu = 1, w_units = "arcsecs")
						image = hdulist[1].data
						vmin, vmax = cutout_scale(image)
						ax.imshow(image, vmin = vmin, vmax = vmax, interpolation = "nearest")
						if bands.index(band) == 0:
							ax.set_title(str(width), fontsize = "10")
						if widths.index(width) == 0:
							ax.set_ylabel(band + "   ", rotation = 0)
						s1 += 1

			else:
				for width in [10,20,40,80,160]:
					ax = fig.add_subplot(gs[r,s1])
					null_image = np.zeros(shape = (10,10), dtype = "int8")
					ax.imshow(null_image)
					ax.annotate("No " + band + " band data", xy = (0.01,0.9), size = "small")
					s1 += 1

			r += 1

		if not small:
			fig.suptitle(id)
			fig.subplots_adjust(wspace = 0.001, left = 0.10, right = 0.95, top = 0.85, bottom = 0.05, hspace = 0)
			fig = plt.gcf()
			fig.set_size_inches(6,6.9, forward = True)

		if small:
			fig.subplots_adjust(wspace = 0.000, left = 0.05, right = 0.94, top = 0.95, bottom = 0.02, hspace = 0)
			fig = plt.gcf()
			fig.set_size_inches(5.0,4.8, forward = True)

		gc.collect()
		gc.collect()
		gc.collect()

		if save:
			fig.savefig(outdir + "/Cutouts_zoom_" + str(id) + ".png")
			plt.close(fig)
		else:
			return fig

def cutout_VISTA(outdir, RA, DEC, tile, id, save = False, db = "vista", circle = True, radius = 0.0006*3600.0, size = 30.0, onlytiles = True, fits_out = False, cat_info = False, conf = False, return_val = False):
	
	"""
	Makes cutouts from one of the VISTA surveys
	db = vista - VHS
	db = wfcam - UKIDSS
	db = vst - VST
	onlytile = False - gives non stack images as well
	conf = True - adds the confidence map to the output fits file
	Based off Eduardo's code
	"""

	import matplotlib.pyplot as plt
	import cutout_vista
	import xmlrpclib
	import matplotlib.gridspec as gridspec
	import subprocess
	import astropy.io.fits as fits
	from astropy import wcs
	from matplotlib.patches import Ellipse, PathPatch
	import numpy as np
	import os
	import sqlutil
	import ConfigParser

	config_file = "/home/sr525/Python_Code/srpylib/srpylib.cfg"
	config = ConfigParser.RawConfigParser()
	config.read(config_file)
	uname = config.get("wsdb", "uname")
	pword = config.get("wsdb", "pword")
	host = config.get("wsdb", "host")

	fig = plt.figure()
	files = []
	datas = []
	try:
		files, apm_files, chipnos = cutout_vista.querydb(RA, DEC, size = size, db=db, onlytiles=onlytiles, ctype='fit')
		if conf:
			files_conf, apm_files_conf, chipnos_conf = cutout_vista.querydb(RA, DEC, size = size, db=db, onlytiles=onlytiles, ctype='fit', conf = True)
		#print apm_files
		l = len(files)
		if l > 20:
			files = files[:20]
			apm_files = apm_files[:20]
			chipnos = chipnos[:20]
		l = len(files)
		if l > 4.0:
			r = int(np.ceil(l/4.0))
			w = 4
		else:
			r = 1
			w = l
		if files == []:
			ax = fig.add_subplot(111)
			null_image = np.zeros(shape = (10,10), dtype = "int8")
			ax.imshow(null_image)
			ax.annotate("No data", xy = (0.01,0.9), size = "small")
			ax.set_title("No Match")
		else:
			gs = gridspec.GridSpec(r+1, w)
			ax2 = fig.add_subplot(gs[0,:])
			ax2.axes.get_xaxis().set_visible(False)
			ax2.axes.get_yaxis().set_visible(False)

	except xmlrpclib.Fault as e:
		if e.faultCode == 1:
			#faultCode 1 is IOError
			print "Data not found"
			null_image = np.zeros(shape = (10,10), dtype = "int8")
			ax = fig.add_subplot(111)
			ax.imshow(null_image)
			ax.annotate("No data", xy = (0.01,0.9), size = "small")
			ax.set_title("Data not found")
		else:
			raise e

	r = 1
	s1 = 0

	if not conf:
		files_conf = []
		for file in files:
			files_conf.append("")
	for (file, conf_file) in zip(files, files_conf):
		if fits_out and not conf:
			subprocess.call(["cp", file, outdir + file[4:]])

		im_file = fits.open(file)
		im = im_file[0].data

		if fits_out and conf:
			h_conf = fits.open(conf_file)
			phdu = fits.PrimaryHDU()
			imhdu = fits.ImageHDU(header = im_file[0].header, data = im)
			conf_hdu = fits.ImageHDU(header = h_conf[0].header, data = h_conf[0].data)
			hdulist = fits.HDUList([phdu, imhdu, conf_hdu])
			print outdir + conf_file[5:]
			hdulist.writeto(outdir + conf_file[5:], clobber = True)

		ax2.annotate(file, xy = (0.05, 0.9 - (0.6*(r-1)+0.15*s1)), fontsize = "x-small")
		ax = fig.add_subplot(gs[r,s1])
		plt.tick_params(axis='x',which='both', bottom='off', top='off', labelbottom='off')
		plt.tick_params(axis='y',which='both', left='off', right='off', labelleft='off')
		for i in ax.spines.itervalues():
			i.set_linewidth(3.0)

		vmin, vmax = cutout_scale(im)
		ax.imshow(im, vmax = vmax, vmin = vmin, interpolation = "none")

		if circle:
			a = (radius/30.0)*(len(im)/2.0)*2.0
			b = a
			theta = 0.0
			hdr = im_file[0].header
			w = wcs.WCS(hdr)
			xy = w.wcs_world2pix(np.array([[RA, DEC]], np.float_), 1)
			xy = (xy[0][0], xy[0][1])
			ell = Ellipse(xy = xy, width = a, height = b, angle = theta, fill = False)
			ax.add_patch(ell)

		if cat_info:
			a = (radius/30.0)*(len(im)/2.0)*2.0
			b = a
			theta = 0.0
			R = str(30.0/3600.0)
			query = "select ra, dec, japercormag3, japercormag3_err, " + \
					"hapercormag3, hapercormag3_err, kapercormag3, " + \
					"kapercormag3_err FROM vhs_201409.des " + \
					"WHERE q3c_radial_query(ra, dec, " + str(RA) + ", " + \
					str(DEC) + ", " + R + ")"

			RAs, DECs, J, Jerr, H, Herr, K, Kerr = sqlutil.get(query, db="wsdb", \
				host=host, user=uname, password=pword)
			hdr = im_file[0].header
			w = wcs.WCS(hdr)
			for (ra, dec) in zip(RAs, DECs):
				xy = w.wcs_world2pix(np.array([[ra, dec]], np.float_), 1)
				xy = (xy[0][0], xy[0][1])
				ell = Ellipse(xy = xy, width = a, height = b, angle = theta, fill = False)
				ax.add_patch(ell)

		#med = np.median(im)
		#sigma_MAD
		#med = str(med)
		#sigma_MAD = str(sigma_MAD)
		if "J" in file:
			band = "J"
		elif "K" in file:
			band = "K"
		elif "H" in file:
			band = "H"
		elif "u" in file:
			band = "u"
		elif "Y" in file:
			band = "Y"
		elif "r" in file:
			band = "r"
		elif "i" in file[9:-4]:
			band = "i"
		elif "z" in file:
			band = "z"
		elif "g" in file[9:]:
			band = "g"
		ax.set_title(band)
		os.remove(file)
		s1 += 1
		if s1 >= 4:
			s1 = 0
			r += 1

		if return_val:
			datas.append((im, vmin, vmax, band))

	fig.suptitle(id)
	fig.subplots_adjust(wspace = 0.001, left = 0.01, right = 0.99, top = 0.85, bottom = 0.01)
	fig = plt.gcf()
	#fig.set_size_inches(6,3.0)

	if return_val:
		plt.close(fig)
		return datas

	if save: 
		plt.savefig(outdir + "/Cutouts_" + db + "_" + str(id) + ".png")
		plt.close()
	else:
		return fig



