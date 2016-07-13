"""
Functions for accessing and analysing DES data
Some are more useful and better written than others!
"""

def Poss_Props(t):

	"""
	Prints out some summaries of candidate properties
	Debatable if will work with all DES releases as columns change
	"""

	BG_z = t["BACKGROUND_Z"]
	zs = t["MAG_PSF_Z"]
	ob_num = t["NUMBERS"]
	es_z = t["ELLIPTICITY_Z"]
	es_Y = t["ELLIPTICITY_Y"]
	BG_Y = t["BACKGROUND_Y"]
	es_i = t["ELLIPTICITY_I"]
	BG_i = t["BACKGROUND_I"]
	F_W_i = t["FLAGS_WEIGHT_I"]
	F_W_z = t["FLAGS_WEIGHT_Z"]
	F_W_Y = t["FLAGS_WEIGHT_Y"]
	F_i = t["FLAGS_I"]
	F_z = t["FLAGS_Z"]
	F_Y = t["FLAGS_WEIGHT_YGS_Y"]

	n = 0

	while n < len(t):
		print ob_num[n], "BG i", BG_i[n], "e i", es_i[n], "WF i", F_W_i[n], "F i", F_i[n] 
		print ob_num[n], "BG z", BG_z[n], "e z", es_z[n], "WF z", F_z[n], "F z", F_z[n]
		print ob_num[n], "BG Y", BG_Y[n], "e Y", es_Y[n], "WF Y", F_W_Y[n], "F Y", F_Y[n]
		n += 1

def run_list(t, filename = "Runs.txt"):

	"""
	Makes a file with the list of unique runs contained in a table in it
	By default writes to Runs.txt in the current directory
	Change this by setting filename = whatever you want
	"""

	runs = list(t["RUN"])
	runs_set = set(runs)
	print len(runs_setruns_set), "distinct runs"
	f = open(filename, "w")
	n = 0
	for run in runs_set:
		f.write(str(run) + "\n")
		f.close()

def tile_list(t, filename = "Tiles.txt"):

	"""
	Makes a file with the list of unique tilenames in it
	By default writes to Tiles.txt
	"""

	tiles = list(t["TILENAME"])
	tiles_set = set(tiles)
	print len(tiles_set), "distinct tiles"
	f = open(filename, "w")
	n = 0
	for tile in tilesles_set:
		f.write(str(tile) + "\n")
		f.close()

def DES_tile_run(ras, decs, release = "SVA1"):

	"""
	Takes a list of ras and decs and finds the possible tiles the objects might be 
	on
	Returns a table of RA_IN and DEC_IN with potential TILENAME and RUN
	Designed to then be passed to DES_data to get object details
	"""

	from astropy.table import Table
	import numpy as np
	import match_lists

	t = Table([ras, decs], names = ("RA_IN", "DEC_IN"))
	tiles = Table.read("/data/desardata/" + release + "/InfoTables/" + release + "_Tilename_Run_COADD.fits")
	tile_ras = []
	tile_decs = []

	for tile in tiles["TILENAME"]:
		ra = float(tile[3:5])/24.0*360.0 + float(tile[5:7])/(24.0*60.0)*360.0
		if tile[7] == "-":
			dec = float(tile[7:10]) - float(tile[10:])/60.0
		else:
			dec = float(tile[7:10]) + float(tile[10:])/60.0
		tile_ras.append(ra)
		tile_decs.append(dec)

	dists, inds = match_lists.match_lists(ras, decs, tile_ras, tile_decs, 0.5, numNei = 5)

	l = len(list(set(inds[inds < len(tile_ras)].flatten())))

	ids = np.where( (inds <> len(tile_ras)) )[0]
	l = len(ids)

	t = Table([["DES0000-0000"]*l, ["00000000000000_DES0000-0000"]*l, [0.000000]*l, [0.000000]*l], names = ("TILENAME", "RUN", "RA_IN", "DEC_IN"))

	n = 0
	i = 0
	ids = np.where( (inds <> len(tile_ras)) )[0]
	while n < len(ras):
		ra = ras[n]
		dec = decs[n]
		for (dist, ind) in zip(dists[n], inds[n]):
			if ind <> len(tile_ras):
				t["RA_IN"][i] = ra
				t["DEC_IN"][i] = dec
				t["TILENAME"][i] = tiles["TILENAME"][ind]
				t["RUN"][i] = tiles["RUN"][ind]
				i += 1
		n += 1

	return t

def DES_data(t, release = "SVA1"):

	"""
	Takes a table produced by DES_tile_run of ra, dec and possible TILENAME and RUN
	Uses this to then search each tile for the potential DES match to the object.
	Puts matching objects into a new table
	"""

	import os
	import subprocess
	import match_lists
	from astropy.table import Table, vstack

	if release == "SVA1":
		rn = 0
	elif release == "Y1P1":
		rn = 1
	elif release == "Y1A1":
		rn = 2

	n = 0
	m = 0
	for (ra, dec, tile, run) in zip(t["RA_IN"], t["DEC_IN"], t["TILENAME"], t["RUN"]):
		print m, "out of", len(t)
		if tile <> "DES0000-0000":
			print tile
			if os.path.exists("/data/desardata/" + release + "/" + tile + "/" + tile + ".fits"):
				t_info = Table.read("/data/desardata/" + release + "/" + tile + "/" + tile + ".fits")
			else:
				subprocess.call(["bash", "/home/sr525/bash_scripts/tile_file.bash", release, tile])
				t_info = Table.read("/data/desardata/" + release + "/" + tile + "/" + tile + ".fits")

			dists, inds = match_lists.match_lists([ra], [dec], t_info["RA"], t_info["DEC"], 2.0/3600)
			print dists, inds
			if inds[0] <> len(t_info):
				if n == 0:
					t1 = Table(t_info[inds[0]])
					n += 1
				else:
					if inds[0] <> len(t_info):
						t_info = Table(t_info[inds[0]])
						t1 = vstack([t1, t_info], join_type = "exact")
		m += 1

	try:
		return t1
	except UnboundLocalError as e:
		print e
		print "Returning Nothing"
		return Table()

def DES_Y2_image_file(RA, DEC, bands = ["g", "r", "i", "z", "Y"], download = True):

	"""
	Returns a list of image files and a list of bands for the DES Y2 single epoch 
	images from Y2N1 with a given RA and DEC on them
	Some of these are corrupted and objects can be near/on the edge
	If the object is near/on the edge the image will not be square and the object 
	may not be central
	Can specify just a few bands
	Downloads the images then funzips them, then gzips them.
	"""

	from astropy.table import Table, vstack, hstack
	import numpy as np
	import os
	import subprocess

	info = Table.read("/data/desardata/Y2N/FIRSTCUT/files.fits")

	ids = np.where( ((RA > info["RAC1"]) | (RA > info["RAC2"])) & \
					((RA < info["RAC3"]) | (RA < info["RAC4"])) & \
					((DEC < info["DECC1"]) | (DEC < info["DECC4"])) & \
					((DEC > info["DECC2"]) | (DEC > info["DECC3"])) & \
					(info["RAC1"] <> 0.0) & (info["RAC3"] <> 0.0) )[0]

	info = info[ids]
	files = []
	out_bands = []
	n = 0
	while n < len(info):
		nite = str(info["NITE"][n])
		uname = str(info["UNITNAME"][n][:-11])
		rnum = str(info["REQNUM"][n])
		anum = str(info["ATTNUM"][n])
		if info["CCDNUM"][n] > 9:
			cnum = str(info["CCDNUM"][n])
		else:
			cnum = "0" + str(info["CCDNUM"][n])
		band = str(info["BAND"][n][:-4])

		file = "/data/desardata/Y2N/FIRSTCUT/" + nite + "/" + uname + "/p0" + anum + "/red/" + uname + "_" + band + "_c" + cnum + "_r"+ rnum + "p0" + anum + "_immasked.fits"
		#print file
		if os.path.exists(file + ".gz") and not os.path.exists(file + ".fz"):
			subprocess.call(["gunzip", file + ".gz"])
			subprocess.call(["fpack", file])
		elif os.path.exists(file + ".gz") and os.path.exists(file + ".fz"):
			subprocess.call(["rm", "-f", file + ".gz"])
		if (not os.path.exists(file) and not os.path.exists(file + ".fz")) \
			and download and band in bands:
			#print "Trying to download"
			subprocess.call(["bash", "/home/sr525/bash_scripts/DES_Y2_download.bash", nite, rnum, uname, anum, band, cnum])
		#if not os.path.exists(file + ".gz") and os.path.exists(file + ".fz") and not os.path.exists(file):
		#	try:
		#		subprocess.check_call(["funpack", file + ".fz"])
		#		subprocess.check_call(["gzip", "-f", file])
		#	except subprocess.CalledProcessError:
		#		print "Corrupt File Skipping"
		#if os.path.exists(file):
		#	try:
		#		subprocess.check_call(["gzip", "-f", file])
		#	except subprocess.CalledProcessError:
		#		print "Corrupt File Skipping"
		#if os.path.exists(file + ".gz") and band in bands:
		#	files.append(file + ".gz")
		#	out_bands.append(band)
		if not os.path.exists(file + ".fz") and os.path.exists(file):
			try:
				subprocess.check_call(["fpack", file + ".fz"])
			except subprocess.CalledProcessError:
				print "Corrupt File Skipping"
		if band in bands:
			files.append(file + ".fz")
			out_bands.append(band)

		n += 1

	return files, out_bands

def flux_match(ra, dec, t_all, RA_main, DEC_main, l, width, RA_all, DEC_all, bands):

	"""
	Does something with fluxes and I'm not sure what or why I wrote this
	"""

	from astropy.table import Table
	import match_lists

	t_all.sort([RA_all, DEC_all])
	dist, ind = match_lists.match_lists(ra, dec, t_all[RA_all], t_all[DEC_all], width)
	i = ind[0]
	inds = range(i-l, i+l)
	t_short = t_all[inds]
	mc = 0
	mags = []
	fs = []
	ffs = []
	fmags = []
	for band in bands:
		mag = t_all["MAG_APER_5_" + band][i]
		flux = t_all["FLUX_APER_5_" + band][i]
		mc += 1
		mags.append(mag)
		fs.append(flux)
		if mc == 5:
			fmags = mags
			ffs = fs
	return ffs, fmags

def inside_contour(xs, ys, xs1, ys1, c_graph = True, radius = -0.1, nbins = 1000):

	"""
	Takes two sets of points and sees if one is inside the other
	Might need nbins and radius tweaking
	xs, ys are the points to draw the contours off
	xs1, ys1 are the points to check if in contour
	"""

	import numpy as np
	import matplotlib.pyplot as plt

	xwalls = np.linspace( min(xs) - 5.0, max(xs) + 5.0, nbins + 1 )
	ywalls = np.linspace( min(ys) - 5.0, max(ys) + 5.0, nbins + 1 )

	im, xs_bin, ys_bin, ax = plt.hist2d(xs, ys, bins = (xwalls, ywalls) )
	xs_mids = 0.5*(xs_bin[:-1] + xs_bin[1:])
	ys_mids = 0.5*(ys_bin[:-1] + ys_bin[1:])
	plt.close()
	im[im>0] = 1
	if c_graph:
		#plt.plot(xs[::1000], ys[::1000], "k.", ms = 1)
		plt.plot(xs, ys, "k.", ms = 1)
	conts = plt.contour(xs_mids, ys_mids, im.T, 1)
	if c_graph:
		plt.show()
		#plt.plot(xs1[::1000], ys1[::1000], "k.", ms = 1)
		plt.plot(xs1, ys1, "k.", ms = 1)
		plt.contour(xs_mids, ys_mids, im.T, 1)
		plt.show()
	else:
		plt.close()
	paths = conts.collections[0]
	paths = paths.get_paths()

	z = np.array( [xs1, ys1] ).T
	inside = paths[0].contains_points(z, radius = radius)
	counter = 0
	for p in paths[1:]:
		print counter
		counter += 1
		inside_new = p.contains_points(z, radius = radius)
		inside = inside | inside_new

	#xs1_out = xs1[inside]
	#ys1_out = ys1[inside]

	return inside

def survey_overlap(xs, ys, xs1, ys1, c_graph = True, radius = -0.1, nbins = 1000):

	"""
	Takes two sets of points and sees if one is inside the other
	Might need nbins and radius tweaking
	xs, ys are the points to draw the contours off
	xs1, ys1 are the points to check if in contour
	"""

	import numpy as np
	import matplotlib.pyplot as plt

	xwalls = np.linspace( min(xs) - 5.0, max(xs) + 5.0, nbins + 1 )
	ywalls = np.linspace( min(ys) - 5.0, max(ys) + 5.0, nbins + 1 )

	im, xs_bin, ys_bin, ax = plt.hist2d(xs, ys, bins = (xwalls, ywalls) )
	xs_mids = 0.5*(xs_bin[:-1] + xs_bin[1:])
	ys_mids = 0.5*(ys_bin[:-1] + ys_bin[1:])
	plt.close()
	im[im>0] = 1
	conts = plt.contour(xs_mids, ys_mids, im.T, 1)
	
	xwalls1 = np.linspace( min(xs1) - 5.0, max(xs1) + 5.0, nbins + 1)
	ywalls1 = np.linspace( min(ys1) - 5.0, max(ys1) + 5.0, nbins + 1)
	im1, xs_bin1, ys_bin1, ax1 = plt.hist2d(xs1, ys1, bins = (xwalls1, ywalls1))
	xs_mids1 = 0.5*(xs_bin1[:-1] + xs_bin1[1:])
	ys_mids1 = 0.5*(ys_bin1[:-1] + ys_bin1[1:])
	plt.close()
	im1[im1>0] = 1

	conts = plt.contour(xs_mids, ys_mids, im.T, 1)
	conts1 = plt.contour(xs_mids1, ys_mids1, im1.T, 1)
	plt.show()

def DES_tile(t, release = "Y1A1", ra_col = "ra", dec_col = "dec"):

	"""
	Takes a table and finds out what DES tile each object is on
	"""

	import atpy
	from astropy.table import Table, vstack, hstack
	import astropy.io.fits as fits
	import numpy as np

	t_tiles = atpy.Table("/data/desardata/"+ release+"/InfoTables/" + release + "_Tilename_Run_COADD.fits")
	tiles = t_tiles["TILENAME"]
	bad_tiles = ["DES0001-5705", "DES0000-5248", "DES0319-6456"]
	n = 0
	m = 0
	for tile in tiles:
		if tile not in bad_tiles:
			print m+1, "out of", len(t_tiles), tile
			with fits.open("/data/desardata/Y1A1/" + tile + "/" + tile + "_z.fits.fz") as fhlist:
				hdr = fhlist[1].header
				from astropy import wcs
				w = wcs.WCS(hdr, naxis = 2)
				pix_coord = [[0.0,0.0], [hdr["NAXIS1"],hdr["NAXIS2"]]]
				w_coord = w.wcs_pix2world(pix_coord, 1)
				ra_min = min([w_coord[0][0], w_coord[1][0]])
				ra_max = max([w_coord[0][0], w_coord[1][0]])
				dec_min = min([w_coord[0][1], w_coord[1][1]])
				dec_max = max([w_coord[0][1], w_coord[1][1]])
				print ra_min, ra_max, dec_min, dec_max
				if ra_max > ra_min + 10.0:
					ra_min, ra_max = ra_max, ra_min
	
				ids = np.where( (t[ra_col] > ra_min) & (t[ra_col] < ra_max) & (t[dec_col] > dec_min) & (t[dec_col] < dec_max) )[0]

				t1 = t[ids]
				t1["TILENAME_DES"] = [tile]*len(t1)
				if n == 0 and len(t1) > 0.0:
					t_out = t1
					n += 1
				else:
					if len(t1) > 0.0:
						t_out = vstack([t_out, t1])
						print len(t_out)
						n += 1
			m += 1

	return t_out

def DES_tile_full_table(t, release = "Y1A1", ra_col = "ra", dec_col = "dec"):

	"""
	Takes a table and finds out what DES tile each object is on
	"""

	import atpy
	from astropy.table import Table, vstack, hstack
	import astropy.io.fits as fits
	import numpy as np

	t_tiles = atpy.Table("/data/desardata/"+ release+"/InfoTables/" + release + "_Tilename_Run_COADD.fits")
	tiles = t_tiles["TILENAME"]
	bad_tiles = ["DES0001-5705", "DES0000-5248", "DES0319-6456"]
	n = 0
	m = 0
	t["TILENAME"] = ["DES0000-0000"]*len(t)
	for tile in tiles:
		if tile not in bad_tiles:
			print m+1, "out of", len(t_tiles), tile
			with fits.open("/data/desardata/Y1A1/" + tile + "/" + tile + "_z.fits.fz") as fhlist:
				hdr = fhlist[1].header
				from astropy import wcs
				w = wcs.WCS(hdr, naxis = 2)
				pix_coord = [[0.0,0.0], [hdr["NAXIS1"],hdr["NAXIS2"]]]
				w_coord = w.wcs_pix2world(pix_coord, 1)
				ra_min = min([w_coord[0][0], w_coord[1][0]])
				ra_max = max([w_coord[0][0], w_coord[1][0]])
				dec_min = min([w_coord[0][1], w_coord[1][1]])
				dec_max = max([w_coord[0][1], w_coord[1][1]])
				print ra_min, ra_max, dec_min, dec_max
				if ra_max > ra_min + 10.0:
					ra_min, ra_max = ra_max, ra_min
	
				ids = np.where( (t[ra_col] > ra_min) & (t[ra_col] < ra_max) & (t[dec_col] > dec_min) & (t[dec_col] < dec_max) )[0]

				t["TILENAME"][ids] = tile
		m += 1

	return t

def get_apm_cat(RA, DEC, db = "vst", aper_num = "3"):

	"""
	Finds the APM catalogue for a given RA and DEC
	"""

	import cutout_vista
	import xmlrpclib
	import os
	import subprocess
	import numpy as np
	from astropy.table import Table, vstack, hstack

	try:
		files, apm_files = cutout_vista.querydb(RA, DEC, db=db, onlytiles=True, ctype='fit')
		if files == []:
			"Print no files returned"

	except xmlrpclib.Fault as e:
		if e.faultCode == 1:
			#faultCode 1 is IOError
			print "Data not found"
		else:
			raise e

	n = 0
	for filename in apm_files:
		dir, file = os.path.split(filename[:-6])
		file = "/" + file + "cat.fits"
		subprocess.call(["/home/sr525/bash_scripts/scp_amp.bash", dir, file])
		with fits.open("/tmp/" + file) as hlist:
			th = hlist[1].header
			td = hlist[1].data
			td["RA"] = np.rad2deg(td["RA"])
			td["DEC"] = np.rad2deg(td["DEC"])
			dists, inds = match_lists.match_lists([RA], [DEC], td["RA"], td["DEC"], 10.1)
			print dists, inds
			if "177.A-3011" in hlist[0].header["HIERARCH ESO OBS PROG ID"]:
				band = files[n][-8:-7]
				print band
				apcor = th["APCOR" + aper_num]
				
				#pho_cal.
			else:
				print "Can't use"
				return
		n += 1

		os.remove("/tmp" + file)
	for file in files:
		os.remove(file)

def Y2_bands(RA, DEC, bands = ["g", "r", "i", "z", "Y"]):

	files, out_bands = DES_Y2_image_file(RA, DEC, bands = ["g", "r", "i", "z", "Y"], download = False)

	if len(set(out_bands)) >= len(bands):
		return True

	else:
		return False

def Y2_bands_Table(t, bands = ["g", "r", "i", "z", "Y"]):

	n = 0
	keeps = []
	while n < len(t):
		print n+1, "out of", len(t)
		RA = t["RA"][n]
		DEC = t["DEC"][n]
		all_bands = Y2_bands(RA, DEC, bands = bands)
		if all_bands:
			keeps.append(n)
		n += 1

	return t[keeps]

def DES_name(RA, DEC, length = "short"):

	from astropy import coordinates as coord
	from astropy import units as u
	from astropy.coordinates import SkyCoord
	import numpy as np

	c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='icrs')
	if length == "short":
		if c.dec.dms[0] < 0.0:
			name = "VDES%02d%02d%02d%02d" % (c.ra.hms[0], c.ra.hms[1], c.dec.dms[0], np.fabs(c.dec.dms[1]))
		else:
			name = "VDES%02d%02d+%02d%02d" % (c.ra.hms[0], c.ra.hms[1], c.dec.dms[0], np.fabs(c.dec.dms[1]))

	elif length == "long":
		if c.dec.dms[0] < 0.0:
			name = "VDES%02d%02d%05.2f%02d%02d%05.2f" % (c.ra.hms[0], c.ra.hms[1], c.ra.hms[2], c.dec.dms[0], np.fabs(c.dec.dms[1]), np.fabs(c.dec.dms[2]))
		else:
			name = "VDESJ%02d%02d%05.2f+%02d%02d%05.2f" % (c.ra.hms[0], c.ra.hms[1], c.ra.hms[2], c.dec.dms[0], np.fabs(c.dec.dms[1]), np.fabs(c.dec.dms[2]))

	return name

def get_DES_cat(tile, release = "Y1A1"):

	from astropy.table import Table
	import numpy as np
	import subprocess

	t_info = Table.read("/data/desardata/" + release + "/InfoTables/" + release + "_Tilename_Run_COADD.fits")

	ids = np.where( (t_info["TILENAME"] == tile) )

	runs = t_info["RUN"][ids].data
	
	for run in runs:
		run = run[:-3]
		subprocess.call(["bash", "/home/sr525/bash_scripts/wget_cat.bash", tile, run, release])

	return runs

def unWISE_list(t):

	from astropy.table import Table

	tw = Table.read("/data/rgm/wise/wise_allwise_metadata_thin.fits")

	n = 0
	files = []

	while n < len(t):

		RA = t["ALPHAWIN_J2000_Z"][n]
		DEC = t["DELTAWIN_J2000_Z"][n]

		m = 0
		while m < len(tw):
			ras = sorted([tw["ra1"][n], tw["ra2"][n], tw["ra3"][n], tw["ra4"][n]])
			decs = sorted([tw["dec1"][n], tw["dec2"][n], tw["dec3"][n], tw["dec4"][n]])
			min_ra = ras[1]
			max_ra = ras[2]
			min_dec = decs[1]
			max_dec = decs[2]
			if max_ra > min_ra + 10.0:
				min_ra, max_ra = max_ra, min_ra
			if DEC > min_dec and DEC < max_dec and RA > min_ra and RA < max_ra:
				tile_id = tw["coadd_id"][n][0:8]
				f2 = tile_id[0:3]
				wise_dir = "/data/wiseardata/" + release + "/p3am_cdd/" + f2 + "/" 
				wise_file = wise_dir + tile_id + "/unwise-" + tile_id + "-"
			files.append(wise_file)
			m += 1
		n += 1

	return files
