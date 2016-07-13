"""
Matching functions for various data
Some needs tidying up but it works so I'm ignoring it
None of the functions change the magnitudes from their native systems
"""
def check_graph(t, RA_main, DEC_main, survey, width, gtype = "all"):

	"""
	Makes a diagnostic plot for checking matches.
	Plot can either be square, the square inside the circle it matches it. Or all 
	which has all the points in the matching circle. Square makes the histograms 
	more comparable.
	Compares RA_main and DEC_main columns with RA and Dec columns in the format 
	output by the matching codes. Eg. RA_ + survey.
	Width needs to be in arcsecs
	"""

	import math
	import stats
	import matplotlib.pyplot as plt
	import matplotlib.gridspec as gridspec
	import numpy as np
	from matplotlib.colors import LogNorm

	n = 0
	xs = []
	ys = []
	while n < len(t):
		x = (t[RA_main][n] - t["RA_" + survey.upper()][n])*math.cos((t[DEC_main][n]+t["DEC_" + survey][n])*math.pi/360.0)*3600.0
		y = (t[DEC_main][n] - t["DEC_" + survey.upper()][n])*3600.0

		if not np.isnan(x) and not np.isnan(y):
			xs.append(x)
			ys.append(y)
		n += 1

	n = 0
	xs_s = []
	ys_s = []
	if gtype == "square":
		w = width / math.sqrt(2.0)
		while n < len(xs):
			x = xs[n]
			y = ys[n]
			if x <= w and x >= -w and y <= w and y >= -w:
				xs_s.append(xs[n])
				ys_s.append(ys[n])
			n += 1

		xs = xs_s
		ys = ys_s

	xs1 = list(xs) + []
	ys1 = list(ys) + []

	RA_med = np.median(xs1)
	DEC_med = np.median(ys1)
	RA_MAD = stats.MAD(xs1, RA_med)
	DEC_MAD = stats.MAD(ys1, DEC_med)
	print "Number of points", len(xs)
	print "RA offset", RA_med, "DEC offset", DEC_med
	print "RA MAD", RA_MAD, "DEC MAD", DEC_MAD
	print "RA Sigma MAD", 1.486*RA_MAD, "DEC Sigma DEC", 1.486*DEC_MAD
	print "RA Median Error", 1.486*RA_MAD/math.sqrt(len(xs)), "DEC Median Error", 1.486*DEC_MAD/math.sqrt(len(ys))
	if len(xs) == 0:
		print "No matches"
		return RA_med, DEC_med
	gs = gridspec.GridSpec(2,2, width_ratios = [2,1], height_ratios = [1,2])
	fig = plt.figure()
	ax1 = plt.subplot(gs[0])
	ax1.hist(xs, bins = 100, color = "r")
	ax1.set_xlim(-2.0, 2.0)
	ax1.axes.get_xaxis().set_visible(False)
	ax1.set_ylabel("Number")
	
	ax2 = plt.subplot(gs[2])
	#ax2.plot(xs, ys, "k+")
	if len(xs) > 100:
		plt.hist2d(xs, ys, bins = 100, cmap = "binary", norm = LogNorm())
	else:
		plt.plot(xs, ys, "k.", ms = 2)
	ax2.set_ylim(-2.0, 2.0)
	ax2.set_xlim(-2.0, 2.0)
	ax2.set_xlabel('Delta RA /"')
	ax2.set_ylabel('Delta Dec /"')
	labels1 = ax2.get_xticks()
	ax2.set_xticklabels(labels1, rotation = 270)
	
	fig.suptitle("Errors in matching between DES and " + survey)
	
	ax3 = plt.subplot(gs[3]) 
	ax3.hist(ys, bins = 100, orientation = "horizontal", color = "r")
	ax3.set_ylim(-2.0, 2.0)
	ax3.set_xlabel("Number")
	ax3.axes.get_yaxis().set_visible(False)
	labels2 = ax3.get_xticks()
	ax3.set_xticklabels(labels2, rotation = 270)

	ax4 = plt.subplot(gs[1])
	ax4.annotate("Number of points: " + str(len(xs)), xy = (0.01, 0.1), size = "small")
	ax4.annotate("RA offset: {0:.4f}".format(RA_med) + '"', xy = (0.01,0.90), size = "small")
	ax4.annotate("DEC offset: {0:.4f}".format(DEC_med) + '"', xy = (0.01, 0.8), size = "small")
	ax4.annotate("RA MAD: {0:.4f}".format(RA_MAD) + '"', xy = (0.01, 0.7), size = "small")
	ax4.annotate("DEC MAD: {0:.4f}".format(DEC_MAD) + '"', xy = (0.01, 0.6), size = "small")
	ax4.annotate("RA median error: {0:.4f}".format(1.486*RA_MAD/math.sqrt(len(xs))) + '"', xy = (0.01, 0.5), size = "small")
	ax4.annotate("DEC median error: {0:.4f}".format(1.486*DEC_MAD/math.sqrt(len(ys))) + '"', xy = (0.01, 0.4), size = "small")
	ax4.annotate("RA sigma MAD: {0:.4f}".format(RA_MAD*1.486) + '"', xy = (0.01, 0.3), size = "small")
	ax4.annotate("DEC sigma MAD: {0:.4f}".format(DEC_MAD*1.486) + '"', xy = (0.01, 0.2), size = "small")

	ax4.axes.get_xaxis().set_visible(False)
	ax4.axes.get_yaxis().set_visible(False)

	plt.show()

	return RA_med, DEC_med

def min_max_match(ras, decs, width, DES_box = False):

	"""
	Finds the min and max ra and dec of the input points
	This governs the corners of the box the SQL downloads
	DES_box downloads an area about the size of a DES tile centred on the median
	ra and median dec
	"""

	import numpy as np

	print "Finding max and mins"
	ra_min = str(min(ras) - width)
	ra_max = str(max(ras) + width)
	dec_min = str(min(decs) - width)
	dec_max = str(max(decs) + width)
	if float(ra_min) < 3.0 and float(ra_max) > 357.0:
		ra_min = str(0.0)
		ra_max = str(360.0)
	if DES_box:
		ra_min = str(np.median(ras) - 0.7)
		ra_max = str(np.median(ras) + 0.7)
		dec_min = str(np.median(decs) - 0.7)
		dec_max = str(np.median(decs) + 0.7)

	return ra_min, ra_max, dec_min, dec_max

def VIKING_match(t, RA_main, DEC_main, width = 0.0004, w_units = "degrees", c_graph = True, DES_box = False):

	"""
	Matches to Viking Data from WSDB
	Table used is viking_201207.main
	DES_box returns all objects within a box of about the size of a DES tile centred
	on the median RA and DEC of the table.
	Iterates through a table where the RA column is RA_main and the DEC_column is 
	DEC_main
	Puts NaNs in for no matches
	c_graph plots diagnostic plot of the match - crashes if there are no matches
	"""

	import numpy as np
	import match_lists
	import sqlutil
	import ConfigParser

	config_file = "/home/sr525/Python_Code/srpylib/srpylib.cfg"
	config = ConfigParser.RawConfigParser()
	config.read(config_file)
	uname = config.get("wsdb", "uname")
	pword = config.get("wsdb", "pword")
	host = config.get("wsdb", "host")

	if w_units == "arcsecs":
		width = width / 3600.0

	cols = ["RA_VIKING", "DEC_VIKING", "VIKING_Y", "VIKING_J"]

	try:
		for col in cols:
			t.add_empty_column(col, dtype = np.float64)

	except AttributeError as e:
		for col in cols:
			t[col] = np.zeros(len(t))

	n = 0

	if len(t) > 1000 and len(t) < 1000000:

		ra_min, ra_max, dec_min, dec_max = min_max_match(t[RA_main], t[DEC_main], \
												width, DES_box = DES_box)

		print "Retreiving data"
		query_V = "select ra, dec, mag_y, mag_j FROM viking_201207.main WHERE ra > "\
					+ ra_min + " and ra < " + ra_max + " and dec > " + \
					dec_min + " and dec < " + dec_max
		RA_V, DEC_V, Y_V, J_V = sqlutil.get(query_V, db = "wsdb", host = host,\
								user = uname, password = pword)

		if len(RA_wise) == 0:
			for (i, col) in enumerate(cols):
				t[col] = [np.float64(np.nan)]*len(t)

		else:
			info = [RA_V, DEC_V, Y_V, J_V]

			print "Retreived data, matching lists"
			dists, inds = match_lists.match_lists(t[RA_main], t[DEC_main], \
													RA_V, DEC_V, width, 1)
			print "Matched lists, updating table"

			match = np.where( (inds <> len(RA_V)) )[0]
			no_match = np.where( (inds == len(RA_V)) )[0]

			for (i, col) in enumerate(cols):
				t[col][match] = info[i][inds[match]]
				t[col][no_match] = np.float64(np.nan)
			print "Table updated"

	elif len(t) <= 1000:
		while n < len(t):
			RA = t[RA_main][n]
			DEC = t[DEC_main][n]
			query = "select ra, dec, mag_y, mag_j FROM viking_201207.main WHERE " + \
						"q3c_radial_query(ra, dec, " + str(RA) + \
						", " + str(DEC) + ", 0.1002)"
			RA_V, DEC_V, Y_V, J_V = sqlutil.get(query, db='wsdb',host=host, \
											user = uname, password=pword)

			if len(RA_V) > 0:
				dists, inds = match_lists.match_lists([RA], [DEC], RA_V, DEC_V, width, 1)
				j = inds[0]
				info = [RA_V, DEC_V, Y_V, J_V]

				for (i, col) in enumerate(cols):
					t[col][n] = info[i]

			else:
				for (i, col) in enumerate(cols):
					t[col][n] = np.float64(np.nan)

			n += 1

	elif len(t) > 1000000:
	
		print "Splitting by tile"
		tiles_all = t["TILENAME"]
		tiles = set(tiles_all)
		tiles = list(tiles)
		print len(tiles), "distinct tiles"

		t.sort("TILENAME")
		tiles.sort()

		n = 0
		m = 0
		rows = []
		tile_rows = []

		while n < len(t):
			tile = tiles[m]
			if t["TILENAME"][n] == tile:
				rows.append(n)
				n += 1
			else:
				tile_rows.append(rows)
				rows = []
				m += 1

		for rows in tile_rows:
			t1 = t.rows(rows)
			print "Matching", t1["TILENAME"][0], "with length", len(t1)
			t1 = VIKING_match(t1, RA_main, DEC_main, width, c_graph = False)
			print "Returning table of length:", len(t1)

			if tile_rows.index(rows) == 0:
				t_out = t1
			else:
				t_out.append(t1)
		t = t_out			

	if c_graph:
		width = width *3600.0
		RA_med, DEC_med = check_graph(t, RA_main, DEC_main, "VIKING", width, gtype = "square")
			
	return t

def SDSS_match(t, RA_main, DEC_main, width = 0.0004, w_units = "degrees", c_graph = True, dr = "dr9", DES_box = False):

	"""
	Matches to SDSS Data from WSDB
	Table used is dr + .phototag
	DES_box returns all objects within a box of about the size of a DES tile centred
	on the median RA and DEC of the table.
	Iterates through a table where the RA column is RA_main and the DEC_column is 
	DEC_main
	Puts NaNs in for no matches
	c_graph plots diagnostic graphs of the match - crashes if there are no matches
	"""

	import numpy as np
	import match_lists
	import sqlutil
	from astropy.table import Table, vstack
	import ConfigParser

	config_file = "/home/sr525/Python_Code/srpylib/srpylib.cfg"
	config = ConfigParser.RawConfigParser()
	config.read(config_file)
	uname = config.get("wsdb", "uname")
	pword = config.get("wsdb", "pword")
	host = config.get("wsdb", "host")

	if w_units == "arcsecs":
		width = width / 3600.0

	cols = ["RA_SDSS", "DEC_SDSS", "MAG_PSF_U_SDSS", "MAGERR_PSF_U_SDSS", \
			"MAG_PSF_G_SDSS", "MAGERR_PSF_G_SDSS", "MAG_PSF_R_SDSS", \
			"MAGERR_PSF_R_SDSS", "MAG_PSF_I_SDSS", "MAGERR_PSF_I_SDSS", \
			"MAG_PSF_Z_SDSS", "MAGERR_PSF_Z_SDSS", "MAG_MODEL_U_SDSS", \
			"MAGERR_MODEL_U_SDSS", "MAG_MODEL_G_SDSS", "MAGERR_MODEL_G_SDSS", \
			"MAG_MODEL_R_SDSS", "MAGERR_MODEL_R_SDSS", "MAG_MODEL_I_SDSS", \
			"MAGERR_MODEL_I_SDSS", "MAG_MODEL_Z_SDSS", "MAGERR_MODEL_Z_SDSS", "TYPE"]

	try:
		for col in cols:
			t.add_empty_column(col, dtype = np.float64)

	except AttributeError as e:
		for col in cols:
			t[col] = np.zeros(len(t))

	n = 0

	if len(t) > 100 and len(t) < 1000000:
		print "Finding max and mins"
		ra_min, ra_max, dec_min, dec_max = min_max_match(t[RA_main], t[DEC_main], \
														width, DES_box = DES_box)

		print "Retreiving data"

		query_sdss = "select ra, dec, psfMag_u, psfMagErr_u, psfMag_g, " + \
					"psfMagErr_g, psfMag_r, psfMagErr_r, psfMag_i, psfMagErr_i, " + \
					"psfMag_z, psfMagErr_z, modelMag_u, modelMagErr_u, " + \
					"modelMag_g, modelMagErr_g, modelMag_r, modelMagErr_r, " + \
					"modelMag_i, modelMagErr_i, modelMag_z, modelMagErr_z, " + \
					"type FROM sdss" + dr + ".phototag WHERE ra > " + ra_min + \
					" and ra < " + ra_max + " and dec > " + dec_min + \
					" and dec < " + dec_max

		RA_sdss, DEC_sdss, psfu, psfuerr, psfg, psfgerr, psfr, psfrerr, psfi, \
		psfierr, psfz, psfzerr, modelu, modeluerr, modelg, modelgerr, modelr, \
		modelrerr, modeli, modelierr, modelz, modelzerr, types = \
		sqlutil.get(query_sdss, db = "wsdb", host = host, user = uname, password = pword)

		info = [RA_sdss, DEC_sdss, psfu, psfuerr, psfg, psfgerr, psfr, psfrerr, \
				psfi, psfierr, psfz, psfzerr, modelu, modeluerr, modelg, modelgerr, \
				modelr, modelrerr, modeli, modelierr, modelz, modelzerr, types]

		if len(RA_sdss) == 0:
			for (i, col) in enumerate(cols):
				t[col] = [np.float64(np.nan)]*len(t)

		else:
			print "Retreived data, matching lists"
			dists, inds = match_lists.match_lists(t[RA_main], t[DEC_main], RA_sdss, DEC_sdss, width, 1)
			print "Matched lists, updating table"

			match = np.where( (inds <> len(RA_sdss)) )[0]
			no_match = np.where( (inds == len(RA_sdss)) )[0]

			for (i, col) in enumerate(cols):
				t[col][match] = info[i][inds[match]]
				t[col][no_match] = np.float64(np.nan)

			print "Table updated"

	elif len(t) <= 100:
		while n < len(t):
			RA = t[RA_main][n]
			DEC = t[DEC_main][n]
			print RA, DEC
			query_sdss = "select ra, dec, psfMag_u, psfMagErr_u, psfMag_g, " + \
						"psfMagErr_g, psfMag_r, psfMagErr_r, psfMag_i, " + \
						"psfMagErr_i, psfMag_z, psfMagErr_z, modelMag_u, " + \
						"modelMagErr_u, modelMag_g, modelMagErr_g, modelMag_r, " + \
						"modelMagErr_r, modelMag_i, modelMagErr_i, modelMag_z, " + \
						"modelMagErr_z, type FROM sdss" + dr + \
						".phototag WHERE q3c_radial_query(ra, dec, " + str(RA) + \
						", " + str(DEC) + ", 0.1002)"

			RA_sdss, DEC_sdss, psfu, psfuerr, psfg, psfgerr, psfr, psfrerr, psfi, \
			psfierr, psfz, psfzerr, modelu, modeluerr, modelg, modelgerr, modelr, \
			modelrerr, modeli, modelierr, modelz, modelzerr, stype = \
			sqlutil.get(query_sdss, db = "wsdb", host = host, user = uname, password = pword)

			if len(RA_sdss) > 0:
				dist, inds = match_lists.match_lists([RA], [DEC], RA_sdss, DEC_sdss, width, 1)
				j = inds[0]
				info = [RA_sdss[j], DEC_sdss[j], psfu[j], psfuerr[j], psfg[j], \
						psfgerr[j], psfr[j], psfrerr[j], psfi[j], psfierr[j], \
						psfz[j], psfzerr[j], modelu[j], modeluerr[j], modelg[j], \
						modelgerr[j], modelr[j], modelrerr[j], modeli[j], \
						modelierr[j], modelz[j], modelzerr[j], stype[j]]

				for (i, col) in enumerate(cols):
					t[col][n] = info[i]

			else:
				for (i, col) in enumerate(cols):
					t[col][n] = np.float64(np.nan)

			n += 1


	elif len(t) > 1000000:
	
		print "Splitting by tile"
		tiles_all = t["TILENAME"]
		tiles = set(tiles_all)
		tiles = list(tiles)
		print len(tiles), "distinct tiles"

		t.sort("TILENAME")
		tiles.sort()

		n = 0
		m = 0
		rows = []
		tile_rows = []

		while n < len(t):
			tile = tiles[m]
			if t["TILENAME"][n] == tile:
				rows.append(n)
				n += 1
			else:
				print len(rows)
				tile_rows.append(rows)
				rows = []
				m += 1

		for rows in tile_rows:
			try:
				t1 = t.rows(rows)
			except AttributeError:
				t1 = t[rows]
			print "Matching", t1["TILENAME"][0], "with length", len(t1)
			t1 = SDSS_match(t1, RA_main, DEC_main, width, c_graph = False)

			if tile_rows.index(rows) == 0:
				t_out = t1
			else:
				try:
					t_out.append(t1)
				except AttributeError:
					t_out = vstack([t_out, t1])
		t = t_out			

	if c_graph:
		width = width *3600.0
		RA_med, DEC_med = check_graph(t, RA_main, DEC_main, "SDSS", width, gtype = "square")
			
	return t

def WISE_match(t, RA_main, DEC_main, width = 0.0004, w_units = "degrees", c_graph = True, DES_box = False):

	"""
	Matches to WISE Data from WSDB
	Table used is allwise.main
	DES_box returns all objects within a box of about the size of a DES tile centred
	on the median RA and DEC of the table.
	Iterates through a table where the RA column is RA_main and the DEC_column is 
	DEC_main
	Puts NaNs in for no matches
	c_graph plots diagnostic graphs of the match - crashes if there are no matches
	"""

	import numpy as np
	import match_lists
	import sqlutil
	import ConfigParser

	config_file = "/home/sr525/Python_Code/srpylib/srpylib.cfg"
	config = ConfigParser.RawConfigParser()
	config.read(config_file)
	uname = config.get("wsdb", "uname")
	pword = config.get("wsdb", "pword")
	host = config.get("wsdb", "host")

	if w_units == "arcsecs":
		width = width / 3600.0

	cols = ["RA_WISE", "DEC_WISE", "W1", "W1_ERR", "W2", "W2_ERR", "W3", "W3_ERR", "W4", "W4_ERR"]

	try:
		for col in cols:
			t.add_empty_column(col, dtype = np.float64)

	except AttributeError as e:
		for col in cols:
			t[col] = np.zeros(len(t))

	n = 0

	if len(t) > 5000 and len(t) < 100000000:

		ra_min, ra_max, dec_min, dec_max = min_max_match(t[RA_main], t[DEC_main], \
														width, DES_box = DES_box)

		print "Retreiving data"
		query_wise = "select ra, dec, w1mpro, w1sigmpro, w2mpro, w2sigmpro, " + \
					"w3mpro, w3sigmpro, w4mpro, w4sigmpro FROM " + \
					"allwise.main WHERE ra > " + ra_min + " and ra < " + \
					ra_max + " and dec > " + dec_min + " and dec < " + dec_max

		RA_wise, DEC_wise, W1, W1err, W2, W2err, W3, W3err, W4, W4err = \
		sqlutil.get(query_wise, db="wsdb", host=host, user=uname, password=pword)

		info = [RA_wise, DEC_wise, W1, W1err, W2, W2err, W3, W3err, W4, W4err]

		if len(RA_wise) == 0:
			for (i, col) in enumerate(cols):
				t[col] = [np.float64(np.nan)]*len(t)

		else:
			print "Retreived data, matching lists"
			dists, inds = match_lists.match_lists(t[RA_main], t[DEC_main], RA_wise, DEC_wise, width, 1)
			print "Matched lists, updating table"
			match = np.where( (inds <> len(RA_wise)) )[0]
			no_match = np.where( (inds == len(RA_wise)) )[0]

			for (i, col) in enumerate(cols):
				t[col][match] = info[i][inds[match]]
				t[col][no_match] = np.float64(np.nan)

		print "Table updated"

	elif len(t) <= 5000:
		while n < len(t):
			RA = t[RA_main][n]
			DEC = t[DEC_main][n]
			query = "select ra, dec, w1mpro, w1sigmpro, w2mpro, w2sigmpro, " + \
					"w3mpro, w3sigmpro, w4mpro, w4sigmpro FROM allwise.main " + \
					"WHERE q3c_radial_query(ra, dec, " + str(RA) + ", " + \
					str(DEC) + ", 0.1002)"

			RA_W, DEC_W, W1, W1err, W2, W2err, W3, W3err, W4, W4err = \
			sqlutil.get(query, db='wsdb',host=host, user=uname, password=pword)

			if len(RA_W) > 0:
				dists, inds = match_lists.match_lists([RA], [DEC], RA_W, DEC_W, width, 1)
				j = inds[0]

				if j <> len(RA_W):
					info = [RA_W[j], DEC_W[j], W1[j], W1err[j], W2[j], W2err[j], \
							W3[j], W3err[j], W4[j], W4err[j]]

					for (i, col) in enumerate(cols):
						t[col][n] = info[i]

				else:
					for (i, col) in enumerate(cols):
						t[col][n] = np.float64(np.nan)

			else:
				for (i, col) in enumerate(cols):
					t[col][n] = np.float64(np.nan)

			n += 1

	elif len(t) > 1000000000:

		print "Splitting by tile"
		tiles_all = t["TILENAME"]
		tiles = set(tiles_all)
		tiles = list(tiles)
		print len(tiles), "distinct tiles"

		t.sort("TILENAME")
		tiles.sort()

		n = 0
		m = 0
		rows = []
		tile_rows = []

		while n < len(t):
			tile = tiles[m]
			if t["TILENAME"][n] == tile:
				rows.append(n)
				n += 1
			else:
				tile_rows.append(rows)
				rows = []
				m += 1

		for rows in tile_rows:
			t1 = t.rows(rows)
			print "Matching", t1["TILENAME"][0], "with length", len(t1)
			t1 = WISE_match(t1, RA_main, DEC_main, width, c_graph = False)

			if tile_rows.index(rows) == 0:
				t_out = t1
			else:
				t_out.append(t1)	
		t = t_out			

	if c_graph:
		width = width *3600.0
		RA_med, DEC_med = check_graph(t, RA_main, DEC_main, "WISE", width, gtype = "square")
			
	return t

def VST_match(t, RA_main, DEC_main, c_graph = True, width = 0.0004, w_units = "degrees", DES_box = False):
	
	"""
	Matches to VST Data from WSDB
	Table used is vst_201310.atlas_rgi
	DES_box returns all objects within a box of about the size of a DES tile centred
	on the median RA and DEC of the table.
	Iterates through a table where the RA column is RA_main and the DEC_column is 
	DEC_main
	Puts NaNs in for no matches
	c_graph plots diagnostic graphs of the match - crashes if there are no matches
	"""

	import numpy as np
	import match_lists
	import sqlutil
	import ConfigParser

	config_file = "/home/sr525/Python_Code/srpylib/srpylib.cfg"
	config = ConfigParser.RawConfigParser()
	config.read(config_file)
	uname = config.get("wsdb", "uname")
	pword = config.get("wsdb", "pword")
	host = config.get("wsdb", "host")

	if w_units == "arcsecs":
		width = width / 3600.0
	
	cols = ["RA_VST", "DEC_VST", "G_VST", "G_VST_ERR", "R_VST", "R_VST_ERR", \
			"I_VST", "I_VST_ERR"]

	try:
		for col in cols:
			t.add_empty_column(col, dtype = np.float64)

	except AttributeError as e:
		for col in cols:
			t[col] = np.zeros(len(t))

	if len(t) > 1000 and len(t) < 1000000000:
		ra_min, ra_max, dec_min, dec_max = min_max_match(t[RA_main], t[DEC_main], \
														width, DES_box = DES_box)

		print "Retreiving data"
		query = "select ra, dec, mag_g, magerr_g, mag_r, magerr_r, mag_i, " + \
				"magerr_i FROM vst_201310.atlas_rgi WHERE ra > " + \
				ra_min + " and ra < " + ra_max + " and dec > " + \
				dec_min + " and dec < " + dec_max

		RA_vst, DEC_vst, g, gerr, r, rerr, i, ierr = \
		sqlutil.get(query, db="wsdb", host=host, user=uname, password=pword)

		print "Retreived data"
		dists, inds = match_lists.match_lists(t[RA_main], t[DEC_main], RA_vst, DEC_vst, width, 1)
		print "Matched lists, updating table"

		info = [RA_vst, DEC_vst, g, gerr, r, rerr, i, ierr]

		match = np.where( (inds <> len(RA_vst)) )[0]
		no_match = np.where( (inds == len(RA_vst)) )[0]

		for (i, col) in enumerate(cols):
			t[col][match] = info[i][inds[match]]
			t[col][no_match] = np.float64(np.nan)

		print "Table updated"

	elif len(t) <= 1000:
		n = 0
		while n < len(t):
			RA = t[RA_main][n]
			DEC = t[DEC_main][n]
			query_vhs = "select ra, dec, mag_g, magerr_g, mag_r, magerr_r, " + \
						"mag_i, magerr_i FROM vst_201310.atlas_rgi " + \
						"WHERE q3c_radial_query(ra, dec, " + str(RA) + ", " + \
						str(DEC) + ", 0.1102)"

			RA_vst, DEC_vst, g, gerr, r, rerr, i, ierr = \
			sqlutil.get(query_vhs, db="wsdb", host=host, user=uname, password=pword)

			if len(RA_vst) > 0:
				dists, inds = match_lists.match_lists([RA], [DEC], RA_vst, DEC_vst, width, 1)
				j = inds[0]
				
				if j <> len(RA_vst):
					info = [RA_vst[j], DEC_vst[j], g[j], gerr[j], r[j], rerr[j], \
							i[j], ierr[j]]

					for (k, col) in enumerate(cols):
						t[col][n] = info[k]

				else:
					for (k, col) in enumerate(cols):
						t[col][n] = np.float64(np.nan)

			else:
				for (k, col) in enumerate(cols):
					t[col][n] = np.float64(np.nan)

			n += 1

	elif len(t) > 1000000000:

		print "Splitting by tile"
		tiles_all = t["TILENAME"]
		tiles = set(tiles_all)
		tiles = list(tiles)
		print len(tiles), "distinct tiles"

		t.sort("TILENAME")
		tiles.sort()

		n = 0
		m = 0
		rows = []
		tile_rows = []

		while n < len(t):
			tile = tiles[m]
			if t["TILENAME"][n] == tile:
				rows.append(n)
				n += 1
			else:
				tile_rows.append(rows)
				rows = []
				m += 1

		for rows in tile_rows:
			t1 = t.rows(rows)
			print "Matching", t1["TILENAME"][0], "with length", len(t1)
			t1 = VST_match(t1, RA_main, DEC_main, width, c_graph = False)

			if tile_rows.index(rows) == 0:
				t_out = t1
			else:
				t_out.append(t1)
		t = t_out

	if c_graph:
		width = width *3600.0
		RA_med, DEC_med = check_graph(t, RA_main, DEC_main, "VST", width)

	return t

def GALEX_match(t, RA_main, DEC_main, width = 0.0004, w_units = "degrees", c_graph = True, table = "ais"):

	"""
	Matches to GALEX table
	DES_box returns all objects within a box of about the size of a DES tile centred
	on the median RA and DEC of the table.
	Puts NaNs in for no matches
	c_graph plots diagnostic graphs of the match
	"""

	import numpy as np
	import match_lists
	import sqlutil
	import ConfigParser

	config_file = "/home/sr525/Python_Code/srpylib/srpylib.cfg"
	config = ConfigParser.RawConfigParser()
	config.read(config_file)
	uname = config.get("wsdb", "uname")
	pword = config.get("wsdb", "pword")
	host = config.get("wsdb", "host")

	if w_units == "arcsecs":
		width = width / 3600.0

	cols = ["RA_GALEX", "DEC_GALEX", "NUV_MAG", "NUV_MAGERR", "FUV_MAG", \
			"FUV_MAGERR", "E_BV", "NUV_CLASS_STAR", "FUV_CLASS_STAR"]

	try:
		for col in cols:
			t.add_empty_column(col, dtype = np.float64)

	except AttributeError as e:
		for col in cols:
			t[col] = np.zeros(len(t))

	if len(t) > 1000:

		ra_min, ra_max, dec_min, dec_max = min_max_match(t[RA_main], t[DEC_main], \
														width, DES_box = DES_box)

		print "Retreiving data"

		query = "select ra, dec, nuv_mag, nuv_magerr, fuv_mag, fuv_magerr, e_bv, \
			nuv_class_star, fuv_class_star FROM galexgr6." + table + "primary \
			WHERE ra > " + ra_min + " and ra < " + ra_max + " and dec > " + \
			dec_min + " and dec < " + dec_max

		RA_galex, DEC_galex, nuv_mag, nuv_magerr, fuv_mag, fuv_magerr, e_bv, \
		nuv_class_star, fuv_class_star = \
		sqlutil.get(query, db = "wsdb", host=host, user=uname, password=pword)

		print "Retreived data"

		if len(RA_galex) == 0:
			for (i, col) in enumerate(cols):
				t[col] = [np.float64(np.nan)]*len(t)

		else:
			dists, inds = match_lists.match_lists(t[RA_main], t[DEC_main], RA_vhs, DEC_vhs, width, 1)

			info = [ra_galex, dec_galex, nuv_mag, nuv_magerr, fuv_mag, fuv_magerr, \
					e_bv, nuv_class_star, fuv_class_star]

			match = np.where( (inds <> len(RA_galex)) )[0]
			no_match = np.where( (inds == len(RA_galex)) )[0]

			for (i, col) in enumerate(cols):
				t[col][match] = info[i][inds[match]]
				t[col][no_match] = np.float64(np.nan)

	elif len(t) <= 1000:
		n = 0

		while n < len(t):
			RA = t[RA_main][n]
			DEC = t[DEC_main][n]

			query_galex = "select ra, dec, nuv_mag, nuv_magerr, fuv_mag, \
							fuv_magerr, e_bv, nuv_class_star, fuv_class_star FROM \
							galexgr6." + table + "primary WHERE \
							q3c_radial_query(ra, dec, " + str(RA) + ", " + str(DEC) + ", 0.1102)"

			RA_galex, DEC_galex, nuv_mag, nuv_magerr, fuv_mag, fuv_magerr, e_bv, \
			nuv_class_star, fuv_class_star = sqlutil.get(query_galex, db = "wsdb", \
			host=host, user=uname, password=pword)

			if len(RA_galex) > 0:
				dists, inds = match_lists.match_lists([RA], [DEC], RA_galex, DEC_galex, width, 1)
				j = inds[0]

				if j <> len(RA_galex):
					info = [RA_galex[j], DEC_galex[j], nuv_mag[j], nuv_magerr[j], \
							fuv_mag[j], fuv_magerr[j], e_bv[j], nuv_class_star[j], \
							fuv_class_star[j]]

					for (k, col) in enumerate(cols):
						t[col][n] = info[rki]

				else:
					for (k, col) in enumerate(cols):
						t[col][n] = np.float64(np.nan)

			else:
				for (k, col) in enumerate(cols):
					t[col][n] = np.float64(np.nan)

			n += 1

		print "Done GALEX"

	elif len(t) > 1000000000:

		print "Splitting by tile"
		tiles_all = t["TILENAME"]
		tiles = set(tiles_all)
		tiles = list(tiles)
		print len(tiles), "distinct tiles"

		t.sort("TILENAME")
		tiles.sort()

		n = 0
		m = 0
		rows = []
		tile_rows = []

		while n < len(t):
			tile = tiles[m]
			if t["TILENAME"][n] == tile:
				rows.append(n)
				n += 1
			else:
				tile_rows.append(rows)
				rows = []
				m += 1

		for rows in tile_rows:
			t1 = t.rows(rows)
			print "Matching", t1["TILENAME"][0], "with length", len(t1)
			t1 = GALEX_match(t1, RA_main, DEC_main, width, c_graph = False)

			if tile_rows.index(rows) == 0:
				t_out = t1
			else:
				t_out.append(t1)
		t = t_out

	if c_graph:
		width = width *3600.0
		RA_med, DEC_med = check_graph(t, RA_main, DEC_main, "GALEX", width)

	return t

def VHS_match(t, RA_main, DEC_main, c_graph = True, width = 0.0004, w_units = "degrees", survey = "des", bands = ["Y", "J", "H", "K"], DES_box = False):

	"""
	Matches to VHS Data from WSDB
	Tables used are vhs_1504. + survey
	DES_box returns all objects within a box of about the size of a DES tile centred
	on the median RA and DEC of the table.
	Iterates through a table where the RA column is RA_main and the DEC_column is 
	DEC_main
	Puts NaNs in for no matches
	c_graph plots diagnostic graphs of the match - crashes if there are no matches
	Can specify not the full number of bands
	#29/03/16 - Changed from 1504 to 1603
	"""

	import numpy as np
	import match_lists
	import sqlutil
	from astropy.table import Table, vstack
	import ConfigParser

	config_file = "/home/sr525/Python_Code/srpylib/srpylib.cfg"
	config = ConfigParser.RawConfigParser()
	config.read(config_file)
	uname = config.get("wsdb", "uname")
	pword = config.get("wsdb", "pword")
	host = config.get("wsdb", "host")

	if w_units == "arcsecs":
		width = width / 3600.0

	if survey == "des":
		cols = ["RA_VHS", "DEC_VHS", "J", "J_ERR", "H", "H_ERR", "K", "K_ERR"]

	elif survey == "atlas":
		cols = ["RA_VHS", "DEC_VHS", "Y", "Y_ERR", "J", "J_ERR", "H", "H_ERR", \
				"K", "K_ERR"]

	try:
		for col in cols:
			t.add_empty_column(col, dtype = np.float64)

	except AttributeError as e:
		for col in cols:
			t[col] = np.zeros(len(t))

	ra_min, ra_max, dec_min, dec_max = min_max_match(t[RA_main], t[DEC_main], \
													width, DES_box = DES_box)

	area = (float(ra_max)-float(ra_min))*(float(dec_max)-float(dec_min))
	print area

	if len(t) > 1000 and area < 2000:

		ra_min, ra_max, dec_min, dec_max = min_max_match(t[RA_main], t[DEC_main], \
														width, DES_box = DES_box)

		print "Retreiving data"

		if survey == "atlas":
			query_vhs = "select ra, dec, yapercormag3, yapercormag3_err, " + \
						"japercormag3, japercormag3_err, hapercormag3, " + \
						"hapercormag3_err, kapercormag3, kapercormag3_err " + \
						"FROM vhs_1603." + survey + " WHERE ra > " + ra_min + \
						" and ra < " + ra_max + " and dec > " + dec_min + \
						" and dec < " + dec_max

			RA_vhs, DEC_vhs, Y, Yerr, J, Jerr, H, Herr, K, Kerr  = \
			sqlutil.get(query_vhs, db="wsdb", host=host, user=uname, password=pword)
		
		elif survey == "des":
			query_vhs = "select ra, dec, japercormag3, japercormag3_err, " + \
						"hapercormag3, hapercormag3_err, kapercormag3, " + \
						"kapercormag3_err FROM vhs_1603." + survey + \
						" WHERE ra > " + ra_min + " and ra < " + ra_max + \
						" and dec > " + dec_min + " and dec < " + dec_max

			RA_vhs, DEC_vhs, J, Jerr, H, Herr, K, Kerr  = \
			sqlutil.get(query_vhs, db="wsdb", host=host, user=uname, password=pword)

		print "Retreived data"

		if len(RA_vhs) == 0:
			for (i, col) in enumerate(cols):
				t[col] = [np.float64(np.nan)]*len(t)

		else:
			dists, inds = match_lists.match_lists(t[RA_main], t[DEC_main], RA_vhs, DEC_vhs, width, 1)

			if survey == "des":
				info = [RA_vhs, DEC_vhs, J, Jerr, H, Herr, K, Kerr]

			elif survey == "atlas":
				info = [RA_vhs, DEC_vhs, Y, Yerr, J, Jerr, H, Herr, K, Kerr]

			match = np.where( (inds <> len(RA_vhs)) )[0]
			no_match = np.where( (inds == len(RA_vhs)) )[0]

			for (i, col) in enumerate(cols):
				t[col][match] = info[i][inds[match]]
				t[col][no_match] = np.float64(np.nan)

	if len(t) <= 1000:
		n = 0

		while n < len(t):
			RA = t[RA_main][n]
			DEC = t[DEC_main][n]

			if survey == "atlas":
				query_vhs = "select ra, dec, yapercormag3, yapercormag3_err, " + \
							"japercormag3, japercormag3_err, hapercormag3, " + \
							"hapercormag3_err, kapercormag3, kapercormag3_err " + \
							"FROM vhs_1504." + survey + \
							" WHERE q3c_radial_query(ra, dec, " + str(RA) + \
							", " + str(DEC) + ", 0.1102)"

				RA_vhs, DEC_vhs, Y, Yerr, J, Jerr, H, Herr, K, Kerr = \
				sqlutil.get(query_vhs, db="wsdb", host=host, user=uname, password=pword)

			elif survey == "des":
				query_vhs = "select ra, dec, japercormag3, japercormag3_err, " + \
							"hapercormag3, hapercormag3_err, kapercormag3, " + \
							"kapercormag3_err FROM vhs_1504." + survey + \
							" WHERE q3c_radial_query(ra, dec, " + str(RA) + \
							", " + str(DEC) + ", 0.1102)"

				RA_vhs, DEC_vhs, J, Jerr, H, Herr, K, Kerr = \
				sqlutil.get(query_vhs, db="wsdb", host=host, user=uname, password=pword)

			if len(RA_vhs) > 0:
				dists, inds = match_lists.match_lists([RA], [DEC], RA_vhs, DEC_vhs, width, 1)
				j = inds[0]

				if j <> len(RA_vhs):
					if survey == "des":
						info = [RA_vhs[j], DEC_vhs[j], J[j], Jerr[j], H[j], Herr[j], \
								K[j], Kerr[j]]

					if survey == "atlas":
						info = [RA_vhs[j], DEC_vhs[j], Y[j], Yerr[j], J[j], Jerr[j], \
								H[j], Herr[j], K[j], Kerr[j]]

					for (k, col) in enumerate(cols):
						t[col][n] = info[k]

				else:
					for (k, col) in enumerate(cols):
						t[col][n] = np.float64(np.nan)

			else:
				for (k, col) in enumerate(cols):
					t[col][n] = np.float64(np.nan)

			n += 1

		print "Done VHS"

	if area > 2000 and len(t) > 1000:

		t.sort(RA_main)
		l1 = np.ceil(area/1001.0)
		l = len(t)/l1
		i = 0

		while i < l1:
			if (i+1)*l < len(t):
				t1 = t[i*l:(i+1)*l]
			else:
				t1 = t[i*l:]

			t1 = VHS_match(t1, RA_main, DEC_main, c_graph = c_graph, width = width, \
							w_units = w_units, survey = survey, bands = bands, \
							DES_box = DES_box)

			if i == 0:
				t_out = t1
			else:
				t_out = vstack([t_out, t1])

			i += 1

		t = t_out

	if c_graph:
		width = width *3600.0
		RA_med, DEC_med = check_graph(t, RA_main, DEC_main, "VHS", width)

	return t

def VHS_match_all_cols(t, RA_main, DEC_main, survey = "des"):

	"""
	Matches to VHS Data from WSDB
	Table used is vhs_201306.atlas
	Iterates through a table where the RA column is RA_main and the DEC_column is 
	DEC_main
	Puts NaNs in for no matches
	c_graph plots diagnostic graphs of the match - crashes if there are no matches
	Best to do for small tables only - othersise takes a long time and uses lots of 
	memory
	"""

	import numpy
	import match_lists
	import sqlutil
	import atpy
	import numpy as np
	from astropy.table import Table, hstack
	import ConfigParser

	config_file = "/home/sr525/Python_Code/srpylib/srpylib.cfg"
	config = ConfigParser.RawConfigParser()
	config.read(config_file)
	uname = config.get("wsdb", "uname")
	pword = config.get("wsdb", "pword")
	host = config.get("wsdb", "host")

	n = 0
	any_matches = False
	while n < len(t):
		RA = t[RA_main][n]
		DEC = t[DEC_main][n]

		print n+1, "out of", len(t)

		query = "select * FROM vhs_1504." + survey + " WHERE q3c_radial_query(ra, \
				dec, " + str(RA) + ", " + str(DEC) + ", 0.05)"

		try:
			t1 = atpy.Table("postgres", query=query, database="wsdb", host=host, \
							user=uname, password=pword)
			if not any_matches:
				t2 = t1
				any_matches = True
			else:
				t2.append(t1)

		except Exception as e:
			print e

		n += 1

	if not any_matches:
		print "No matches"
		return t

	t_vhs = Table(t2[:])

	dists, inds = match_lists.match_lists(t[RA_main], t[DEC_main], t_vhs["ra"], t_vhs["dec"], 0.0006)
	ms = np.where( (inds < len(t_vhs)) )[0]

	t_vhs = t_vhs[inds[ms]]
	try:
		t = t.rows([ms])
		for col in t2.columns:
			t.add_column(col, t2[col])
	except AttributeError as e:
		t = t[ms]
		t = hstack([t, t_vhs])

	return t

def UKIDSS_match(t, RA_main, DEC_main, c_graph = True, width = 0.0004, w_units = "degrees", DES_box = False):

	"""
	Matches to VHS Data from WSDB
	Tables used are vhs_1504. + survey
	DES_box returns all objects within a box of about the size of a DES tile centred
	on the median RA and DEC of the table.
	Iterates through a table where the RA column is RA_main and the DEC_column is 
	DEC_main
	Puts NaNs in for no matches
	c_graph plots diagnostic graphs of the match - crashes if there are no matches
	Can specify not the full number of bands
	"""

	import numpy as np
	import match_lists
	import sqlutil
	from astropy.table import Table, vstack
	import ConfigParser

	config_file = "/home/sr525/Python_Code/srpylib/srpylib.cfg"
	config = ConfigParser.RawConfigParser()
	config.read(config_file)
	uname = config.get("wsdb", "uname")
	pword = config.get("wsdb", "pword")
	host = config.get("wsdb", "host")

	if w_units == "arcsecs":
		width = width / 3600.0

	cols = ["RA_UKIDSS", "DEC_UKIDSS", "Y_UKIDSS", "Y_ERR_UKIDSS", "J_UKIDSS", \
			"J_ERR_UKIDSS", "H_UKIDSS", "H_ERR_UKIDSS", "K_UKIDSS", "K_ERR_UKIDSS"]

	try:
		for col in cols:
			t.add_empty_column(col, dtype = np.float64)

	except AttributeError as e:
		for col in cols:
			t[col] = np.zeros(len(t))

	ra_min, ra_max, dec_min, dec_max = min_max_match(t[RA_main], t[DEC_main], \
													width, DES_box = DES_box)

	area = (float(ra_max)-float(ra_min))*(float(dec_max)-float(dec_min))
	print area

	if len(t) > 1000 and area < 2000:

		ra_min, ra_max, dec_min, dec_max = min_max_match(t[RA_main], t[DEC_main], \
														width, DES_box = DES_box)

		print "Retreiving data"

		query_ukidss = "select ra, dec, yAperMag3, yAperMag3Err, j_1AperMag3, \
						j_1AperMag3Err, hAperMag3, hAperMag3Err, kAperMag3, \
						kAperMag3Err FROM ukidssdr10.lassource \
						WHERE ra > " + ra_min + " and ra < " + ra_max + " and \
						dec > " + dec_min + " and dec < " + dec_max
		RA_u, DEC_u, Y, Yerr, J, Jerr, H, Herr, K, Kerr  = \
		sqlutil.get(query_ukidss, db="wsdb", host=host, user=uname, password=pword)

		print "Retreived data"

		if len(RA_u) == 0:
			for (i, col) in enumerate(cols):
				t[col] = [np.float64(np.nan)]*len(t)

		else:
			dists, inds = match_lists.match_lists(t[RA_main], t[DEC_main], RA_u, DEC_u, width, 1)

			info = [RA_u, DEC_u, Y, Yerr, J, Jerr, H, Herr, K, Kerr]

			match = np.where( (inds <> len(RA_u)) )[0]
			no_match = np.where( (inds == len(RA_u)) )[0]

			for (i, col) in enumerate(cols):
				t[col][match] = info[i][inds[match]]
				t[col][no_match] = np.float64(np.nan)

	if len(t) <= 1000:
		n = 0

		while n < len(t):
			RA = t[RA_main][n]
			DEC = t[DEC_main][n]

			query_ukidss = "select ra, dec, yAperMag3, yAperMag3Err, j_1AperMag3, \
							j_1AperMag3Err, hAperMag3, hAperMag3Err, kAperMag3, \
							kAperMag3Err FROM ukidssdr10.lassource WHERE \
							q3c_radial_query(ra, dec, " + str(RA) + ", " + \
							str(DEC) + ", 0.1102)"

			RA_u, DEC_u, Y, Yerr, J, Jerr, H, Herr, K, Kerr = \
			sqlutil.get(query_ukidss, db="wsdb", host=host, user=uname, password=pword)

			if len(RA_u) > 0:
				dists, inds = match_lists.match_lists([RA], [DEC], RA_u, DEC_u, width, 1)
				j = inds[0]

				if j <> len(RA_u):
					info = [RA_u[j], DEC_u[j], Y[j], Yerr[j], J[j], Jerr[j], \
							H[j], Herr[j], K[j], Kerr[j]]

					for (k, col) in enumerate(cols):
						t[col][n] = info[k]

				else:
					for (k, col) in enumerate(cols):
						t[col][n] = np.float64(np.nan)

			else:
				for (k, col) in enumerate(cols):
					t[col][n] = np.float64(np.nan)

			n += 1

		print "Done UKIDSS"

	if area > 2000 and len(t) > 1000:

		t.sort(RA_main)
		l1 = np.ceil(area/1001.0)
		l = len(t)/l1
		i = 0

		while i < l1:
			if (i+1)*l < len(t):
				t1 = t[i*l:(i+1)*l]
			else:
				t1 = t[i*l:]

			t1 = UKIDSS_match(t1, RA_main, DEC_main, c_graph = c_graph, width = width, \
							w_units = w_units, DES_box = DES_box)

			if i == 0:
				t_out = t1
			else:
				t_out = vstack([t_out, t1])

			i += 1

		t = t_out

	if c_graph:
		width = width *3600.0
		RA_med, DEC_med = check_graph(t, RA_main, DEC_main, "UKIDSS", width)

	return t


def DES_match(t, RA_main, DEC_main, c_graph = True, width = 0.0004, w_units = "degrees", release = "Y1A1"):

	"""
	Matches to DES data ingested into WSDB
	Tables used are des_y1a1.coadd_objects and des_y2q1.objects
	Iterates through a table where the RA column is RA_main and the DEC_column is 
	DEC_main
	Puts NaNs in for no matches
	c_graph plots diagnostic graphs of the match - crashes if there are no matches
	"""

	import numpy as np
	import match_lists
	import sqlutil
	from astropy.table import Table, vstack
	import ConfigParser

	config_file = "/home/sr525/Python_Code/srpylib/srpylib.cfg"
	config = ConfigParser.RawConfigParser()
	config.read(config_file)
	uname = config.get("wsdb", "uname")
	pword = config.get("wsdb", "pword")
	host = config.get("wsdb", "host")

	if w_units == "arcsecs":
		width = width / 3600.0

	if release == "Y1A1":
		cols = ["COADD_OBJECTS_ID", "RA_DES_Y1A1", "DEC_DES_Y1A1", "TILENAME_Y1A1", \
				"RUN_Y1A1", "MPSF_G_Y1A1", "MPSF_ERR_G_Y1A1", "MPSF_R_Y1A1", \
				"MPSF_ERR_R_Y1A1", "MPSF_I_Y1A1", "MPSF_ERR_I_Y1A1", "MPSF_Z_Y1A1", \
				"MPSF_ERR_Z_Y1A1", "MPSF_Y_Y1A1", "MPSF_ERR_Y_Y1A1", "MAPER3_G_Y1A1", \
				"MAPER3_ERR_G_Y1A1", "MAPER3_R_Y1A1", "MAPER3_ERR_R_Y1A1", \
				"MAPER3_I_Y1A1", "MAPER3_ERR_I_Y1A1", "MAPER3_Z_Y1A1", \
				"MAPER3_ERR_Z_Y1A1", "MAPER3_Y_Y1A1", "MAPER3_ERR_Y_Y1A1"]
		str_cols = ["TILENAME_Y1A1", "RUN_Y1A1"]

	elif release == "Y2Q1":
		cols = ["QUICK_OBJECT_ID", "RA_DES_Y2Q1", "DEC_DES_Y2Q1", "MPSF_G_Y2Q1", \
				"MPSF_ERR_G_Y2Q1", "MPSF_R_Y2Q1", "MPSF_ERR_R_Y2Q1", "MPSF_I_Y2Q1", \
				"MPSF_ERR_I_Y2Q1", "MPSF_Z_Y2Q1", "MPSF_ERR_Z_Y2Q1", "MPSF_Y_Y2Q1", \
				"MPSF_ERR_Y_Y2Q1", "MAUTO_G_Y2Q1", "MAUTO_ERR_G_Y2Q1", \
				"MAUTO_R_Y2Q1", "MAUTO_ERR_R_Y2Q1", "MAUTO_I_Y2Q1", \
				"MAUTO_ERR_I_Y2Q1", "MAUTO_Z_Y2Q1", "MAUTO_ERR_Z_Y2Q1", \
				"MAUTO_Y_Y2Q1", "MAUTO_ERR_Y_Y2Q1"]
		str_cols = []

	try:
		for col in cols:
			if col not in str_cols:
				t.add_empty_column(col, dtype = np.float64)
			else:
				t.add_empty_col(col, dtype = np.string)

	except AttributeError as e:
		for col in cols:
			if col not in str_cols:
				t[col] = np.zeros(len(t))
			else:
				t[col] = ["AAAAAAAAAAAAAAAAAAAAAAAAAA"]


	ra_min, ra_max, dec_min, dec_max = min_max_match(t[RA_main], t[DEC_main], \
													width, DES_box = False)

	area = (float(ra_max)-float(ra_min))*(float(dec_max)-float(dec_min))
	print area

	if len(t) > 1000 and area < 2000:

		ra_min, ra_max, dec_min, dec_max = min_max_match(t[RA_main], t[DEC_main], \
														width, DES_box = False)

		print "Retreiving data"

		if release == "Y1A1":
			query_des = "select COADD_OBJECTS_ID, RA, DEC, TILENAME, RUN, MAG_PSF_G, \
						MAGERR_PSF_G, MAG_PSF_R, MAGERR_PSF_R, MAG_PSF_I, MAGERR_PSF_I, \
						MAG_PSF_Z, MAGERR_PSF_Z, MAG_PSF_Y, MAGERR_PSF_Y, MAG_APER_3_G, \
						MAGERR_APER_3_G, MAG_APER_3_R, MAGERR_APER_3_R, MAG_APER_3_I, \
						MAGERR_APER_3_I, MAG_APER_3_Z, MAGERR_APER_3_Z, MAG_APER_3_Y, \
						MAGERR_APER_3_Y FROM des_y1a1.coadd_objects WHERE \
						ra > " + ra_min + " and ra < " + ra_max + " and dec > " + \
						dec_min + " and dec < " + dec_max

			co_id, RA_des, DEC_des, tile, run, G, Gerr, R, Rerr, I, Ierr, Z, Zerr, Y, \
			Yerr, Ga, Gaerr, Ra, Raerr, Ia, Iaerr, Za, Zaerr, Ya, Yaerr  = \
			sqlutil.get(query_des, db="wsdb", host=host, user=uname, password=pword)
		
		elif release == "Y2Q1":
			query_des = "select QUICK_OBJECT_ID, RA, DEC, MAG_PSF_G, MAGERR_PSF_G, \
						MAG_PSF_R, MAGERR_PSF_R, MAG_PSF_I, MAGERR_PSF_I, MAG_PSF_Z, \
						MAGERR_PSF_Z, MAG_PSF_Y, MAGERR_PSF_Y, MAG_AUTO_G, \
						MAGERR_AUTO_G, MAG_AUTO_R, MAGERR_AUTO_R, MAG_AUTO_I, \
						MAGERR_AUTO_I, MAG_AUTO_Z, MAGERR_AUTO_Z, MAG_AUTO_Y, \
						MAGERR_AUTO_Y FROM des_y2q1.objects WHERE ra > " + ra_min \
						+ " and ra < " + ra_max + " and dec > " + dec_min + \
						" and dec < " + dec_max

			q_id, RA_des, DEC_des, G, Gerr, R, Rerr, I, Ierr, Z, Zerr, Y, Yerr, Ga, \
			Gaerr, Ra, Raerr, Ia, Iaerr, Za, Zaerr, Ya, Yaerr = \
			sqlutil.get(query_des, db="wsdb", host=host, user=uname, password=pword)

		print "Retreived data"

		if len(RA_des) == 0:
			for (i, col) in enumerate(cols):
				t[col] = [np.float64(np.nan)]*len(t)

		else:
			dists, inds = match_lists.match_lists(t[RA_main], t[DEC_main], RA_des, DEC_des, width, 1)

			if release == "Y1A1":
				info = [co_id, RA_des, DEC_des, tile, run, G, Gerr, R, Rerr, I, \
				Ierr, Z, Zerr, Y, Yerr, Ga, Gaerr, Ra, Raerr, Ia, Iaerr, Za, Zaerr, Ya, Yaerr]

			elif release == "Y2Q1":
				info = [q_id, RA_des, DEC_des, G, Gerr, R, Rerr, I, Ierr, Z, Zerr, \
				Y, Yerr, Ga, Gaerr, Ra, Raerr, Ia, Iaerr, Za, Zaerr, Ya, Yaerr]

			match = np.where( (inds <> len(RA_des)) )[0]
			no_match = np.where( (inds == len(RA_des)) )[0]

			for (i, col) in enumerate(cols):
				t[col][match] = info[i][inds[match]]
				t[col][no_match] = np.float64(np.nan)

	if len(t) <= 1000:
		n = 0

		while n < len(t):
			RA = t[RA_main][n]
			DEC = t[DEC_main][n]

			if release == "Y1A1":
				query_des = "select COADD_OBJECTS_ID, RA, DEC, TILENAME, RUN, MAG_PSF_G, \
							MAGERR_PSF_G, MAG_PSF_R, MAGERR_PSF_R, MAG_PSF_I, MAGERR_PSF_I, \
							MAG_PSF_Z, MAGERR_PSF_Z, MAG_PSF_Y, MAGERR_PSF_Y, MAG_APER_3_G, \
							MAGERR_APER_3_G, MAG_APER_3_R, MAGERR_APER_3_R, MAG_APER_3_I, \
							MAGERR_APER_3_I, MAG_APER_3_Z, MAGERR_APER_3_Z, MAG_APER_3_Y, \
							MAGERR_APER_3_Y FROM des_y1a1.coadd_objects WHERE \
							q3c_radial_query(ra, dec, " + str(RA) + ", " + str(DEC) + ", 0.1102)"

				co_id, RA_des, DEC_des, tile, run, G, Gerr, R, Rerr, I, Ierr, Z, \
				Zerr, Y, Yerr, Ga, Gaerr, Ra, Raerr, Ia, Iaerr, Za, Zaerr, Ya, Yaerr  = \
				sqlutil.get(query_des, db="wsdb", host=host, user=uname, password=pword)

			elif release == "Y2Q1":
				query_des = "select QUICK_OBJECT_ID, RA, DEC, MAG_PSF_G, MAGERR_PSF_G, MAG_PSF_R, \
							MAGERR_PSF_R, MAG_PSF_I, MAGERR_PSF_I, MAG_PSF_Z, \
							MAGERR_PSF_Z, MAG_PSF_Y, MAGERR_PSF_Y,  MAG_AUTO_G, \
							MAGERR_AUTO_G, MAG_AUTO_R, MAGERR_AUTO_R, MAG_AUTO_I, \
							MAGERR_AUTO_I, MAG_AUTO_Z, MAGERR_AUTO_Z, MAG_AUTO_Y, \
							MAGERR_AUTO_Y FROM des_y2q1.objects WHERE \
							q3c_radial_query(ra, dec, " + str(RA) + ", " + str(DEC) + ", 0.1102)"
				
				q_id, RA_des, DEC_des, G, Gerr, R, Rerr, I, Ierr, Z, Zerr, Y, Yerr, \
				Ga, Gaerr, Ra, Raerr, Ia, Iaerr, Za, Zaerr, Ya, Yaerr = \
				sqlutil.get(query_des, db="wsdb", host=host, user=uname, password=pword)

			if len(RA_des) > 0:
				dists, inds = match_lists.match_lists([RA], [DEC], RA_des, DEC_des, width, 1)
				j = inds[0]

				if j <> len(RA_des):
					if release == "Y1A1":
						info = [co_id[j], RA_des[j], DEC_des[j], tile[j], run[j], G[j], \
								Gerr[j], R[j], Rerr[j], I[j], Ierr[j], Z[j], \
								Zerr[j], Y[j], Yerr[j], Ga[j], Gaerr[j], Ra[j], \
								Raerr[j], Ia[j], Iaerr[j], Za[j], Zaerr[j], Ya[j], Yaerr[j]]

					if release == "Y2Q1":
						info = [q_id[j], RA_des[j], DEC_des[j], G[j], Gerr[j], R[j], Rerr[j], \
								I[j], Ierr[j], Z[j], Zerr[j], Y[j], Yerr[j], Ga[j], \
								Gaerr[j], Ra[j], Raerr[j], Ia[j], Iaerr[j], Za[j], \
								Zaerr[j], Ya[j], Yaerr[j]]

					for (k, col) in enumerate(cols):
						t[col][n] = info[k]

				else:
					for (k, col) in enumerate(cols):
						t[col][n] = np.float64(np.nan)

			else:
				for (k, col) in enumerate(cols):
					t[col][n] = np.float64(np.nan)

			n += 1

		print "Done DES"

	if area > 2000 and len(t) > 1000:

		t.sort(RA_main)
		l1 = np.ceil(area/1001.0)
		l = len(t)/l1
		i = 0

		while i < l1:
			if (i+1)*l < len(t):
				t1 = t[i*l:(i+1)*l]
			else:
				t1 = t[i*l:]

			t1 = DES_match(t1, RA_main, DEC_main, c_graph = False, width = width, \
							w_units = w_units, release = release)

			if i == 0:
				t_out = t1
			else:
				t_out = vstack([t_out, t1])

			i += 1

		t = t_out

	if c_graph:
		width = width *3600.0
		if release == "Y1A1":
			survey = "DES_Y1A1"
		elif release == "Y2Q1":
			survey = "DES_Y2Q1"
		RA_med, DEC_med = check_graph(t, RA_main, DEC_main, survey, width)

	return t


