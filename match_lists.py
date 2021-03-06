#returns two arrays one with the distances and inf if not matched and one with the indices of the match with list
#length + 1 if not matched
#Doesn't remove duplicates if they are both within the radius

import scipy.spatial.kdtree, numpy
from scipy import version
from numpy import sin,cos,deg2rad,rad2deg,arcsin
import time

scipy_version  = ('.'.join(version.version.split('.')[0:2])).split('.')[0:2]

def match_lists(ra1, dec1, ra2, dec2, dist, numNei=1):
	"""crossmatches the list of objects (ra1,dec1) with
	another list of objects (ra2,dec2) with the matching radius "dist"
	The routines searches for up to numNei closest neighbors
	the routine returns the distance to the neighbor and the list 
	of indices of the neighbor. Everything is in degrees.
	if no match is found the distance is NaN.
	Example: 
	> dist, ind = match_lists(ra1,dec1,ra2,dec2, 1./3600)
	> goodmatch_ind = numpy.isfinite(dist)
	> plot(ra1[goodmatch_ind],ra2[ind[goodmatch_ind]])
	Another example:
	> print match_lists( [1,1], [2,3], [1.1,1.2,6,], [2.1,2.2,10], 0.3,numNei=2)
	    (array([[ 0.1413761 ,  0.28274768],
	            [        inf,         inf]]),
         array([[0, 1],
	            [3, 3]]))
	"""
	cosd = lambda x : cos(deg2rad(x))
	sind = lambda x : sin(deg2rad(x))
	mindist = 2 * sind(dist/2.)	
	getxyz = lambda r, d: [cosd(r)*cosd(d), sind(r)*cosd(d), sind(d)]
	xyz1 = numpy.array(getxyz(ra1, dec1))
	xyz2 = numpy.array(getxyz(ra2, dec2))

	if (int(scipy_version[0])==0) and (int(scipy_version[1])<8):
	# If old scipy version is detected then we use KDTree instead of 
	# cKDTtree because there is a bug in the cKDTree
	# http://projects.scipy.org/scipy/ticket/1178
		tree2 = scipy.spatial.KDTree(xyz2.T)
	else:
		tree2 = scipy.spatial.cKDTree(xyz2.T)
	del xyz2
	ret = tree2.query(xyz1.T, numNei, 0, 2, mindist)
	del xyz1
	dist, ind = ret
	finite = numpy.isfinite(dist)
	dist[finite] = rad2deg(2*arcsin(dist[finite]/2))

	return dist, ind

def match_ids(ids, id1s, numNei = 1):
	ids = numpy.array(ids)
	id1s = numpy.array(id1s)
	if (int(scipy_version[0])==0) and (int(scipy_version[1])<8):
		id_tree = scipy.spatial.KDTree((id1s.T).reshape((len(id1s), 1)))
	else:
		id_tree = scipy.spatial.cKDTree((id1s.T).reshape((len(id1s), 1)))
	dist, ind = id_tree.query((ids.T).reshape((len(ids),1)), k = numNei, eps=0, p=2, distance_upper_bound = 0.1)
	return dist, ind
