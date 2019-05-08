#! /usr/bin/env python

# Script for converting a file of rows/columns, with rows 2,3,4 storing x/y/z of electrode positions
# on a subject's fsnative brain, into fsaverage coordinates.

import os, sys, pandas, pimms, neuropythy as ny, numpy as np, scipy as sp, scipy.spatial as space

def note(s, *argv):
    if len(argv) > 0: s = s % tuple(argv)
    sys.stdout.write(s + '\n')
    sys.stdout.flush()
    return s
def die(s, *argv):
    if len(argv) > 0: s = s % tuple(argv)
    sys.stderr.write(s + '\n')
    sys.stderr.flush()
    sys.exit(1)

if sys.argv != 4:
    dir('SYNTAX: fsnative_to_fsaverage.py <freesurfer_subject_id_or_path> <coordinates_file> <out>')
    
# load the freesurfer subject
sid = sys.argv[1]
try: sub = ny.freesurfer_subject(sid)
except Exception: die('Could not load freesurfer subject: %s', sid)

# load the coords file
coords = pandas.read_csv(sys.argv[2], sep=' ', header=None)
coords = coords[[1,2,3]].values

# find nearest for each of these points on the subject's pial surface
lpial = sub.lh.pial_surface
rpial = sub.rh.pial_surface
(nl,nr) = (lpial.vertex_count, rpial.vertex_count)
n = nl + nr
pial = np.vstack([lpial.coordinates.T, rpial.coordinates.T])
try: sh = space.cKDTree(pial)
except Exception: sh = space.KDTree(pial)

(d,ii) = sh.query(coords)
lpts = (ii < nl)
rpts = ~lpts
(lpts,rpts) = [np.where(pts)[0] for pts in (lpts,rpts)]

fsa = ny.freesurfer_subject('fsaverage')
res = np.full(coords.shape, np.nan)
for (h,pts) in zip(['lh','rh'], [lpts,rpts]):
    if len(pts) == 0: continue
    idcs  = ii[pts]
    if h == 'rh': idcs = idcs - nl
    fshem = fsa.hemis[h]
    sbhem = sub.hemis[h]
    fsmsh = fshem.registrations['fsaverage']
    addrs = fsmsh.address(sbhem.registrations['fsaverage'].coordinates[:,idcs].T)
    res[pts,:] = fshem.pial_surface.unaddress(addrs).T

df = ny.to_dataframe({k:v for (k,v) in zip(['x','y','z'], res.T)})
ny.save(sys.argv[3], df)

note('Success!')
sys.exit(0)

