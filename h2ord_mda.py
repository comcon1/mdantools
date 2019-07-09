#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, getopt
from os import environ
from read_ndx import GmxIndex
from struct import *
from numpy import *
from scipy.integrate import trapz

from MDAnalysis import *


'''                                                                          
MDAnTOOLs: H2ORD_MDA
Analogous of the 'gmx h2order' but calculating both 1st and 2nd order
parameters (Aman, 2003, BJ). The anlge is calculated betweeen OZ axis and 
the dipole vector of every water molecule. The reference plane corresponds to the 
mean position of some lipid atom. Note that here you must specify names of
oxygen and hydrogens of water molecule in your system.
'''

def main():
    # starting params
    infile = ''
    tprfile = ''
    oufile = ''
    _STIME = 0
    _ETIME = 1000000
    _DTIME = 300
    _DSLICE = 0.05
    _NSLICE = 50
    _FSLICE = -1.0
    __RESNM = ''
    __REFAT = ''
    _NDXNAME = ''
    __SOLNM = 'SOL'
    __SOLATS = ['OW', 'HW1', 'HW2']
    __SOLN = 3
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hi:o:d:t:s:r:w:n:')
    except getopt.error, msg:
        print msg
        print 'for help use -h'
        sys.exit(2)
    # process options
    for o, a in opts:
        if o == '-h':
            print './h2ord.py -ssd0.tpr -iinput.xtc -oout.xvg -t0:300:10000 '+\
                  '-d0.05:50:-1.0 -rDPC:C20 -wSOL:OW:HW1:HW2 -nindex.ndx'
            sys.exit(0)
        elif o == '-i':
            infile = a
        elif o == '-o':
            oufile = a
        elif o == '-s':
            tprfile = a
        elif o == '-t':
            __ar = a.split(':')
            _STIME = int(__ar[0])
            _DTIME = int(__ar[1])
            _ETIME = int(__ar[2])
        elif o == '-d':
            __ar = a.split(':')
            _DSLICE = float(__ar[0])
            _NSLICE = int(__ar[1])
            _FSLICE = float(__ar[2])
        elif o == '-r':
            __ar = a.split(':')
            __RESNM = __ar[0]
            __REFAT = __ar[1]
        elif o == '-w':
            __ar = a.split(':')
            __SOLNM = __ar[0]
            __SOLATS = __ar[1:]
        elif o == '-n':
	    _NDXNAME = a
           
    if infile == '' or oufile == '' or tprfile == '' \
        or __RESNM == '' or __REFAT == '' or _NDXNAME == '':
        print 'for help use -h'
        sys.exit(2)

    # main objects
    print 'Loading Universe 1. wait..'
    un = Universe(tprfile, infile)
    

    # work with trajectory
    maxframes = int(_ETIME/un.trajectory.dt)
    dframes   = int(_DTIME/un.trajectory.dt)
    skipframes = int(_STIME/un.trajectory.dt)
    
    ndxObj = GmxIndex(_NDXNAME)
    ndx  = ndxObj.getRNdx(__SOLNM) # index array 4 sol. group
    print 'Group ',__SOLNM,', ',len(ndx),' elems'    
    
    # skipping
    for i in range(skipframes):
      un.trajectory.next()


    print '''
      Starting h2o order calculations over %d slices: %8.3f-%-8.3f [nm]
      Trajectory borders: %d - %d [ps] | %d - %d [fr]
      Statistical frame: %d [ps] | %d [fr]
    ''' % ( _NSLICE, _FSLICE, (_FSLICE+_DSLICE*_NSLICE), _STIME, _ETIME, \
        skipframes, maxframes, _DTIME, dframes )


    s = ''

    # initial Z-box define
    Z0_= array(un.trajectory.ts.dimensions[2])
    
    # pointer to useful
    natm = un.atoms.n_atoms
    
    # refatm
    _zh12   = un.dimensions[2]/2.
    refat1 = un.selectAtoms( 'name %s and resname %s and prop z < %f' % \
        (__REFAT, __RESNM, _zh12) )
    refat2 = un.selectAtoms( 'name %s and resname %s and prop z > %f' % \
        (__REFAT, __RESNM, _zh12) )

    # scanning trajectory
    FI_STAT = []
    finished = False
    # useful array
    ztmp = zeros((len(ndx)/__SOLN,3))
    ztmp[:,2] = 1.
    
    # main cycle
    while not finished:
        print '\nOne! - ',un.trajectory.ts.time
        Z0_= array(un.trajectory.ts.dimensions[2])

        frames = 0
        _cum = zeros((4,_NSLICE)) # angles of statistical period
        _frm = zeros((4,_NSLICE)) # density of single frames [ fi | cos(fi) | 1.5cos2(fi)-0.5 ]
        mtrx = zeros((len(ndx)/__SOLN,3,3)) # MOL ATOM AXIS
        
        # over one statistical period
        for curfr in range(dframes):
	    # common trajectory stuff
	    _ts = un.trajectory.ts
	    box = _ts.dimensions[:3]
            frames += 1            
            z0 = refat1.center_of_geometry()[2]
            z1 = refat2.center_of_geometry()[2]
            zc = (z0 + z1)/2.
            xxx = un.trajectory.ts._pos
            _frm *= 0.
            
            # adding anlges
            mtrx = reshape(ravel(xxx[ndx]), (len(ndx)/__SOLN,__SOLN,3))
            mtrx = mtrx[:,:3,:]
            dipoles = (mtrx[:,1,:] + mtrx[:,2,:] - 2*mtrx[:,0,:])*.5 # calculate all water dipole vectors
            # normalize all water dipole vectors and compute dot to z-axis, then sum trace and acos
            dipoles = (multiply(dipoles, array([1./sqrt(sum(dipoles * dipoles, axis=1))]).transpose())*ztmp)[:,2]
            dipoles = arccos(dipoles)
            
            # start summarize over slices
	    _zm = zeros((len(ndx)/__SOLN,2))
	    _zm *= 0.
	    _zm[:,0] = xxx[ndx[::__SOLN],2]
	    _zm[:,1] = dipoles
	    # sort for Z
	    _as = _zm.argsort(axis=0)
	    _zm = _zm[_as[:,0],:]
	    # determine divider
	    _dv = _zm[:,0].searchsorted(zc)
	    # adding both histogramms for both monolayers
	    # FI
	    _mdplus, _edges = histogramdd(_zm[_dv:,0], bins=_NSLICE,\
            range=[(z1+_FSLICE, z1+(_FSLICE+_DSLICE*_NSLICE))], weights=_zm[_dv:,1])
	    _frm[0,:] += _mdplus
	    _mdplus, _edges = histogramdd(_zm[:_dv,0], bins=_NSLICE,\
            range=[(z0-(_FSLICE+_DSLICE*_NSLICE), z0-_FSLICE)], weights=-_zm[:_dv,1]) #minus sign appears because of acos assymmetry
	    _frm[0,:] += _mdplus[::-1] # invert histogramm
	    # COS(FI)
	    _mdplus, _edges = histogramdd(_zm[_dv:,0], bins=_NSLICE,\
            range=[(z1+_FSLICE, z1+(_FSLICE+_DSLICE*_NSLICE))], weights=cos(_zm[_dv:,1]))
	    _frm[1,:] += _mdplus
	    _mdplus, _edges = histogramdd(_zm[:_dv,0], bins=_NSLICE,\
            range=[(z0-(_FSLICE+_DSLICE*_NSLICE), z0-_FSLICE)], weights=cos(_zm[:_dv,1])) # minus is nothing for cosine
	    _frm[1,:] += _mdplus[::-1] # invert histogramm
	    # 1.5 cos^2 (fi) - 0.5
	    _mdplus, _edges = histogramdd(_zm[_dv:,0], bins=_NSLICE,\
            range=[(z1+_FSLICE, z1+(_FSLICE+_DSLICE*_NSLICE))], weights=1.5*cos(_zm[_dv:,1])**2-0.5 )
	    _frm[2,:] += _mdplus
	    _mdplus, _edges = histogramdd(_zm[:_dv,0], bins=_NSLICE,\
            range=[(z0-(_FSLICE+_DSLICE*_NSLICE), z0-_FSLICE)], weights=1.5*cos(_zm[:_dv,1])**2-0.5) # minus is nothing for cosine
	    _frm[2,:] += _mdplus[::-1] # invert histogramm
	    # MOL NUMBER
	    _mdplus, _edges = histogramdd(_zm[_dv:,0], bins=_NSLICE,\
            range=[(z1+_FSLICE, z1+(_FSLICE+_DSLICE*_NSLICE))] )
	    _frm[3,:] += _mdplus
	    _mdplus, _edges = histogramdd(_zm[:_dv,0], bins=_NSLICE,\
            range=[(z0-(_FSLICE+_DSLICE*_NSLICE), z0-_FSLICE)] )
	    _frm[3,:] += _mdplus[::-1] # invert histogramm	    
	    # add to global stats
	    _cum += _frm
	
	    if (not un.trajectory.next() or _ts.frame > maxframes):
	      finished = True
	      break
	    print '.',
	    sys.stdout.flush()
	
	# stat
	#print _cum
	#print multiply(_cum[:3,:], array([1.0/_cum[3,:]]) )
        FI_STAT.append( multiply(_cum[:3,:], array([1.0/_cum[3,:]]) ) )
        
    FI_STAT = array(FI_STAT)
    FI_STAT[:,0,:] *= 180./pi
    print FI_STAT
    print FI_STAT.shape
	    
    
    f = open(oufile, 'w')
    f.write ('''
# Starting h2o order calculations over %d slices: %8.3f-%-8.3f [nm]
# Trajectory borders: %d - %d [ps] | %d - %d [fr]
# Statistical frame: %d [ps] | %d [fr] *** NOT WORKING NOW ***
# layer        <fi>      <cos fi>      <1.5cos2(fi)-0.5>
    ''' % ( _NSLICE, _FSLICE, (_FSLICE+_DSLICE*_NSLICE), _STIME, _ETIME, \
        skipframes, maxframes, _DTIME, dframes ) ) 

    for i in range(FI_STAT.shape[2]):
      x = _FSLICE + _DSLICE*i
      f.write('%8.3f    %10g %10g    %10g %10g    %10g %10g' % \
      ( x, average(FI_STAT[:,0,i]), std(FI_STAT[:,0,i])/sqrt(FI_STAT.shape[0]),\
           average(FI_STAT[:,1,i]), std(FI_STAT[:,1,i])/sqrt(FI_STAT.shape[0]),\
           average(FI_STAT[:,2,i]), std(FI_STAT[:,2,i])/sqrt(FI_STAT.shape[0])) )
      f.write('\n')

    f.close()

    print '\nBye!\n'


if __name__ == '__main__':
    main()
