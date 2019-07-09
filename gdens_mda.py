#!/usr/bin/python
# coding: utf-8

'''              
MDAnTOOLs: GDENS_MDA
Analogous of the 'gmx density' classical utility. The difference is the reference plane: here we attach the plane to the mean position of some lipid atom. The rest is the same: space is sliced with selected distance gap and mass density is calculated in every  slice.

WARNING: Pay special attention to the key '-r' which set the atom of the reference plane.
'''

import sys, getopt
from os import environ
from read_ndx import GmxIndex
from struct import *
from numpy import *

from MDAnalysis import *

aa2sgs = 1.67262158

def main():
    # starting params
    infile = ''
    tprfile = ''
    oufile = ''
    oufileD = ''
    __skiptime = 0
    __maxtime  = 10000
    __dtime  = 1000
    _DSLICE = 0.05
    _NSLICE = 50
    _FSLICE = -1.0
    __RESNM = 'DPC'
    __REFAT = 'C20'
    _NDXNAME = ''
    _NDXGRPS = []
    # parse command line options
    try:
      opts, args = getopt.getopt(sys.argv[1:], 'hi:o:d:b:s:x:t:r:')
    except getopt.error, msg:
        print msg
        print 'for help use -h'
        sys.exit(2)
    # process options
    for o, a in opts:
        if o == '-h':
          print './gdens.py -ssd0.pqr -iinput.xtc -oout.xvg -b0:300:10000 '+\
              ' -d0.05:50:-1.0 -xindex.ndx:SOL:DPC [ -rDPC:C20 -h ]'
          sys.exit(0)
        elif o == '-i':
            infile = a
        elif o == '-o':
            oufile = a
        elif o == '-t':
            oufileD = a
        elif o == '-s':
            tprfile = a
        elif o == '-b':
            __ar = a.split(':')
            __skiptime = int(__ar[0])
            __dtime = int(__ar[1])
            __maxtime = int(__ar[2])
        elif o == '-d':
            __ar = a.split(':')
            _DSLICE = float(__ar[0])
            _NSLICE = int(__ar[1])
            _FSLICE = float(__ar[2])
        elif o == '-x':
            __ar = a.split(':')
            _NDXNAME = __ar[0]
            _NDXGRPS = __ar[1:]
        elif o == '-r':
            __ar = a.split(':')
            __RESNM = __ar[0]
            __REFAT = __ar[1]

            
            
    if infile == '' or oufile == '' or tprfile == '':
        print 'for help use -h'
        sys.exit(2)
    
    un = Universe(tprfile, infile)
    
    maxframes = int(__maxtime/un.trajectory.dt)
    dframes   = int(__dtime/un.trajectory.dt)
    skipframes = int(__skiptime/un.trajectory.dt)

    ndx  = [  ] # index array 4 e. group
    ndxM = [  ] # masses array 4 e. group
    _allM = un.atoms.masses

    # turn off SYSTEM default group
    if _NDXNAME != '':
      ndxObj = GmxIndex(_NDXNAME)
      for __v in _NDXGRPS:
        ndx.append( ndxObj.getRNdx(__v) )
        ndxM.append(_allM[ndx[-1]])
        print 'Group ',__v,', ',len(ndx[-1]),' elems'


    print '''
      Starting density calculations over %d slices: %8.3f-%-8.3f [Ang]
      Trajectory borders: %d - %d [fr]
      Statistical frame: %d [fr]
      Index: %s
      GROUPS: %s
    ''' % ( _NSLICE, _FSLICE, (_FSLICE+_DSLICE*_NSLICE), skipframes, maxframes, \
        dframes, _NDXNAME, ', '.join(_NDXGRPS))


    DENS_STAT = []

    finished = False
    absolute = skipframes
    
    # pointer to useful
    natm = un.atoms.n_atoms
    
    # refatm
    _zh2   = un.dimensions[2]/2.
    refat1 = un.select_atoms( 'name %s and resname %s and prop z < %f' % \
        (__REFAT, __RESNM, _zh2) )
    refat2 = un.select_atoms( 'name %s and resname %s and prop z > %f' % \
        (__REFAT, __RESNM, _zh2) )

    rngs = range(skipframes, maxframes-skipframes, dframes)

    # skipping
    for i in range(skipframes):
      un.trajectory.next()

    while not finished:
        print '\nOne! - ',un.trajectory.ts.time
        Z0_= array(un.trajectory.ts.dimensions[2])

        frames = 0
        _cum = zeros((len(ndx),_NSLICE)) # density of statistical period
        _frm = zeros((len(ndx),_NSLICE)) # density of single frame
        for curfr in range(dframes):
            _ts = un.trajectory.ts
            box = _ts.dimensions[:3]
            frames += 1
            absolute += 1
            z0 = refat1.center_of_geometry()[2]
            z1 = refat2.center_of_geometry()[2]
            zc = (z0 + z1)/2.
            xxx = un.trajectory.ts._pos
            _frm *= 0.
            # adding density
            for g in range(len(ndx)):
              _zm = zeros((len(ndx[g]),2))
              _zm *= 0.
              _zm[:,0] = xxx[ndx[g],2]
              _zm[:,1] = ndxM[g]
              # sort for Z
              _as = _zm.argsort(axis=0)
              _zm = _zm[_as[:,0],:]
              # determine divider
              _dv = _zm[:,0].searchsorted(zc)
              # adding both histogramms for both monolayers
              if (_dv < _zm.shape[0]):
                _mdplus, _edges = histogramdd(_zm[_dv:,0], bins=_NSLICE,\
                    range=[(z1+_FSLICE, z1+(_FSLICE+_DSLICE*_NSLICE))], weights=_zm[_dv:,1])
                _frm[g,:] += _mdplus
              if (_dv > 0):
                _mdplus, _edges = histogramdd(_zm[:_dv,0], bins=_NSLICE,\
                    range=[(z0-(_FSLICE+_DSLICE*_NSLICE), z0-_FSLICE)], weights=_zm[:_dv,1])
                _frm[g,:] += _mdplus[::-1] # invert histogramm

            _cum += _frm * (0.5/box[0]/box[1]/_DSLICE)
	    try:
		    if (not un.trajectory.next() or _ts.frame > maxframes):
		      finished = True
		      break
	    except IOError:
		    finished = True
                    break
            print '.',
            sys.stdout.flush()
        # a.u/ang^3
        DENS_STAT.append( _cum*(1.0/frames)*aa2sgs )
    # end working with trajectory 

    # converting array to numpy
    DENS_STAT = array(DENS_STAT)

    f = open(oufile, 'w')
    f.write ('''
# Starting density calculations over %d slices: %8.3f-%-8.3f [nm]
# Trajectory borders: %d - %d [ps]
# Statistical frame: %d [ps]
# GROUPS: %s
#   x   |    F,V   |  +-m(F)
''' % ( _NSLICE, _FSLICE, (_FSLICE+_DSLICE*_NSLICE), skipframes, maxframes,\
      dframes, ' '.join(_NDXGRPS) ) )

    for i in range(DENS_STAT.shape[2]):
      x = _FSLICE + _DSLICE*i
      f.write('%8.3f  ' % x)
      for g in range(len(ndx)):
        f.write('%10g %10g ' % ( average(DENS_STAT[:,g,i]), \
            std(DENS_STAT[:,g,i])/sqrt(DENS_STAT.shape[0]) ) )
      f.write('\n')
    f.close()


    print '\nBye!\n'


if __name__ == '__main__':
    main()
