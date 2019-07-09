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
MDAnTOOLs: FISTAT_MDA
Analogous of the 'gmx potential' classical utility. The difference is the reference plane: here we attach the plane to the mean position of some lipid atom. The rest is the same: space is sliced with selected distance gap and electron density is calculated in every  slice. Finally, the electron density is integrated twice to get the potential.

WARNING: Pay special attention to the key '-r' which set the atom of the reference plane.          
'''

gmx2mv_pot  = 4.8e-10/1e-8*3e5*4*pi # sgse
gmx2sgs_fie = 4.8e-10/1e-8**2*4*pi  # sgse
gmx2mel_den = 1e+9**3/6.022e23 # mol-e/l

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
    __RESNM = 'DMS'
    __REFAT = 'C17'
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
          print './fi-stat.py -ssd0.pqr -iinput.xtc -oout.xvg -b0:300:10000 '+\
              ' -d0.05:50:-1.0 [-toutdens.xvg] -xindex.ndx:SOL:DPC [-rC20:DPC]'
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
            __REFAT = __ar[0]
            __RESNM = __ar[1]
            
            
    if infile == '' or oufile == '' or tprfile == '':
        print 'for help use -h'
        sys.exit(2)

    un = Universe(tprfile, infile)
    print un.trajectory.dt

    maxframes = int(__maxtime/un.trajectory.dt)
    dframes   = int(__dtime/un.trajectory.dt)
    skipframes = int(__skiptime/un.trajectory.dt)
    KP = 18.07 # e/eps0*1e9

    ndx  = [  ] # index array 4 e.group
    ndxC = [  ] # charges array 4 e.group
    _allC = un.atoms.charges
    # system group
    ndx.append( array(range(un.atoms.n_atoms) ) )
    ndxC.append( _allC )
    # user goups
    if _NDXNAME != '':
      ndxObj = GmxIndex(_NDXNAME)
      for __v in _NDXGRPS:
        ndx.append( ndxObj.getRNdx(__v) )
        ndxC.append(_allC[ndx[-1]])
        print 'Group ',__v,', ',len(ndx[-1]),' elems'


    print '''
      Starting potential calculations over %d slices: %8.3f-%-8.3f [nm]
      Trajectory borders: %d - %d [fr]
      Statistical frame: %d [fr]
      Index: %s
      GROUPS: %s
    ''' % ( _NSLICE, _FSLICE, (_FSLICE+_DSLICE*_NSLICE), skipframes, maxframes, \
        dframes, _NDXNAME, ', '.join(_NDXGRPS))


    POT_STAT = []
    CD_STAT = []

    finished = False
    absolute = skipframes
    
    # pointer to useful
    natm = un.atoms.n_atoms 
    
    # refatm
    refat0 = un.select_atoms( 'name %s and resname %s' % \
        (__REFAT, __RESNM) )
    _zh2 = refat0.center_of_geometry()[2]
    del refat0
    # _zh2   = un.dimensions[2]/2. -- does not suite for non-centered bilayer
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
        # initial Z-box define

        frames = 0
        _cum = zeros((2,len(ndx),_NSLICE)) # density of statistical period
        _frm = zeros((2,len(ndx),_NSLICE)) # density of single frame
        _edges = None
        for curfr in range(dframes):
            _ts = un.trajectory.ts
            box = _ts.dimensions[:3]
            frames += 1
            z0 = refat1.center_of_geometry()[2]
            z1 = refat2.center_of_geometry()[2]
            zc = (z0 + z1)/2.
            xxx = un.trajectory.ts._pos
            _frm *= 0.
            # adding e-density
            for g in range(len(ndx)):
              _zm = zeros((len(ndx[g]),2))
              _zm *= 0.
              _zm[:,0] = xxx[ndx[g],2]
              _zm[:,1] = ndxC[g]
              # sort for Z
              _as = _zm.argsort(axis=0)
              _zm = _zm[_as[:,0],:]
              # determine divider
              _dv = _zm[:,0].searchsorted(zc)
              # adding both histogramms for both monolayers
              _mdplus, _edges = histogramdd(_zm[_dv:,0], bins=_NSLICE,\
                  range=[(z1+_FSLICE, z1+(_FSLICE+_DSLICE*_NSLICE))], weights=_zm[_dv:,1])
              _frm[0,g,:] += _mdplus

              _mdplus, _edges = histogramdd(_zm[:_dv,0], bins=_NSLICE,\
                  range=[(z0-(_FSLICE+_DSLICE*_NSLICE), z0-_FSLICE)], weights=_zm[:_dv,1])
              _frm[1,g,:] += _mdplus[::-1] # invert histogramm


            _cum += _frm * (1./box[0]/box[1]/_DSLICE)

            if (not un.trajectory.next() or _ts.frame > maxframes):
              finished = True
              break
            print '.',
            sys.stdout.flush()

        _cum *= 1.0/frames
        DENS = _cum

        edge = array(_edges[0])
        zpos = 0.5*(edge[:-1]+edge[1:])         
        # Arbitrary zero reference point
        zpos -= min(zpos)
        _xzp = array((zpos,)*len(ndx))
	print 'xzp ',_xzp.shape
        x_zpos = zeros((2,_xzp.shape[0], _xzp.shape[1]))
	x_zpos[0,:,:] = _xzp
	x_zpos[1,:,:] = _xzp
	print 'xzp ',x_zpos.shape
	print 'cum ',_cum.shape

 
        # Now integrate charge distribution twice (see derivation in Sachs et al paper)
        sum_q_binwidth = add.accumulate(_cum*_DSLICE, axis=2)
        sum_z_q_binwidth = add.accumulate(x_zpos*_cum*_DSLICE, axis=2)
 
        phi = x_zpos*sum_q_binwidth - sum_z_q_binwidth

        # Enforce periodic boundary and convert units to mV
        # pbc_factor = sum_z_q_binwidth[-1]/zbox

        # converting charge distribution value to mV units
        CD_STAT.append(DENS*gmx2mel_den)
        POT_STAT.append(phi*gmx2mv_pot*-1.0)


    # converting array to numpy
    POT_STAT2 = array(POT_STAT)

    POT_STAT = POT_STAT2[0,:,:,:]
    f = open(oufile+'ml1', 'w')
    f.write ('''
# Starting potential calculations over %d slices: %8.3f-%-8.3f [nm]
# Trajectory borders: %d - %d [fr]
# Statistical frame: %d [fr]
# GROUPS: %s
#   x   |    F,V   |  +-m(F)
''' % ( _NSLICE, _FSLICE, (_FSLICE+_DSLICE*_NSLICE), skipframes, maxframes,\
      dframes, ' '.join(_NDXGRPS) ) )

    for i in range(POT_STAT.shape[2]):
      x = _FSLICE + i*_DSLICE
      f.write('%8.3f  ' % x)
      for g in range(len(ndx)):
        f.write('%10g %10g ' % ( average(POT_STAT[:,g,i]), \
            std(POT_STAT[:,g,i])/sqrt(POT_STAT.shape[0]) ) )
      f.write('\n')
    f.close()

    POT_STAT = POT_STAT2[1,:,:,:]
    f = open(oufile+'ml2', 'w')
    f.write ('''
# Starting potential calculations over %d slices: %8.3f-%-8.3f [nm]
# Trajectory borders: %d - %d [fr]
# Statistical frame: %d [fr]
# GROUPS: %s
#   x   |    F,V   |  +-m(F)
''' % ( _NSLICE, _FSLICE, (_FSLICE+_DSLICE*_NSLICE), skipframes, maxframes,\
      dframes, ' '.join(_NDXGRPS) ) )

    for i in range(POT_STAT.shape[2]):
      x = _FSLICE + i*_DSLICE
      f.write('%8.3f  ' % x)
      for g in range(len(ndx)):
        f.write('%10g %10g ' % ( average(POT_STAT[:,g,i]), \
            std(POT_STAT[:,g,i])/sqrt(POT_STAT.shape[0]) ) )
      f.write('\n')
    f.close()

    if (oufileD != ''):
      # converting array to numpy
      CD_STAT2 = array(CD_STAT)

      CD_STAT = CD_STAT2[0,:,:,:]
      f = open(oufileD+'ml1', 'w')
      f.write ('''
# Starting potential calculations over %d slices: %8.3f-%-8.3f [nm]
# Trajectory borders: %d - %d [fr]
# Statistical frame: %d [fr]
# GROUPS: %s
#   x   |    F,V   |  +-m(F)
  ''' % ( _NSLICE, _FSLICE, (_FSLICE+_DSLICE*_NSLICE), skipframes, maxframes,\
        dframes, ' '.join(_NDXGRPS) ) )

      for i in range(CD_STAT.shape[2]):
        x = _FSLICE + i*_DSLICE
        f.write('%8.3f  ' % x)
        for g in range(len(ndx)):
          f.write('%10g %10g ' % ( average(CD_STAT[:,g,i]), \
              std(CD_STAT[:,g,i])/sqrt(CD_STAT.shape[0]) ) )
        f.write('\n')
      f.close()

      CD_STAT = CD_STAT2[1,:,:,:]
      f = open(oufileD+'ml2', 'w')
      f.write ('''
# Starting potential calculations over %d slices: %8.3f-%-8.3f [nm]
# Trajectory borders: %d - %d [fr]
# Statistical frame: %d [fr]
# GROUPS: %s
#   x   |    F,V   |  +-m(F)
  ''' % ( _NSLICE, _FSLICE, (_FSLICE+_DSLICE*_NSLICE), skipframes, maxframes,\
        dframes, ' '.join(_NDXGRPS) ) )

      for i in range(CD_STAT.shape[2]):
        x = _FSLICE + i*_DSLICE
        f.write('%8.3f  ' % x)
        for g in range(len(ndx)):
          f.write('%10g %10g ' % ( average(CD_STAT[:,g,i]), \
              std(CD_STAT[:,g,i])/sqrt(CD_STAT.shape[0]) ) )
        f.write('\n')
      f.close()



    print '\nBye!\n'


if __name__ == '__main__':
    main()
