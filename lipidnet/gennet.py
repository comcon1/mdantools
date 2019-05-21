#!/usr/bin/python

import sys, getopt
sys.path.append('/home/comcon1/DEVEL/pylib')
from os import environ
from read_ndx import GmxIndex
from struct import *
from numpy import *

from MDAnalysis import *

'''                                                                          
MDAnTOOLs: GENNET
The script calculates average positions of every hydrocarbon chain of every lipid and averages its position over the trajectory. 

The calculation is performed only for residues indicated as a 'LIPID RESIDUE'. Reference atoms are using only for the determination of which monolayer every lipid belongs to. 

The tail index is the atom number in reference residue listed in groups named RESNAME_tail1 and RESNAME_tail2 (see tails.ndx for example).
'''

def main():
    # starting params
    infile = ''
    tprfile = ''
    oufile = ''
    oufileD = ''
    __skiptime = 0
    __maxtime  = 10000
    __dtime  = 1000
    __RESNM = 'HL1'
    __REFAT = 'O3'
    _NDXNAME = ''
    # parse command line options
    try:
      opts, args = getopt.getopt(sys.argv[1:], 'hi:o:b:s:x:t:r:')
    except getopt.error, msg:
        print msg
        print 'for help use -h'
        sys.exit(2)
    # process options
    for o, a in opts:
        if o == '-h':
          print './gdens.py -ssd0.pqr -iinput.xtc -oout.xvg -b0:300:10000 '+\
              ' -xindex.ndx'
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
        elif o == '-x':
            _NDXNAME = a
        elif o == '-r': # LIPID RESIDUE : REFERENCE ATOM
            __ar = a.split(':')
            __RESNM = __ar[0]
            __REFAT = __ar[1]
            
            
    if infile == '' or oufile == '' or tprfile == '' or _NDXNAME == '':
        print 'for help use -h'
        sys.exit(2)
        
    un = Universe(tprfile, infile)
    
    maxframes = int(__maxtime/un.trajectory.delta)
    dframes   = int(__dtime/un.trajectory.delta)
    skipframes = int(__skiptime/un.trajectory.delta)

    ndx  = [  ] # index array 4 e. group
    _allM = un.atoms.masses()
    
    print '''
      Starting OXY net calculations.
      Trajectory borders: %d - %d [fr]
      Statistical frame: %d [fr]
      Index: %s
    ''' % ( skipframes, maxframes, dframes, _NDXNAME )


    DENS_STAT = []

    finished = False
    absolute = skipframes
    
    # pointer to useful
    natm = un.atoms.numberOfAtoms() 

    rngs = range(skipframes, maxframes-skipframes, dframes)

    # skipping
    for i in range(skipframes):
      un.trajectory.next()


    refat0 = un.selectAtoms( 'name %s and resname %s' % (__REFAT, __RESNM) )
    blcenter = refat0.centerOfGeometry()[2]
    
    refat1 = un.selectAtoms( 'name %s and resname %s and prop z > %.1f' % \
            (__REFAT, __RESNM, blcenter) )
    refat2 = un.selectAtoms( 'name %s and resname %s and prop z < %.1f' % \
            (__REFAT, __RESNM, blcenter) )
    
    ndxObj = GmxIndex(_NDXNAME)
    print 'Monolayer 1: ',refat1.numberOfAtoms(),' lipids'
    print 'Monolayer 2: ',refat2.numberOfAtoms(),' lipids'

    print 'Lipid has 2 hydrocarbon tails:'
    ndx.append( ndxObj.getNdx(__RESNM+'_tail1') )
    print '==Tail 1: ',len(ndx[-1]),' atoms'
    ndx.append( ndxObj.getNdx(__RESNM+'_tail2') )
    print '==Tail 2: ',len(ndx[-1]),' atoms'
    
    # generate indexes for every tail1
    print 'Generating indexes for every tail..',
    TL_1 = []
    for at in refat1.atoms:
      vrtx = (at.residue.indices())[0]-1
      TL_1.append(map(lambda x: x+vrtx, ndx[0]) )
      TL_1.append(map(lambda x: x+vrtx, ndx[1]) )
    
    TL_2 = []
    for at in refat2.atoms:
      vrtx = (at.residue.indices())[0]-1
      TL_2.append(map(lambda x: x+vrtx, ndx[0]) )
      TL_2.append(map(lambda x: x+vrtx, ndx[1]) )
    print 'done.'


    frames = 0
    _ml1 = zeros((len(TL_1),2)) # coordinates of hc tails of 1st ML
    _ml2 = zeros((len(TL_2),2)) # coordinates of hc tails of 2st ML
    
    while not finished:
        print '\nOne! - ',un.trajectory.ts.time
        Z0_= array(un.trajectory.ts.dimensions[2])
        
	_ts = un.trajectory.ts
	box = _ts.dimensions[:3]
        frames += 1
        
        xxx = un.trajectory.ts._pos
	# finding com of every tail
	for i in range(len(TL_1)):
	  TLM = _allM[TL_1[i]]
	  _ml1[i,:] += sum( xxx[TL_1[i],0:2] * transpose(vstack((TLM,TLM))), axis=0 ) / sum(TLM)
	
	for i in range(len(TL_2)):
	  TLM = _allM[TL_2[i]]
	  _ml2[i,:] += sum( xxx[TL_2[i],0:2] * transpose(vstack((TLM,TLM))), axis=0 ) / sum(TLM)

#	print _ml1
#	sys.exit(0)

	if (not un.trajectory.next() or _ts.frame > maxframes):
	  finished = True
	  break
    # end working with trajectory 
    
    print '\n%d frames were processed!' % (frames)
    _ml1 *= 1./frames
    _ml2 *= 1./frames
    
    shp = min(_ml1.shape[0], _ml2.shape[0])
    savetxt(oufile, hstack((_ml1[:shp,:], _ml2[:shp,:])) ) 
    print 'Data successfully saved to ', oufile, '!'
    print 'Bye!'
    
if __name__ == '__main__':
    main()
