#!/usr/bin/python

'''
The example script reads atom positions from XVG files,
replicates periodic conditions and connects the vertexes 
with Delaunay algorithm.
'''

from numpy import *
from matplotlib.pyplot import *
from matplotlib.mlab import *

bo = loadtxt('box.xvg')
i0 = bo[:,0].searchsorted(990000)
i1 = bo[:,0].searchsorted(995000)
box = bo[i0:i1,1:3].mean(axis=0)*10.
print box

def spreadar(a, dim):
    parr = array(a)
    parr2 = array(parr)
    parr2[:,0] += dim[0]
    parr = vstack((parr,parr2))
    parr2 = array(parr)
    parr2[:,1] += dim[1]
    parr = vstack((parr,parr2))
    return parr


oe = loadtxt('netOPE.xvg')
os = loadtxt('netOPS.xvg')
#cl = loadtxt('netCHL.xvg')

figure(figsize=(7,7))
plot(oe[:,0], oe[:,1], 'sr')
plot(os[:,0], os[:,1], 'ob')
#plot(cl[:,0], cl[:,1], 'xk')

oee = spreadar(oe[:,0:2], box)
oss = spreadar(os[:,0:2], box)
#ocl = spreadar(cl[:,0:2], box)

figure(figsize=(7,7))
plot(oee[:,0], oee[:,1], 'sr')
plot(oss[:,0], oss[:,1], 'ob')
#plot(ocl[:,0], ocl[:,1], 'xk')

parr = vstack((oee,oss))

from scipy.spatial import Delaunay

dl = Delaunay(parr)
convex = dl.convex_hull

edge_points = []
edges = set()

def add_edge(i, j):
    """Add a line between the i-th and j-th points, if not in the list already"""
    if (i, j) in edges or (j, i) in edges:
        # already added
        return
    edges.add( (i, j) )
    edge_points.append(parr[ [i, j] ])

# loop over triangles: 
# ia, ib, ic = indices of corner points of the triangle
for ia, ib, ic in dl.vertices:
    add_edge(ia, ib)
    add_edge(ib, ic)
    add_edge(ic, ia)

from matplotlib.collections import LineCollection

lines = LineCollection(edge_points)
gca().add_collection(lines)

xlim(10,110)
ylim(10,110)

show()
