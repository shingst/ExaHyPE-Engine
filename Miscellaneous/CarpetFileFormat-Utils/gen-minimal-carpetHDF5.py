#!/usr/bin/python
#
# Minimal working example to create a correct CarpetHDF5 file.
# Script provided by Roland Haas in 2017.
#
# Further comments:
#
# Momentan geben sie vor, von 4 MPI ranks geschrieben worden zu sein. Der
# c=... steht bei Carpet fuer "Component" was also in ExaHyPE der
# Octantzelle entsprechen duerfte. Die Components sollten (muessen aber
# glaube ich nicht) ohne Luecken durchnummeriert sein, aber ein Rank kann
# mehr als eine Component in eine Datei schreiben (was mein Beispiel jetzt
# nicht macht). rl waere refinement level aber man kann das auch immer
# auf 0 lassen, selbst wenn mesh refinement da ist, solange es fuer jede
# Stelle im Raum nur einen HDF5 Datensatz gibt der sie abdeckt (sonst
# sollte da das refinement level hin damit VisIt die Daten des feinsten
# Levels nimmt).
# 
# Das Beispiel ist auch nur 2d, 3d lasse auch Euch selbst machen. Eine
# Warnung: HDF5 sortiert Eintraege in z.B. shape als [z,y,x] aber Cactus
# z.B. in origin als [x,y,z] (also gerade andersrum). Evtl. also etwas
# ausprobieren, bis alles stimmt (h5py stellt es evtl. fuer mich um). Das
# ist im Beispiel meine ich falsch ist aber egal weil sowieso alles
# quadratisch ist.
# 
# Das irorigin attribute (integer Koordinate der Ecke der Component) wird
# offenbar nicht von VisIt verwendet, denn momentan ist's falsch und es
# klappt trotzdem.
#

from math import *
import h5py
import numpy

ranks = 4 # number of MPI ranks (one file per rank)
for fn in range(0,ranks):
    origin=[-200.+fn*(400./ranks),-200]
    shape=[(160/ranks)+1,161]
    dx = 2.5
    x,y = numpy.meshgrid(numpy.arange(origin[0],origin[0]+shape[0]*dx,dx),
                         numpy.arange(origin[1],origin[1]+shape[1]*dx,dx))
    r = numpy.sqrt(x**2+y**2)
    with h5py.File("psi4.file_%d.h5" % fn, "w") as outfh:
        g = outfh.create_group("Parameters and Global Attributes")
        g.attrs.create("nioprocs", ranks, dtype="int32")
        for it,t in enumerate(numpy.arange(0.,200.,42.)):
            print it
            # could write multiple components per rank but the c=... number
            # must be unique during each timestep and the sequence should
            # not have holes
            dset = outfh.create_dataset("GRID::r it=%d tl=0 m=0 rl=0 c=%d" % (it,fn),
                                        shape=r.shape,
                                        dtype=r.dtype,
                                        data=r)
            dset.attrs["origin"] = origin
            dset.attrs["iorigin"] = numpy.array([0,0],dtype='int32')
            dset.attrs.create("level", 0, dtype='int32')
            dset.attrs.create("timestep", it, dtype='int32')
            dset.attrs.create("time", t, dtype='float')
            dset.attrs["delta"] = numpy.array([dx,dx], dtype='float')
            dset.attrs["name"] = numpy.string_("GRID::r\0")
            #plt.pcolormesh(x, y, numpy.real(h_ret));
            #plt.savefig("radius%05d.png" % it);
            #plt.clf();
