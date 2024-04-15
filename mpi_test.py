

from petsc4py import PETSc
from mpi4py import MPI


from fenics import *
import numpy as np
import PeriodicBoundary as pb
import Homogenization as hm
import GeometryBoxInABox3D as Geometry 

#comm = PETSc.Comm(MPI.COMM_SELF)

common.cpp.log.set_log_level(1) 

#parameters["mesh_partitioner"] = "SCOTCH";
parameters["mesh_partitioner"] = "ParMETIS"

comm = MPI.comm_world
mpiRank = MPI.rank(comm)
print('MPI PROCESS RANK ', mpiRank)


geom = Geometry.Geometry(2., 1., resolution=40, ds=40., sc=20.)
xi_phi = hm.reference_cell_solver_potential(geom)
#xi_phi = hm.reference_cell_solver_diffusion(geom)

# start with:
# mpirun -np 3 python myScript.py

# TODO test if failure of MPI is caused by periodic boundaries
