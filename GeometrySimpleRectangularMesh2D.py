#!/usr/bin/python

import numpy as np
from fenics import *
from Permittivity import Eps
from GeometryBase import GeometryBase
from NeumannBoundaries import *
from PeriodicBoundary import *

class SimpleRectangularMesh2D(GeometryBase):

    def __init__(self, p1, p2, resolution, eps_l, eps_s, res_y=None, 
                 mesh_diagonal='right', tol=DOLFIN_EPS*100.):
        
        super().__init__()
        
        # tolerance to postions/coordinates
        self.tol = tol
        
        # permittivity of liquid and solid phase
        self.eps_l = eps_l
        self.eps_s = eps_s
        
        
        self.p_min = p1
        self.p_max = p2
        self.res = resolution
        if res_y == None:
            self.res_y = self.res
        else:
            self.res_y = res_y
        self.res_z = 0  #  third dim does not exist

        # full mesh over solid and liquid domains
        self.mesh = RectangleMesh(self.p_min, self.p_max, self.res, self.res_y, diagonal=mesh_diagonal)
        
        # mark domains:
        # * Liquid   : 1
        # * Solid    : 2
        # * Unmarked : 0
        self.materials = MeshFunction("size_t", self.mesh, self.geometric_dimension())
        try:            
            self.materials.set_all(0)
            self.subdomain_liquid.mark(self.materials, 1)
            self.subdomain_solid.mark(self.materials, 2)
        except:
            self.materials.set_all(0)
            print('Marking the subdomains was not possible. Domains are not defined')
        
        self.epsilon = self.compute_permittivity()
              
        # mesh over liquid domain
        self.liquid_mesh = SubMesh(self.mesh, self.materials, 1) 
        
        # BOUNDARIES ##################################################
        # Notaion for boundarys
        # first two letter are given by adjacent domains
        # LL: liquid - liquid
        # LS: liquid - solid
        
        self.ll_boundaries = LLBoundaries2D(self.p_min.x(), self.p_min.y(), self.p_max.x(), self.p_max.y(), self.tol)
        self.ls_boundaries = LSBoundaries2D(self.p_min.x(), self.p_min.y(), self.p_max.x(), self.p_max.y(), self.tol)
        
        # mark boundaries:
        # * LL : 1
        # * LS : 2
        # * SS : 3
        self.boundaries = MeshFunction("size_t", self.liquid_mesh, self.geometric_dimension()-1)
        self.boundaries.set_all(0)
        self.ll_boundaries.mark(self.boundaries, 1)
        self.ls_boundaries.mark(self.boundaries, 2)
        
        # apply periodic boudaries 
        self.mesh_pb = PeriodicBoundary2D(self.mesh)
        self.liquid_mesh_pb = PeriodicBoundary2D(self.liquid_mesh)

         
