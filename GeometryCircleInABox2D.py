#!/usr/bin/python

from fenics import *
from GeometrySimpleRectangularMesh2D import SimpleRectangularMesh2D
from PeriodicBoundary import *
import numpy as np

# DOMAINS ########################################################
        
class Liquid(SubDomain):

    def __init__(self, tol, ds, r, **kwargs):
        self.tol = tol
        self.ds = ds
        self.r = r
        super().__init__(**kwargs)
  
    def inside(self, x, on_boundary):
        
        if ((x[1]-0.5*self.ds) * (x[1]-0.5*self.ds) + (x[0]-0.5*self.ds) * (x[0]-0.5*self.ds))  >= (self.r * self.r - self.tol):
            return True
        else: 
            return True
 
class Solid(SubDomain):

    def __init__(self, tol, ds, r, **kwargs):
        self.tol = tol
        self.ds = ds
        self.r = r
        super().__init__(**kwargs)
        
    def inside(self, x, on_boundary):

        if ((x[1]-0.5*self.ds) * (x[1]-0.5*self.ds) + (x[0]-0.5*self.ds) * (x[0]-0.5*self.ds))  <= (self.r * self.r + self.tol):
            return True
        else: 
            return False                
            


##############################################################

class Geometry(SimpleRectangularMesh2D):

    def __init__(self, eps_l, eps_s,  resolution=200, ds=18., sc=4.51,**kwargs):
        
        self.domain_size = ds
        self.solid_scale = sc
        
        p1 = Point(0., 0.)
        p2 = Point(self.domain_size, self.domain_size)
        
        self.tol=DOLFIN_EPS * 100.
        
        self.subdomain_liquid = Liquid(self.tol, ds=self.domain_size, r=self.solid_scale)
        self.subdomain_solid = Solid(self.tol, ds=self.domain_size, r=self.solid_scale)
        
        SimpleRectangularMesh2D.__init__(self, p1, p2, resolution, eps_l, eps_s, 
        mesh_diagonal='crossed')
        
        self.path = './../circle_in_a_box_2D/'
        
        
        
