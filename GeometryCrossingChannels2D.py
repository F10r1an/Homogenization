#!/usr/bin/python

from fenics import *
from GeometrySimpleRectangularMesh2D import SimpleRectangularMesh2D
from PeriodicBoundary import *

# DOMAINS ########################################################
        
class Liquid(SubDomain):

    def __init__(self, tol, ds, bs, **kwargs):
        self.tol = tol
        self.ds = ds  # domain size
        self.bs = bs  # box size
        super().__init__(**kwargs)
  
    def inside(self, x, on_boundary):
        
        if (x[0] >= (0.5*self.bs - self.tol) and x[0] <= (self.ds-0.5*self.bs + self.tol) ):
            return True
        elif (x[1] >= (0.5*self.bs - self.tol) and x[1] <= (self.ds-0.5*self.bs + self.tol) ):
            return True
        else: 
            return False
 
class Solid(SubDomain):

    def __init__(self, tol, ds, bs, **kwargs):
        self.tol = tol
        self.ds = ds  # domain size
        self.bs = bs  # box size
        super().__init__(**kwargs)
        
    def inside(self, x, on_boundary):
        
        if (x[0] <= (0.5*self.bs + self.tol) and x[1] <= (0.5*self.bs + self.tol)):
            return True
        elif (x[0] >= (self.ds-0.5*self.bs - self.tol) and x[1] <= (0.5*self.bs + self.tol)):
            return True
        elif (x[0] <= (0.5*self.bs + self.tol) and x[1] >= (self.ds-0.5*self.bs - self.tol)):
            return True
        elif (x[0] >= (self.ds-0.5*self.bs - self.tol) and x[1] >= (self.ds-0.5*self.bs - self.tol)):
            return True
        else: 
            return False         

##############################################################

class Geometry(SimpleRectangularMesh2D):

    def __init__(self, eps_l, eps_s, resolution=36, ds=18., sc=4., **kwargs):
        
        self.domain_size = ds
        self.solid_scale = sc
        
        p1 = Point(0., 0.)
        p2 = Point(self.domain_size, self.domain_size)
        
        self.tol=DOLFIN_EPS * 100.
        
        self.subdomain_liquid = Liquid(self.tol, self.domain_size, self.solid_scale)
        self.subdomain_solid = Solid(self.tol, self.domain_size, self.solid_scale)
        
        SimpleRectangularMesh2D.__init__(self, p1, p2, resolution, eps_l, eps_s)
        
        self.path = './../crossing_channels_2D/'
        
       
