#!/usr/bin/python

from fenics import *
from GeometrySimpleRectangularMesh2D import SimpleRectangularMesh2D
from PeriodicBoundary import *

# DOMAINS ########################################################
        
class Liquid(SubDomain):

    def __init__(self, tol, ds, bs, **kwargs):
        self.tol = tol
        self.ds = ds * np.sqrt(2.)  # domain size
        self.bs = bs * np.sqrt(2.)  # box size
        self.ls = (ds - bs) * np.sqrt(2.) # liquid size
        super().__init__(**kwargs)
  
    def inside(self, x, on_boundary):
        """
        if (x[1] >= (x[0] + 2.5 - self.tol) and x[1] <= (x[0] + 7.5 + self.tol)):
            return True
        elif (x[1] >= (x[0] - 7.5 - self.tol) and x[1] <= (x[0] - 2.5 + self.tol)):
            return True
        
        else: 
            return False
        """
        # always return True, let solid domain overwrite parts of liquid domain
        # in GeometrySimpleRectangularMesh.py
        return True
 
class Solid(SubDomain):

    def __init__(self, tol, ds, bs, **kwargs):
        self.tol = tol
        self.ds = ds * np.sqrt(2.)  # domain size
        self.bs = bs * np.sqrt(2.)  # box size
        self.ls = (ds - bs) * np.sqrt(2.) # liquid size
        super().__init__(**kwargs)
        
    def inside(self, x, on_boundary):
        
        if x[1] >= (x[0] + 0.5*self.ds+0.5*self.ls - self.tol):
            return True
        elif (x[1] >= (x[0] - 0.5*self.ds+0.5*self.ls - self.tol) and x[1] <= (x[0]  + 0.5*self.ds-0.5*self.ls + self.tol)):
            return True  
        elif x[1] <= (x[0] - 0.5*self.ds-0.5*self.ls + self.tol):
            return True
        else: 
            return False
            
##############################################################

class Geometry(SimpleRectangularMesh2D):

    def __init__(self, eps_l, eps_s, resolution=200, ds=10., sc=5., **kwargs):
        
        self.domain_size = ds
        self.solid_scale = sc
        
        p1 = Point(0., 0.)
        p2 = Point(self.domain_size*np.sqrt(2.), self.domain_size*np.sqrt(2.))
        
        self.tol=DOLFIN_EPS * 100.
        
        self.subdomain_liquid = Liquid(self.tol, self.domain_size, self.solid_scale)
        self.subdomain_solid = Solid(self.tol, self.domain_size, self.solid_scale)
        
        SimpleRectangularMesh2D.__init__(self, p1, p2, resolution, eps_l, eps_s,
        mesh_diagonal='crossed')
        
        self.path = './../diagonal_channels_2D/'
        
        
        
        
