#!/usr/bin/python

from fenics import *
from GeometrySimpleCubicMesh3D import SimpleCubicMesh3D
from PeriodicBoundary import *

# DOMAINS ########################################################

#    TODO: CAN BE REMOVED
"""      
class Liquid(SubDomain):

    def __init__(self, tol, **kwargs):
        self.tol = tol
        super().__init__(**kwargs)
  
    def inside(self, x, on_boundary):
        
        if (x[0] >= (0. - self.tol) and x[0] <= (5. + self.tol) and 
            x[1] >= (0. - self.tol) and x[1] <= (18. + self.tol)):
            return True
        elif (x[0] >= (5. - self.tol) and x[0] <= (13. + self.tol) and 
              x[1] >= (0. - self.tol) and x[1] <= (5. + self.tol)):
            return True
        elif (x[0] >= (5. - self.tol) and x[0] <= (13. + self.tol) and 
              x[1] >= (13. - self.tol) and x[1] <= (18. + self.tol)):
            return True
        elif (x[0] >= (13. - self.tol) and x[0] <= (18. + self.tol) and 
              x[1] >= (0. - self.tol) and x[1] <= (18. + self.tol)):
            return True
        else: 
            return False
"""
 
class Solid(SubDomain):

    def __init__(self, tol, ds, sc, **kwargs):
        self.tol = tol
        self.ds = ds
        self.sc = sc
        super().__init__(**kwargs)
        
    def inside(self, x, on_boundary):
        
        r = 0.5 * self.ds
        m = 0.5 * self.sc
        
        if (x[0] >= (r - m - self.tol) and x[0] <= (r + m + self.tol) and 
            x[1] >= (r - m - self.tol) and x[1] <= (r + m + self.tol)):
            return True
        elif (x[0] >= (r - m - self.tol) and x[0] <= (r + m + self.tol) and 
              x[2] >= (r - m - self.tol) and x[2] <= (r + m + self.tol)):
            return True
        elif (x[1] >= (r - m - self.tol) and x[1] <= (r + m + self.tol) and 
              x[2] >= (r - m - self.tol) and x[2] <= (r + m + self.tol)):
            return True
        else: 
            return False         
            
##############################################################

class Geometry(SimpleCubicMesh3D):

    def __init__(self, eps_l, eps_s, resolution=18, ds=18., sc=8., **kwargs):
        
        self.domain_size = ds
        self.solid_scale = sc
        
        p1 = Point(0., 0., 0.)
        p2 = Point(self.domain_size, self.domain_size, self.domain_size)
        
        self.tol = DOLFIN_EPS * 100.
        
        # self.subdomain_liquid = Liquid(self.tol)
        self.subdomain_solid = Solid(tol = self.tol, ds=self.domain_size, sc=self.solid_scale)
        
        SimpleCubicMesh3D.__init__(self, p1, p2, resolution, eps_l, eps_s)
        
        self.path = './../crossing_filaments_3d/'
        
        
        
