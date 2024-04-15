#!/usr/bin/python

from fenics import *
from GeometrySimpleCubicMesh3D import SimpleCubicMesh3D
from PeriodicBoundary import *

# DOMAINS ########################################################

# TODO can be removed
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
        
        middle = self.ds * 0.5
        r = self.sc* 0.5
        
        if (x[0] >= (middle - r - self.tol) and x[0] <= (middle + r + self.tol) and 
            x[1] >= (middle - r - self.tol) and x[1] <= (middle + r + self.tol)):
            return True
        elif (x[0] >= (middle - r - self.tol) and x[0] <= (middle + r + self.tol) and 
              x[2] >= (middle - r - self.tol) and x[2] <= (middle + r + self.tol)):
            return True
        elif (x[1] >= (middle - r - self.tol) and x[1] <= (middle + r + self.tol) and 
              x[2] >= (middle - r - self.tol) and x[2] <= (middle + r + self.tol)):
            return True
        
        # corners
        elif (x[0] <= (r + self.tol) and x[1] <= (r + self.tol) and x[2] <= (r + self.tol)):
            return True
        elif (x[0] <= (r + self.tol) and x[1] <= (r + self.tol) and x[2] >= (self.ds - r - self.tol)):
            return True
        elif (x[0] <= (r + self.tol) and x[1] >= (self.ds - r - self.tol) and x[2] <= (r + self.tol)):
            return True
        elif (x[0] <= (r + self.tol) and x[1] >= (self.ds - r - self.tol) and x[2] >= (self.ds - r - self.tol)):
            return True
        elif (x[0] >= (self.ds - r - self.tol) and x[1] <= (r + self.tol) and x[2] <= (r + self.tol)):
            return True
        elif (x[0] >= (self.ds - r - self.tol) and x[1] <= (r + self.tol) and x[2] >= (self.ds - r - self.tol)):
            return True
        elif (x[0] >= (self.ds - r - self.tol) and x[1] >= (self.ds - r - self.tol) and x[2] <= (r + self.tol)):
            return True
        elif (x[0] >= (self.ds - r - self.tol) and x[1] >= (self.ds - r - self.tol) and x[2] >= (self.ds - r - self.tol)):
            return True
        else: 
            return False    
        """
        # edges   
        elif x[0] <= (r + self.tol) and x[1] <= (r + self.tol):
            return True
        elif x[0] <= (r + self.tol) and x[1] >= (self.ds - r - self.tol):
            return True
        elif x[0] >= (self.ds - r - self.tol) and x[1] <= (r + self.tol):
            return True
        elif x[0] >= (self.ds - r - self.tol) and x[1] >= (self.ds - r - self.tol):
            return True
        elif x[2] <= (r + self.tol) and x[1] <= (r + self.tol):
            return True
        elif x[2] <= (r + self.tol) and x[1] >= (self.ds - r - self.tol):
            return True
        elif x[2] >= (self.ds - r - self.tol) and x[1] <= (r + self.tol):
            return True
        elif x[2] >= (self.ds - r - self.tol) and x[1] >= (self.ds - r - self.tol):
            return True
        elif x[0] <= (r + self.tol) and x[2] <= (r + self.tol):
            return True
        elif x[0] <= (r + self.tol) and x[2] >= (self.ds - r - self.tol):
            return True
        elif x[0] >= (self.ds - r - self.tol) and x[2] <= (r + self.tol):
            return True
        elif x[0] >= (self.ds - r - self.tol) and x[2] >= (self.ds - r - self.tol):
            return True
        """
             
            
##############################################################

class Geometry(SimpleCubicMesh3D):

    def __init__(self, eps_l, eps_s, resolution=32, ds=16., sc=4., **kwargs):
        
        self.domain_size = ds
        self.solid_scale = sc
        
        p1 = Point(0., 0., 0.)
        p2 = Point(self.domain_size, self.domain_size, self.domain_size)
        
        self.tol=DOLFIN_EPS * 100.
        
        #self.subdomain_liquid = Liquid(self.tol)
        self.subdomain_solid = Solid(self.tol, self.domain_size, self.solid_scale)
        
        SimpleCubicMesh3D.__init__(self, p1, p2, resolution, eps_l, eps_s)
        
        self.path = './../actin_filaments_3d/'
        
        
        
