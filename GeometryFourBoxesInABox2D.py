#!/usr/bin/python

from fenics import *
from GeometrySimpleRectangularMesh2D import SimpleRectangularMesh2D
from PeriodicBoundary import *

# DOMAINS ########################################################
        
class Liquid(SubDomain):

    def __init__(self, tol, ds, bs, **kwargs):
        self.tol = tol
        self.ds = ds
        self.bs = bs
        super().__init__(**kwargs)
  
    def inside(self, x, on_boundary):
        
        """
        # Box 1    
        if (x[0] >= (0. - self.tol) and x[0] <= (2.5 + self.tol) and 
            x[1] >= (0. - self.tol) and x[1] <= (18. + self.tol)):
            return True
        elif (x[0] >= (6.5 - self.tol) and x[0] <= (11.5 + self.tol) and 
              x[1] >= (0. - self.tol) and x[1] <= (18. + self.tol)):
            return True
        elif (x[0] >= (15.5 - self.tol) and x[0] <= (18. + self.tol) and 
              x[1] >= (0. - self.tol) and x[1] <= (18. + self.tol)):
            return True
        elif (x[0] >= (0. - self.tol) and x[0] <= (18. + self.tol) and 
              x[1] >= (0. - self.tol) and x[1] <= (2.5 + self.tol)):
            return True
        elif (x[0] >= (0. - self.tol) and x[0] <= (18. + self.tol) and 
              x[1] >= (6.5 - self.tol) and x[1] <= (11.5 + self.tol)):
            return True
        elif (x[0] >= (0. - self.tol) and x[0] <= (18. + self.tol) and 
              x[1] >= (15. - self.tol) and x[1] <= (18. + self.tol)):
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
        self.ds = ds*2.  # domain size
        self.bs = bs # box size
        super().__init__(**kwargs)
        
    def inside(self, x, on_boundary):
        if (x[0] >= (0.25*self.ds-0.5*self.bs - self.tol) and x[0] <= (0.25*self.ds+0.5*self.bs + self.tol) and 
            x[1] >= (0.25*self.ds-0.5*self.bs - self.tol) and x[1] <= (0.25*self.ds+0.5*self.bs + self.tol)):
            return True  # BOX 1       
        elif (x[0] >= (0.25*self.ds-0.5*self.bs - self.tol) and x[0] <= (0.25*self.ds+0.5*self.bs + self.tol) and 
              x[1] >= (0.75*self.ds-0.5*self.bs - self.tol) and x[1] <= (0.75*self.ds+0.5*self.bs + self.tol)):
            return True  # BOX 2
        elif (x[0] >= (0.75*self.ds-0.5*self.bs - self.tol) and x[0] <= (0.75*self.ds+0.5*self.bs + self.tol) and 
              x[1] >= (0.25*self.ds-0.5*self.bs - self.tol) and x[1] <= (0.25*self.ds+0.5*self.bs + self.tol)):
            return True  # BOX 3
        elif (x[0] >= (0.75*self.ds-0.5*self.bs - self.tol) and x[0] <= (0.75*self.ds+0.5*self.bs + self.tol) and 
              x[1] >= (0.75*self.ds-0.5*self.bs - self.tol) and x[1] <= (0.75*self.ds+0.5*self.bs + self.tol)):
            return True  # BOX 4
        else: 
            return False         
            

##############################################################

class Geometry(SimpleRectangularMesh2D):

    def __init__(self, eps_l, eps_s, resolution=288, ds=18., sc=8., **kwargs):
        """
        ds : domain size
        bs : box size
        """
        self.domain_size = ds
        self.solid_scale = sc
        
        p1 = Point(0., 0.)
        p2 = Point(2.*self.domain_size, 2.*self.domain_size)
        
        self.tol=DOLFIN_EPS * 100.
        
        self.subdomain_liquid = Liquid(self.tol, self.domain_size, self.solid_scale)
        self.subdomain_solid = Solid(self.tol, self.domain_size, self.solid_scale)
        
        SimpleRectangularMesh2D.__init__(self, p1, p2, resolution, eps_l, eps_s)
        
        self.path = './../four_boxes_in_a_box_2D/'
        
        
        
