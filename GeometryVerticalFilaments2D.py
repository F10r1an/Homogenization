#!/usr/bin/python

from fenics import *
from GeometrySimpleRectangularMesh2D import SimpleRectangularMesh2D
from PeriodicBoundary import *

# DOMAINS ########################################################
        
class Liquid(SubDomain):

    def __init__(self, tol, sc, **kwargs):
        self.tol = tol
        self.sc = 0.5*sc
        super().__init__(**kwargs)
  
    def inside(self, x, on_boundary):
        
        if (x[0] >= (0. - self.tol) and x[0] <= (36. + self.tol) and 
            x[1] >= (0. - self.tol) and x[1] <= (9 - self.sc + self.tol)):
            return True
        elif (x[0] >= (0. - self.tol) and x[0] <= (9. - self.sc + self.tol) and 
              x[1] >= (9. - self.sc - self.tol) and x[1] <= (9. + self.sc + self.tol)):
            return True
        elif (x[0] >= (27. + self.sc - self.tol) and x[0] <= (36. + self.tol) and 
              x[1] >= (9. - self.sc - self.tol) and x[1] <= (9. + self.sc + self.tol)):
            return True
        elif (x[0] >= (0. - self.tol) and x[0] <= (36. + self.tol) and 
              x[1] >= (9. + self.sc - self.tol) and x[1] <= (27. - self.sc + self.tol)):
            return True
        elif (x[0] >= (9. + self.sc - self.tol) and x[0] <= (27. - self.sc + self.tol) and 
              x[1] >= (27. - self.sc - self.tol) and x[1] <= (27. + self.sc + self.tol)):
            return True
        elif (x[0] >= (0. - self.tol) and x[0] <= (36. + self.tol) and 
              x[1] >= (27. + self.sc - self.tol) and x[1] <= (36. + self.tol)):
            return True
        else: 
            return False
 
class Solid(SubDomain):

    def __init__(self, tol, sc, **kwargs):
        self.tol = tol
        self.sc = 0.5 * sc
        super().__init__(**kwargs)
        
    def inside(self, x, on_boundary):
        
        if (x[0] >= (9. - self.sc - self.tol) and x[0] <= (27. + self.sc + self.tol) and 
            x[1] >= (9. - self.sc - self.tol) and x[1] <= (9. + self.sc + self.tol)):
            return True
        elif (x[0] >= (0. - self.tol) and x[0] <= (9. + self.sc + self.tol) and 
            x[1] >= (27 - self.sc - self.tol) and x[1] <= (27. + self.sc + self.tol)):
            return True
        elif (x[0] >= (27. - self.sc - self.tol) and x[0] <= (36. + self.tol) and 
            x[1] >= (27. - self.sc - self.tol) and x[1] <= (27. + self.sc + self.tol)):
            return True
        else: 
            return False         

##############################################################

class Geometry(SimpleRectangularMesh2D):

    def __init__(self, eps_l, eps_s, resolution=144, ds=36., sc=8., **kwargs):
        
        self.domain_size = ds
        self.solid_scale = sc  # thickness of filemantes 
        
        p1 = Point(0., 0.)
        p2 = Point(36., 36.)
        
        self.tol=DOLFIN_EPS * 100.
        
        self.subdomain_liquid = Liquid(self.tol, sc=self.solid_scale)
        self.subdomain_solid = Solid(self.tol, sc=self.solid_scale)
        
        SimpleRectangularMesh2D.__init__(self, p1, p2, resolution, eps_l, eps_s)
            
        self.path = './../vertical_filaments_2d/'
        

        
        
        
