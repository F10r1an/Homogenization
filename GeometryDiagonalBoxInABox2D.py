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
        
        # always return True, let solid domain overwrite parts of liquid domain
        # in GeometrySimpleRectangularMesh.py
        return True
 
class Solid(SubDomain):

    def __init__(self, tol, ds, bs, **kwargs):
        self.tol = tol
        self.ds = np.sqrt(2.)* ds  # ds as in box_in_a_box_2d
        self.bs = np.sqrt(2.)*bs  # bs  as in box_in_a_box_2d
        super().__init__(**kwargs)
        
    def inside(self, x, on_boundary):
        
        #domain_size = np.sqrt(2.) * self.ds
        #solid_box = self.bs
        b = self.bs / 2.
        #c = domain_size - 2. * b
        
        if x[1] <= (-x[0] + b + self.tol):
            return True
        elif x[1] <= (x[0] + b - self.ds + self.tol):
            return True
        elif x[1] >= (x[0] + self.ds - b - self.tol):
            return True
        elif x[1] >= (-x[0] + 2*self.ds -b - self.tol):
            return True
        #elif ((x[1] >= (x[0] - 0.5*self.ds - b - self.tol)) and (x[1] <= (x[0] - 0.5*self.ds + b + self.tol))
        #      """and (x[1] >= (-x[0] +1.5*self.ds -b - self.tol)) and (x[1] <= (-x[0] +1.5*self.ds +b + self.tol))"""):
        #   return True
        elif ((x[1] >= (x[0]  - b - self.tol)) and (x[1] <= (x[0]  + b + self.tol))
              and (x[1] >= (-x[0] +self.ds -b - self.tol)) and (x[1] <= (-x[0] +self.ds +b + self.tol)) ):
            return True
        else: 
            return False         
            

##############################################################

class Geometry(SimpleRectangularMesh2D):

    def __init__(self, eps_l, eps_s, resolution=200, ds=np.sqrt(2)*18., sc=8., **kwargs):
        
        self.domain_size = ds
        self.solid_scale = sc
        
        self.p1 = Point(0., 0.)
        self.p2 = Point(self.domain_size*np.sqrt(2), self.domain_size*np.sqrt(2))
        
        self.tol=DOLFIN_EPS * 100.
        
        self.subdomain_liquid = Liquid(self.tol, self.domain_size, self.solid_scale)
        self.subdomain_solid = Solid(self.tol, self.domain_size, self.solid_scale)
        
        SimpleRectangularMesh2D.__init__(self, self.p1, self.p2, resolution, eps_l, eps_s,
        mesh_diagonal='crossed')
        
        self.path = './../diagonal_box_in_a_box_2D/'
        
        
        
