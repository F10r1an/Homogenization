#!/usr/bin/python

from fenics import *
from GeometrySimpleRectangularMesh2D import SimpleRectangularMesh2D
from PeriodicBoundary import *

# DOMAINS ########################################################
        
class Liquid(SubDomain):

    def __init__(self, tol, **kwargs):
        self.tol = tol
        super().__init__(**kwargs)
  
    def inside(self, x, on_boundary):
        
        if (x[0] >= (0. - self.tol) and x[0] <= (4. + self.tol) and 
            x[1] >= (1.2 - self.tol) and x[1] <= (3.6 + self.tol)):
            return True
        elif (x[0] >= (2. - self.tol) and x[0] <= (4. + self.tol) and 
              x[1] >= (3.6 - self.tol) and x[1] <= (7.6 + self.tol)):
            return True
        elif (x[0] >= (4. - self.tol) and x[0] <= (5. + self.tol) and 
              x[1] >= (5.2 - self.tol) and x[1] <= (7.6 + self.tol)):
            return True
        elif (x[0] >= (5. - self.tol) and x[0] <= (7. + self.tol) and 
              x[1] >= (3.6 - self.tol) and x[1] <= (7.6 + self.tol)):
            return True
        elif (x[0] >= (5. - self.tol) and x[0] <= (9. + self.tol) and 
              x[1] >= (1.2 - self.tol) and x[1] <= (3.6 + self.tol)):
            return True
        else: 
            return False
 
class Solid(SubDomain):

    def __init__(self, tol, **kwargs):
        self.tol = tol
        super().__init__(**kwargs)
        
    def inside(self, x, on_boundary):
        
        if (x[0] >= (0. - self.tol) and x[0] <= (9. + self.tol) and 
            x[1] >= (0. - self.tol) and x[1] <= (1.2 + self.tol)):
            return True
        elif (x[0] >= (4. - self.tol) and x[0] <= (5. + self.tol) and 
            x[1] >= (1.2 - self.tol) and x[1] <= (5.2 + self.tol)):
            return True
        elif (x[0] >= (0. - self.tol) and x[0] <= (2. + self.tol) and 
            x[1] >= (3.6 - self.tol) and x[1] <= (7.6 + self.tol)):
            return True
        elif (x[0] >= (7. - self.tol) and x[0] <= (9. + self.tol) and 
            x[1] >= (3.6 - self.tol) and x[1] <= (7.6 + self.tol)):
            return True
        elif (x[0] >= (0. - self.tol) and x[0] <= (9. + self.tol) and 
            x[1] >= (7.6 - self.tol) and x[1] <= (9. + self.tol)):
            return True
        else: 
            return False         
            
##############################################################

class Geometry(SimpleRectangularMesh2D):

    def __init__(self, eps_l, eps_s, resolution=18, res_y=22):
        
        self.tol = DOLFIN_EPS * 100.
        
        self.domain_size= [9.,8.8]
        self.solid_scale='undefined'
        
        self.p1 = Point(0., 0.)
        self.p2 = Point(9., 8.8)
        
        self.tol=DOLFIN_EPS * 100.
        
        self.subdomain_liquid = Liquid(self.tol)
        self.subdomain_solid = Solid(self.tol)
        
        SimpleRectangularMesh2D.__init__(self, self.p1, self.p2, resolution, eps_l, eps_s, res_y)
        
        self.path = './../perturbed_channel_2D/'
        
        # apply periodic boudary only in x-direction
        self.liquid_mesh_pb = PeriodicXBoundary2D(self.liquid_mesh)
        
        
        
