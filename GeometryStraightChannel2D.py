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
        self.ls = ds - bs # liquid size
        super().__init__(**kwargs)
  
    def inside(self, x, on_boundary):
        
        if (x[1] >= (0.5*self.ds - 0.5*self.ls - self.tol) and x[1] <= (0.5*self.ds + 0.5*self.ls + self.tol) ):
            return True
        else: 
            return False
 
class Solid(SubDomain):

    def __init__(self, tol, ds, bs, **kwargs):
        self.tol = tol
        self.ds = ds  # domain size
        self.bs = bs  # box size
        self.ls = ds - bs # liquid size
        super().__init__(**kwargs)
        
    def inside(self, x, on_boundary):
        
        if (x[1] >= (0. - self.tol) and x[1] <= (0.5*self.ds-0.5*self.ls + self.tol)):
            return True
        elif (x[1] >= (0.5*self.ds + 0.5*self.ls - self.tol) and x[1] <= (self.ds + self.tol)):
            return True
        else: 
            return False         
            
# BOUNDARIES ##################################################
# Notaion for boundarys:
# * first two letter are given by adjacent domains
#  ~ LL: liquid - liquid
#  ~ SS: solid - solid 
#  ~ LS: liquid - solid

class LLBoundaries(SubDomain):
    
    def __init__(self, tol, **kwargs):
        self.tol = tol
        # boundary is given by three numbers: x1, x2_1, x2_2
        # (x1, x2_1) x-------------------------x (x1, x2_2)      
        self.bnd_coords = [[0., 1., 5.5],[9., 1., 5.5]]  # [[x1, x2_1, x2_2]...]
        super().__init__(**kwargs)

    def inside(self, x, on_boundary):  
        for pnt in self.bnd_coords:
            if (near(x[0], pnt[0], self.tol) and x[1] > (pnt[1] - self.tol)
                                             and x[1] < (pnt[2] + self.tol) 
                ):
                return True     
        # if point is not on boundary:
        return False

class LSBoundaries(SubDomain):
    
    def __init__(self, tol, **kwargs):
        self.tol = tol
        # boundary is given by three numbers: x1, x2_1, x2_2
        # (x1, x2_1) x-------------------------x (x1, x2_2)      
        self.bnd_coords_x1 = []  # [[x1, x2_1, x2_2]...]
                                   
        # boundary is given by three numbers: x2, x1_1, x1_2
        # (x1_1, x2) x-------------------------x (x1_2, x2)      
        self.bnd_coords_x2 = [[1., 0., 9.], [5.5, 0., 9.]]  # [[x2, x1_1, x1_2]...]
        super().__init__(**kwargs)

    def inside(self, x, on_boundary):
        # test if point lies on ans LS boundary
        for pnt in self.bnd_coords_x1:
            if (near(x[0], pnt[0], self.tol) and x[1] > (pnt[1] - self.tol)
                                             and x[1] < (pnt[2] + self.tol) 
                ):
                return True    
        for pnt in self.bnd_coords_x2:
            if (near(x[1], pnt[0], self.tol) and x[0] > (pnt[1] - self.tol)
                                             and x[0] < (pnt[2] + self.tol) 
                ):
                return True      
        # if point is not on boundary:
        return False
        
class SSBoundaries(SubDomain):
    def __init__(self, tol, **kwargs):
        self.tol = tol
        # boundary is given by three numbers: x1, x2_1, x2_2
        # (x1, x2_1) x-------------------------x (x1, x2_2)      
        self.bnd_coords_x1 = []  # [[x1, x2_1, x2_2]...]
                                   
        # boundary is given by three numbers: x2, x1_1, x1_2
        # (x1_1, x2) x-------------------------x (x1_2, x2)      
        self.bnd_coords_x2 = []  # [[x2, x1_1, x1_2]...]
        super().__init__(**kwargs)


    def inside(self, x, on_boundary):
        # test if point lies on ans LS boundary
        for pnt in self.bnd_coords_x1:
            if (near(x[0], pnt[0], self.tol) and x[1] > (pnt[1] - self.tol)
                                             and x[1] < (pnt[2] + self.tol) 
                ):
                return True    
        for pnt in self.bnd_coords_x2:
            if (near(x[1], pnt[0], self.tol) and x[0] > (pnt[1] - self.tol)
                                             and x[0] < (pnt[2] + self.tol) 
                ):
                return True      
        # if point is not on boundary:
        return False

##############################################################

class Geometry(SimpleRectangularMesh2D):

    def __init__(self, eps_l, eps_s, resolution=90, ds=9., sc=3., **kwargs):
        
        self.domain_size = ds
        self.solid_scale = sc
        
        p1 = Point(0., 0.)
        p2 = Point(self.domain_size, self.domain_size)
        
        self.tol=DOLFIN_EPS * 100.
        
        self.subdomain_liquid = Liquid(self.tol, self.domain_size, self.solid_scale)
        self.subdomain_solid = Solid(self.tol, self.domain_size, self.solid_scale)
        
        SimpleRectangularMesh2D.__init__(self, p1, p2, resolution, eps_l, eps_s)
        self.path = './../straight_channel_2D/'
        
        # apply periodic boudary only in x-direction
        self.liquid_mesh_pb = PeriodicXBoundary2D(self.liquid_mesh)
        
