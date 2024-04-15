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
        box_min = 0.5 * self.ds - 0.5 * self.bs
        box_max = 0.5 * self.ds + 0.5 * self.bs   
        if (x[0] >= (0. - self.tol) and x[0] <= (box_min + self.tol) and 
            x[1] >= (0. - self.tol) and x[1] <= (self.ds + self.tol)):
            return True
        elif (x[0] >= (box_min - self.tol) and x[0] <= (box_max + self.tol) and 
              x[1] >= (0. - self.tol) and x[1] <= (box_min + self.tol)):
            return True
        elif (x[0] >= (box_min - self.tol) and x[0] <= (box_max + self.tol) and 
              x[1] >= (box_max - self.tol) and x[1] <= (self.ds + self.tol)):
            return True
        elif (x[0] >= (box_max - self.tol) and x[0] <= (self.ds + self.tol) and 
              x[1] >= (0. - self.tol) and x[1] <= (self.ds + self.tol)):
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
        box_min = 0.5 * self.ds - 0.5 * self.bs
        box_max = 0.5 * self.ds + 0.5 * self.bs                
        if (x[0] >= (box_min - self.tol) and x[0] <= (box_max + self.tol) and 
            x[1] >= (box_min - self.tol) and x[1] <= (box_max + self.tol)):
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
        self.bnd_coords_x1 = [[0., 0., 18.],[18., 0., 18.]]  # [[x1, x2_1, x2_2]...]
        
        # boundary is given by three numbers: x2, x1_1, x1_2
        # (x1_1, x2) x-------------------------x (x1_2, x2)      
        self.bnd_coords_x2 = [[0., 0., 18.],[18., 0., 18.]]  # [[x2, x1_1, x1_2]...]
        
        super().__init__(**kwargs)

    def inside(self, x, on_boundary):  
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

class LSBoundaries(SubDomain):
    
    def __init__(self, tol, **kwargs):
        self.tol = tol
        # boundary is given by three numbers: x1, x2_1, x2_2
        # (x1, x2_1) x-------------------------x (x1, x2_2)      
        self.bnd_coords_x1 = [[5., 5., 13.], [13., 5., 13.]]  # [[x1, x2_1, x2_2]...]
                                   
        # boundary is given by three numbers: x2, x1_1, x1_2
        # (x1_1, x2) x-------------------------x (x1_2, x2)      
        self.bnd_coords_x2 = [[5., 5., 13.], [13., 5., 13.]]  # [[x2, x1_1, x1_2]...]
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

    def __init__(self, eps_l, eps_s, resolution=72, ds=18., sc=8., **kwargs):
        """
        ds : domain size
        bs : box size
        """
        self.domain_size = ds
        self.solid_scale = sc
        
        p1 = Point(0., 0.)
        p2 = Point(self.domain_size, self.domain_size)
        
        self.tol=DOLFIN_EPS * 100.
        
        self.subdomain_liquid = Liquid(self.tol, self.domain_size, self.solid_scale)
        self.subdomain_solid = Solid(self.tol, self.domain_size, self.solid_scale)
        
        SimpleRectangularMesh2D.__init__(self, p1, p2, resolution, eps_l, eps_s)
        
        self.path = './../box_in_a_box_2D/'
        
        
        
