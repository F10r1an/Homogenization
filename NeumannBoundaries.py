#!/usr/bin/python

from fenics import *

# BOUNDARIES ##################################################
# Notaion for boundarys:
# * first two letter are given by adjacent domains
#  ~ LL: liquid - liquid
#  ~ LS: liquid - solid

class LLBoundaries2D(SubDomain):
    
    def __init__(self, x0, y0, x1, y1, tol, **kwargs):
        self.x_min = x0
        self.y_min = y0
        self.x_max = x1
        self.y_max = y1
        self.tol = tol
        super().__init__(**kwargs)

    def inside(self, x, on_boundary):  
        return (
                (near(x[0], self.x_min, self.tol) or near(x[0], self.x_max, self.tol) or
                 near(x[1], self.y_min, self.tol) or near(x[1], self.y_max, self.tol)
                 ) and
                on_boundary)
                              
class LSBoundaries2D(SubDomain):
    
    def __init__(self, x0, y0, x1, y1, tol, **kwargs):
        self.x_min = x0
        self.y_min = y0
        self.x_max = x1
        self.y_max = y1
        self.tol = tol
        super().__init__(**kwargs)

    def inside(self, x, on_boundary):  
        return ((x[0] > self.x_min + self.tol) and (x[1] > self.y_min + self.tol) and
                (x[0] < self.x_max - self.tol) and (x[1] < self.y_max - self.tol) and
                on_boundary)

class LLBoundaries3D(SubDomain):
    
    def __init__(self, x0, y0, z0, x1, y1, z1, tol, **kwargs):
        self.x_min = x0
        self.y_min = y0
        self.z_min = z0
        self.x_max = x1
        self.y_max = y1
        self.z_max = z1
        self.tol = tol
        super().__init__(**kwargs)

    def inside(self, x, on_boundary):  
        return (
                (near(x[0], self.x_min, self.tol) or near(x[0], self.x_max, self.tol) or
                 near(x[1], self.y_min, self.tol) or near(x[1], self.y_max, self.tol) or
                 near(x[2], self.z_min, self.tol) or near(x[2], self.z_max, self.tol)
                 ) and
                on_boundary)
                              
class LSBoundaries3D(SubDomain):
    
    def __init__(self, x0, y0, z0, x1, y1, z1, tol, **kwargs):
        self.x_min = x0
        self.y_min = y0
        self.z_min = z0
        self.x_max = x1
        self.y_max = y1
        self.z_max = z1
        self.tol = tol
        super().__init__(**kwargs)

    def inside(self, x, on_boundary):  
        return ((x[0] > self.x_min + self.tol) and (x[1] > self.y_min + self.tol) and
                (x[0] < self.x_max - self.tol) and (x[1] < self.y_max - self.tol) and
                (x[2] > self.z_min + self.tol) and (x[2] < self.z_max - self.tol) and
                on_boundary)
