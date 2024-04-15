#!/usr/bin/python

import numpy as np
from fenics import *

class PeriodicBoundary1D(SubDomain):
    def __init__(self, mesh, tol = 1.e-7, **kwargs):
        super().__init__(**kwargs)
        self.mesh = mesh
        self.x_min, self.x_max = np.min(self.mesh.coordinates().transpose()[0]), np.max(self.mesh.coordinates().transpose()[0])
        self.Lx = self.x_max - self.x_min
        self.tol = tol

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        # return True if on left boundary 
        return bool(near(x[0], self.x_min, self.tol) and on_boundary)
     
    def map(self, x, y):
        if near(x[0], self.x_max, self.tol):
            y[0] = x[0] - self.Lx
        else:   # not on right border
            y[0] = x[0] - 10000. * self.Lx
            

class PeriodicBoundary2D(SubDomain):
    
    def __init__(self, mesh, tol = 1.e-7, **kwargs):
        super().__init__(**kwargs)
        self.mesh = mesh
        self.x_min, self.x_max = np.min(self.mesh.coordinates().transpose()[0]), np.max(self.mesh.coordinates().transpose()[0])
        self.y_min, self.y_max = np.min(self.mesh.coordinates().transpose()[1]), np.max(self.mesh.coordinates().transpose()[1])
        self.Lx = self.x_max - self.x_min
        self.Ly = self.y_max - self.y_min
        self.tol = tol
        
    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two corners (0, 1) and (1, 0)
        return bool((near(x[0], self.x_min) or near(x[1], self.y_min)) and 
                (not ((near(x[0], self.x_min) and near(x[1], self.y_max)) or 
                        (near(x[0], self.x_max) and near(x[1], self.y_min)))) and on_boundary)

    def map(self, x, y):
        
        if near(x[0], self.x_max, self.tol) and near(x[1], self.y_max, self.tol):
            y[0] = x[0] - self.Lx
            y[1] = x[1] - self.Ly
        elif near(x[0], self.x_max, self.tol):
            y[0] = x[0] - self.Lx
            y[1] = x[1]
        else:   # near(x[1], 1)
            y[0] = x[0]
            y[1] = x[1] - self.Ly
       
class PeriodicXBoundary2D(SubDomain):
    def __init__(self, mesh, tol = 1.e-7, **kwargs):
            super().__init__(**kwargs)
            self.mesh = mesh
            self.x_min, self.x_max = np.min(self.mesh.coordinates().transpose()[0]), np.max(self.mesh.coordinates().transpose()[0])
            self.Lx = self.x_max - self.x_min
            self.tol = tol

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        # return True if on left boundary 
        return bool(near(x[0], self.x_min, self.tol) and on_boundary)

    def map(self, x, y):
        if near(x[0], self.x_max, self.tol):
            y[0] = x[0] - self.Lx
            y[1] = x[1]
        else:   # not on right border
            y[0] = x[0] - 100000. * self.Lx
            y[1] = x[1]
       
class PeriodicBoundary3D(SubDomain):

    def __init__(self, mesh, **kwargs):
        super().__init__(**kwargs)
        self.mesh = mesh
        # coordinates of corners of cuboid
        self.x_min, self.x_max = np.min(self.mesh.coordinates().transpose()[0]), np.max(self.mesh.coordinates().transpose()[0])
        self.y_min, self.y_max = np.min(self.mesh.coordinates().transpose()[1]), np.max(self.mesh.coordinates().transpose()[1])
        self.z_min, self.z_max = np.min(self.mesh.coordinates().transpose()[2]), np.max(self.mesh.coordinates().transpose()[2])
        self.Lx = self.x_max - self.x_min
        self.Ly = self.y_max - self.y_min
        self.Lz = self.z_max - self.z_min

    def inside(self, x, on_boundary):
        # return True if on back (x=0) or left (y=0) or bottom (z=0) boundary AND NOT on one of the slave edges
        # master edges are: 1. (x=0, y=0): z_axis, 2. (y=0, z=0): x_axis, 3. (x=0, z=0): y_axis
        # 1. slave edges of (x=0, y=0) are: (x=0, y=1) and (x=1, y=0)
        # 2. slave edges of (y=0, z=0) are: (y=0, z=1) and (y=1, z=0)
        # 3. slave edges of (x=0, z=0) are: (x=0, z=1) and (x=1, z=0)
        # remark: the upper description is true for a unit box with egde length 1 where  x, y, z in [0, 1],
        # however the function works with cuboids of any size
        return bool(
            (near(x[0], self.x_min) or near(x[1], self.y_min) or near(x[2], self.z_min)) and 
            (not (
                (near(x[0], self.x_min) and near(x[1], self.y_max)) or 
                (near(x[0], self.x_max) and near(x[1], self.y_min)) or
                (near(x[1], self.y_min) and near(x[2], self.z_max)) or
                (near(x[1], self.y_max) and near(x[2], self.z_min)) or
                (near(x[0], self.x_min) and near(x[2], self.z_max)) or 
                (near(x[0], self.x_max) and near(x[2], self.z_min))
                 )) and on_boundary)
    
    def map(self, x, y):
        if near(x[0], self.x_max) and near(x[1], self.y_max) and near(x[2], self.z_max):  # map (1,1,1)  to (0,0,0) 
            y[0] = x[0] - self.Lx
            y[1] = x[1] - self.Ly 
            y[2] = x[2] - self.Lz
        elif near(x[0], self.x_max) and near(x[1], self.y_max):  # map (x=1, y=1) to z_axis
            y[0] = x[0] - self.Lx
            y[1] = x[1] - self.Ly 
            y[2] = x[2]
        elif near(x[1], self.y_max) and near(x[2], self.z_max):  # map (y=1, z=1) to x_axis
            y[0] = x[0]
            y[1] = x[1] - self.Ly 
            y[2] = x[2] - self.Lz
        elif near(x[0], self.x_max) and near(x[2], self.z_max):  # map (x=1, z=1) to y_axis
            y[0] = x[0] - self.Lx
            y[1] = x[1]  
            y[2] = x[2] - self.Lz
        elif near(x[0], self.x_max):  # map x=1 to x=0
            y[0] = x[0] - self.Lx
            y[1] = x[1]
            y[2] = x[2]
        elif near(x[1], self.y_max):  # map y=1 to y=0
            y[0] = x[0]
            y[1] = x[1] - self.Ly
            y[2] = x[2]
        else: #  near(x[1], y_max), map z=1 to z=0
            y[0] = x[0]
            y[1] = x[1] 
            y[2] = x[2] - self.Lz
