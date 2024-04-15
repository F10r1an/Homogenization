#!/usr/bin/python

import numpy as np
from fenics import *
from Permittivity import Eps

class GeometryBase():

    def __init__(self):
        
        # tolerance to postions/coordinates
        self.tol = None
        
        # permittivity of liquid and solid phase
        self.eps_l = None
        self.eps_s = None
        
        # full mesh over solid and liquid domains
        self.mesh = None
        
        # mark liquid and solid domains:
        self.materials = None
        
        self.epsilon = None

        # mesh over liquid domain
        self.liquid_mesh = None
                
        # boundaries of liquid mesh
        self.boundaries = None

        # periodic boundaries of domain
        self.mesh_pb = None
        self.liquid_mesh_pb = None
       
        # path where to save files
        self.path = None
 
    def save_boundaries(self, file_name='boundaries.pvd'):
        file = File(self.path + file_name)
        file << self.boundaries
    
    def save_materials(self, file_name='materials.pvd'):
        file = File(self.path + file_name)
        file << self.materials
     
    def get_permittivity(self):
        if self.epsilon == None:
            self.compute_permittivity()
            
        return self.epsilon
    
    def save_permittivity(self, file_name='permittivity.pvd'):
        
        eps = self.get_permittivity()
        V_temp = FunctionSpace(self.mesh, 'DG', 0)
        eps_interp = project(eps, V_temp)
                          
        file = File(self.path + file_name)
        file << eps_interp
        
    def save_reference_cell_solution_diffusion(self):
        pass
        
    def geometric_dimension(self):
        return self.mesh.geometric_dimension()
        
    def get_liquid_mesh_periodic_boundary(self):
        return self.liquid_mesh_pb
    
    def get_mesh_periodic_boundary(self):
        return self.mesh_pb
    
    def get_folder(self):
        return self.path 
        
    def compute_permittivity(self):
        self.epsilon = Eps(self.materials, self.eps_l, self.eps_s, degree=1) 
                                      
    def get_mesh(self):
       return self.mesh
       
    def get_liquid_mesh(self):
        return self.liquid_mesh
       
    def get_boundaries(self):
        return self.boundaries
        
    def get_materials(self):
        return self.materials
    
    def get_mesh_volume(self):
        """
        get mesh volume by integrating "1" over full mesh
        """
        dx = Measure('dx', domain=self.mesh)
        V = FunctionSpace(self.mesh, 'DG', 0)
        c = Constant(1.)
        cp = project(c, V)
        vol = assemble(cp * dx)
        return vol
     
    def get_liquid_mesh_volume(self):
        """
        get mesh volume by integrating "1" over full mesh
        """
        dx = Measure('dx', domain=self.liquid_mesh)
        V = FunctionSpace(self.liquid_mesh, 'DG', 0)
        c = Constant(1.)
        cp = project(c, V)
        vol = assemble(cp * dx)
        return vol
        
        
    def test_grid_periodicity(self, axis=''):
        """
        test periodicity of grid on boundarys for a rectancular (2D) or 
        cuboid (3D) mesh
        mesh: fenics mesh
        axis: periodicity along which axis to be tested
        returns 1 if test was successfull
        returns 0 if test was not successfull
        """
        
        mesh = self.mesh
        
        if axis == 'x' or axis == 'x1' or axis == 'x_1':
            axis_a = 0
        elif axis == 'y' or axis == 'x2' or axis == 'x_2':
            axis_a = 1
        elif axis == 'z' or axis == 'x3' or axis == 'x_3':
            axis_a = 2
        else:
            print('wrong input argument for "axis". Argument as to be "x_1", "x_2" or "x_3".')
            return 0
        coords = mesh.coordinates().transpose()
        dim = mesh.geometric_dimension()
        
        x_min, x_max = np.min(coords[axis_a]), np.max(coords[axis_a])
        eps = np.abs((x_max-x_min)*1.e-15)  # 1.e-15 : approx double floating point number
        x_min_indices = np.where(np.abs(coords[axis_a] - x_min) < eps)
       
        x_max_indices = np.where(np.abs(coords[axis_a] - x_max) < eps)
        
        if dim == 1: 
            return 1  # grid is periodic as it only contains 1 boundary point
        
        elif dim == 2:
            axis_b = (axis_a + 1) % dim  #  axis of the other spatial dimensions

            b_min_values = coords[axis_b][x_min_indices]  # corresponds to axis b
            b_max_values = coords[axis_b][x_max_indices]
            # test if coordinates are the same
            if np.all(  np.abs(np.sort(b_min_values, axis=0)-np.sort(b_max_values, axis=0)) < eps):
                return 1 
            else:
                return 0
        if dim == 3:
            axis_b = (axis_a + 1) % dim  #  axis of the oter spacial dimensions
            axis_c = (axis_a + 2) % dim  #

            
            b_min_values = coords[axis_b][x_min_indices]  # corresponds to axis b
            b_max_values = coords[axis_b][x_max_indices]
            c_min_values = coords[axis_c][x_min_indices]  # corresponds to axis b
            c_max_values = coords[axis_c][x_max_indices]
            
            min_points = np.array([b_min_values, c_min_values]).transpose()
            max_points = np.array([b_max_values, c_max_values]).transpose()
            
            # test if coordinates are the same
            if np.all(  np.abs(np.sort(min_points, axis=0)-np.sort(max_points, axis=0)) < eps):
                return 1 
            else:
                return 0
         
