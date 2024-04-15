#!/usr/bin/python

from fenics import *

class Eps(UserExpression):
    def __init__(self, materials, eps_l, eps_s, **kwargs):
        self.materials = materials
        self.eps_l = eps_l
        self.eps_s = eps_s
        super().__init__(**kwargs)
        
    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.materials.mesh(), ufc_cell.index)
        if self.materials[cell.index()] == 2:
            values[0] = self.eps_s
        else:
            values[0] = self.eps_l

    def value_shape(self):
        return ()  # scalar elements
        

