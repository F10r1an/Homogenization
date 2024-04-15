#!/usr/bin/python

from fenics import *
from Permittivity import Eps

class Solid(SubDomain):

    def __init__(self, N, tol, **kwargs):
        self.N = N  # number of reference cells in x-direction
        self.tol = tol
        
        super().__init__(**kwargs)
        
    def inside(self, x, on_boundary):
        
        x1 = x[0] % 9.
        x2 = x[1]
        
        if (x1 >= (0. - self.tol) and x1 <= (9. + self.tol) and 
            x2 >= (0. - self.tol) and x2 <= (1.2 + self.tol)):
            return True
        elif (x1 >= (4. - self.tol) and x1 <= (5. + self.tol) and 
            x2 >= (1.2 - self.tol) and x2 <= (5.2 + self.tol)):
            return True
        elif (x1 >= (0. - self.tol) and x1 <= (2. + self.tol) and 
            x2 >= (3.6 - self.tol) and x2 <= (7.6 + self.tol)):
            return True
        elif (x1 >= (7. - self.tol) and x1 <= (9. + self.tol) and 
            x2 >= (3.6 - self.tol) and x2 <= (7.6 + self.tol)):
            return True
        elif (x1 >= (0. - self.tol) and x1 <= (9. + self.tol) and 
            x2 >= (7.6 - self.tol) and x2 <= (9. + self.tol)):
            return True
        else: 
            return False     
            



"""
TOL = DOLFIN_EPS * 100.
N = 1  # number of reference cells in x-direction
res_scale = 1  # scale factor increases number of gridpoints in rectangular mesh
nx = 18 * N * res_scale  # number of grid points in x-direction
ny = 22 * res_scale  # number of grid points in y-direction
eps_l = 80.
eps_s = 1.
path_to_folder = './../verify_homogenization/'

# Create mesh and define funcion space
mesh = RectangleMesh(Point(0.,0.), Point(9.*N,8.8), nx, ny)
File(path_to_folder + "mesh.pvd") << mesh
V = FunctionSpace(mesh, 'CG', 1)

# mark solid and liquid subdomains of mesh
subdomain_solid = Solid(N, TOL)
materials = MeshFunction("size_t", mesh, mesh.geometric_dimension())
materials.set_all(1)  # liquid
subdomain_solid.mark(materials, 2)  # overwrite solid domainv
# permittivity of subdomains
epsilon = Eps(materials, eps_l, eps_s, degree=1) 
File(path_to_folder + "permittivty.pvd") << project(epsilon, FunctionSpace(mesh, 'DG', 0))


# Create Facetfunction to mark boundaries
boundary_markers = MeshFunction("size_t", mesh, mesh.geometric_dimension()-1)
bnd_left = DirichletBoundary(0., TOL)
bnd_left.mark(boundary_markers, 0)
bnd_right = DirichletBoundary(9. * N, TOL)
bnd_right.mark(boundary_markers, 1)
bnd_upper = NeumannBoundary(0., TOL)
bnd_upper.mark(boundary_markers, 2)
bnd_lower = NeumannBoundary(8.8, TOL)
bnd_lower.mark(boundary_markers, 3)

# DIrichlet boundary conditions
# Define Dirichlet boundary conditions at top and bottom boundaries
bcs = [DirichletBC(V, 5.0, boundary_markers, 0),
       DirichletBC(V, 0.0, boundary_markers, 1)]


u = TrialFunction(V)
v = TestFunction(V)

ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
dx = Measure('dx', domain=mesh, subdomain_data=materials)


g_U = Expression("- 10*exp(- pow(x[0] - 0.5, 2))", degree=2)
g_L = Constant(1.0)
f = Constant(1.0)


# Define variational form
F = (inner(eps_l*grad(u), grad(v))*dx(1) + inner(eps_s*grad(u), grad(v))*dx(2)
     - g_L*v*ds(2) - g_U*v*ds(3))
      

# Separate left and right hand sides of equation
a, L = lhs(F), rhs(F)

# Solve problem
w = Function(V)
solve(a == L, w, bcs)

File(path_to_folder + "solution.pvd") << w
"""


from dolfin import *

# Create classes for defining parts of the boundaries and the interior
# of the domain
class Left(SubDomain):
    
    def __init__(self, tol, **kwargs):
        self.tol = tol
        super().__init__(**kwargs)          
    
    def inside(self, x, on_boundary):
        return near(x[0], 0.0, self.tol)

class Right(SubDomain):

    def __init__(self, tol, N, **kwargs):
        self.tol = tol
        self.N = N
        super().__init__(**kwargs)    
    
    def inside(self, x, on_boundary):
        return near(x[0], N*9., self.tol)

class Bottom(SubDomain):

    def __init__(self, tol, **kwargs):
        self.tol = tol
        super().__init__(**kwargs)          
    
    def inside(self, x, on_boundary):
        return near(x[1], 0.0, self.tol)

class Top(SubDomain):
    
    def __init__(self, tol, **kwargs):
        self.tol = tol
        super().__init__(**kwargs)      
    
    def inside(self, x, on_boundary):
        return near(x[1], 8.8, self.tol)

class Obstacle(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[1], (3.5, 4.7)) and between(x[0], (2.2, 5.0)))

TOL = DOLFIN_EPS * 100.
N = 100  # number of reference cells in x-direction
res_scale = 1  # scale factor increases number of gridpoints in rectangular mesh
nx = 18 * N * res_scale  # number of grid points in x-direction
ny = 22 * res_scale  # number of grid points in y-direction
eps_l = 80.
eps_s = 1.
path_to_folder = './../verify_homogenization/'


# Initialize sub-domain instances
left = Left(TOL)
top = Top(TOL)
right = Right(TOL, N)
bottom = Bottom(TOL)
#obstacle = Obstacle()
solid = Solid(N, TOL)

# Define mesh
mesh = RectangleMesh(Point(0.,0.), Point(9.*N,8.8), nx, ny)

# Initialize mesh function for interior domains
domains = MeshFunction("size_t", mesh, mesh.topology().dim())
domains.set_all(0)
solid.mark(domains, 1)

# Initialize mesh function for boundary domains
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)

###########################
# mark solid and liquid subdomains of mesh
subdomain_solid = Solid(N, TOL)
materials = MeshFunction("size_t", mesh, mesh.geometric_dimension())
materials.set_all(1)  # liquid
subdomain_solid.mark(materials, 2)  # overwrite solid domainv
# permittivity of subdomains
epsilon = Eps(materials, eps_l, eps_s, degree=1) 
File(path_to_folder + "permittivty.pvd") << project(epsilon, FunctionSpace(mesh, 'DG', 0))
##############################

# Define input data
g_L = Expression("- 1*exp(- pow(x[1] - 0.5, 2))", degree=2)
g_R = Constant(0.0)
f = Constant(0.0)

# Define function space and basis functions
V = FunctionSpace(mesh, "CG", 2)
u = TrialFunction(V)
v = TestFunction(V)

# Define Dirichlet boundary conditions at top and bottom boundaries
bcs = [DirichletBC(V, 1.0, boundaries, 1),
       DirichletBC(V, 0.0, boundaries, 3)]

# Define new measures associated with the interior domains and
# exterior boundaries
dx = Measure('dx', domain=mesh, subdomain_data=domains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define variational form
F = (inner(eps_l*grad(u), grad(v))*dx(0) + inner(eps_s*grad(u), grad(v))*dx(1)
     - g_R*v*ds(2) - g_R*v*ds(4)
     - f*v*dx(0) - f*v*dx(1))

# Separate left and right hand sides of equation
a, L = lhs(F), rhs(F)

# Solve problem
u = Function(V)
solve(a == L, u, bcs)

File(path_to_folder + "solution.pvd") << u

# Evaluate integral of normal gradient over top boundary
n = FacetNormal(mesh)
m1 = dot(grad(u), n)*ds(1)
v1 = assemble(m1) * N
print("\int grad(u) * n ds(2) = ", v1)

#W = VectorFunctionSpace(mesh, "CG", 2)
#e_field = project(grad(u), W)
#File(path_to_folder + "field.pvd") << e_field


