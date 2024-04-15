#!/usr/bin/python

import numpy as np
from fenics import *
from PeriodicBoundary import *

def reference_cell_solver_potential(geometry):
    
    Eps = geometry.get_permittivity()
    mesh = geometry.get_mesh()
    
    domain_markers = geometry.get_materials()
    
    # set periodic boundary condition
    dim = mesh.geometric_dimension()
    pbc = geometry.get_mesh_periodic_boundary()
    
    #if dim == 1:
    #    pbc = PeriodicBoundary1D(mesh)
    #elif dim == 2:
    #    pbc = PeriodicBoundary2D(mesh)
    #elif dim == 3:
    #    pbc = PeriodicBoundary3D(mesh)
    #else:
    #    raise Exception('Spatial dimension of mesh must be wrong! dim=',dim)
    
    u = []
    
    # Build function space with Lagrange multiplier
    P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    R = FiniteElement("Real", mesh.ufl_cell(), 0)
    W = FunctionSpace(mesh, P1 * R, constrained_domain=pbc)
    
    # unit vectors
    U = VectorFunctionSpace(mesh, 'Lagrange', 1)
    V = FunctionSpace(mesh, 'CG', 1)
    
    if dim == 1:
        raise Exception('not implemented yet for 1D case')
    elif dim == 2:
        C_x1 = Constant([1.,0.])   
        C_x2 = Constant([0.,1.]) 
 
    elif dim == 3:
        C_x1 = Constant([1.,0.,0.])  
        C_x2 = Constant([0.,1.,0.]) 
        C_x3 = Constant([0.,0.,1.])  
     
    # Define variational problem
    (u1, c1) = TrialFunction(W)
    (v1, d1) = TestFunctions(W)
    
    a1 = Eps * inner(grad(u1), grad(v1))*dx + c1*v1*dx + u1*d1*dx
    L1 = inner(grad(project(Eps,V)), C_x1) * v1 * dx 
    
    # Compute solution
    w1 = Function(W)
    solve(a1 == L1, w1)
    (u1, c1) = w1.split()
    
    u.append(u1)
    
    
    if dim > 1:
        # Define variational problem
        (u2, c2) = TrialFunction(W)
        (v2, d2) = TestFunctions(W)
        
        a2 = Eps * inner(grad(u2), grad(v2))*dx + c2*v2*dx + u2*d2*dx
        L2 = inner(grad(project(Eps,V)), C_x2)*v2*dx
        
        # Compute solution
        w2 = Function(W)
        solve(a2 == L2, w2)
        (u2, c2) = w2.split()
        
        u.append(u2)
     
    if dim > 2:

        # Define variational problem
        (u3, c3) = TrialFunction(W)
        (v3, d3) = TestFunctions(W)
        
        a3 = Eps * inner(grad(u3), grad(v3))*dx + c3*v3*dx + u3*d3*dx
        L3 = inner(grad(project(Eps,V)), C_x2)*v3*dx
        
        # Compute solution
        w3 = Function(W)
        solve(a3 == L3, w3)
        (u3, c3) = w3.split()
        
        u.append(u3) 

    return u
    
    
def reference_cell_solver_diffusion(geometry):
    
    mesh = geometry.get_liquid_mesh()
    boundary_markers = geometry.get_boundaries()
    dim = mesh.geometric_dimension()
    pbc = geometry.get_liquid_mesh_periodic_boundary()
    
    u = []

    # Build function space with Lagrange multiplier
    P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    R = FiniteElement("Real", mesh.ufl_cell(), 0)
    W = FunctionSpace(mesh, P1 * R, constrained_domain=pbc)
    
    ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
    dx = Measure('dx', domain=mesh)
    
    if dim == 1:
        raise Exception('not implemented yet for 1D case')
    elif dim == 2:
        c_1 = Constant([1.0, 0.0])
        c_2 = Constant([0.0, 1.0])
    elif dim == 3:
        c_1 = Constant([1.0, 0.0, 0.0])
        c_2 = Constant([0.0, 1.0, 0.0])
        c_3 = Constant([0.0, 0.0, 1.0])
    # Define variational problem
    (u1, c1) = TrialFunction(W)
    (v1, d1) = TestFunctions(W)
    
    n =FacetNormal(mesh)

    a1 = -inner(grad(u1), grad(v1))*dx + c1*v1*dx + u1*d1*dx   
    L1 = dot(c_1, n)* v1 * ds(2)
    
    # Compute solution
    w1 = Function(W)
    solve(a1 == L1, w1)
    (u1, c1) = w1.split()
    
    u.append(u1)
    
    if dim > 1:
        # Define variational problem
        (u2, c2) = TrialFunction(W)
        (v2, d2) = TestFunctions(W)
        
        a2 = -inner(grad(u2), grad(v2))*dx + c2*v2*dx + u2*d2*dx
        L2 = dot(c_2, n)*v1* ds(2)
        
        # Compute solution
        w2 = Function(W)
        solve(a2 == L2, w2)
        (u2, c2) = w2.split()
        
        u.append(u2)
     
    if dim > 2:
       # Define variational problem
        (u3, c3) = TrialFunction(W)
        (v3, d3) = TestFunctions(W)
        
        a3 = -inner(grad(u3), grad(v3))*dx + c3*v3*dx + u3*d3*dx
        L3 = dot(c_3, n)*v1* ds(2)
        
        # Compute solution
        w3 = Function(W)
        solve(a3 == L3, w3)
        (u3, c3) = w3.split()
        
        u.append(u3) 
    
    return u
    
    
def test_solution_periodicity(mesh, xi_phi, res=7, dev=1.e-7):
    """
    test periodicity of solution to guarantee that boundary condtions 
    are fullfilled.
    therefore evaluate function of grid on boundary surfaces and test
    for periodicity
    returns 1 if solution is periodic and 0 otherwise
    
    res: number of grid points in one direction (res^2 points on each surface)
    dev: accepted deviation of function values on opposite sides of mesh
    """
    # coordinates of grid nodes
    coords = mesh.coordinates().transpose()
    # spatial dimension of mesh domain
    dim = mesh.geometric_dimension()
    
    x1_min, x1_max = np.min(mesh.coordinates().transpose()[0]), np.max(mesh.coordinates().transpose()[0])
    x2_min, x2_max = np.min(mesh.coordinates().transpose()[1]), np.max(mesh.coordinates().transpose()[1])
    x3_min, x3_max = np.min(mesh.coordinates().transpose()[2]), np.max(mesh.coordinates().transpose()[2])
    x_min = [x1_min, x2_min, x3_min]
    x_max = [x1_max, x2_max, x3_max]
    print(x_min, x_max)
    
    if dim == 3:
        for cf in xi_phi:
            for a in range(dim):
                xa_min, xa_max = x_min[a], x_max[a]*0.999
                
                b = (a +1)%dim
                c = (a+2)%dim
                
                xb_min, xb_max = x_min[b], x_max[b]
                xc_min, xc_max = x_min[c], x_max[c]
                b_points = np.linspace(xb_min, xb_max, res+2)[1:-1]
                c_points = np.linspace(xc_min, xc_max, res+2)[1:-1]
                b_grid, c_grid = np.meshgrid(b_points, c_points)
                x_min_values = []
                x_max_values = []
                #from IPython import embed
                #embed()
                for i in range(res):
                    for j in range(res):
                        p1 = np.zeros(3)
                        p2 = np.zeros(3)
                        p1[a] = xa_min
                        p2[a] = xa_max
                        p1[b], p2[b] = b_grid[i, j], b_grid[i, j]
                        p1[c], p2[c] = c_grid[i, j], c_grid[i, j]
                        
                        x_min_values.append(cf(p1))
                        x_max_values.append(cf(p2))
                #from IPython import embed
                #embed()
                periodic = (np.all(np.abs(np.array(x_min_values)-np.array(x_max_values))<dev))
                if not periodic==True:
                    return 0
    else:
        raise Exception('2D case not yet implemented')
    
    return 1

def diffusion_tensor(xi_c, geometry):
    
    mesh = geometry.get_liquid_mesh()
    volume = geometry.get_liquid_mesh_volume()
    dx = Measure('dx', domain=mesh)
    dim = mesh.geometric_dimension()
    
    if dim == 2:
        unit_vectors = [Constant([1., 0.]), Constant([0., 1.])]
    elif dim == 3:
        unit_vectors = [Constant([1., 0., 0.]), Constant([0., 1., 0.]), Constant([0., 0., 1.])]
    
    tensor = np.zeros([dim,dim])
    
    for i in range(dim):
        for j in range(dim):
            print('Computing tensor component i,j = ', i,j)
            
            # compute delta_kl
            if i == j:
                delta_ij = Constant(1.)
            else:
                delta_ij = Constant(0.)
            
            # compute dx_k xi^l
            xi_j = xi_c[j]
            e_i = unit_vectors[i]

            tensor[i,j] = assemble((dot(grad(xi_j), e_i) + delta_ij)*dx) / volume
        
    return tensor     
    
def permittivity_tensor(xi_phi, geometry):
    
    mesh = geometry.get_mesh()
    domain_markers = geometry.get_materials()
    volume = geometry.get_mesh_volume()
    dx = Measure('dx', domain=mesh, subdomain_data=domain_markers)
    dim = mesh.geometric_dimension()
    eps = geometry.get_permittivity()
    
    if dim == 2:
        unit_vectors = [Constant([1., 0.]), Constant([0., 1.])]
    elif dim == 3:
        unit_vectors = [Constant([1., 0., 0.]), Constant([0., 1., 0.]), Constant([0., 0., 1.])]
    
    tensor = np.zeros([dim,dim])
    
    for i in range(dim):
        for j in range(dim):
            print('Computing tensor component i,j = ', i,j)
            
            # compute delta_ij
            if i == j:
                delta_ij = Constant(1.)
            else:
                delta_ij = Constant(0.)

            # compute dx_k xi^l
            xi_j = xi_phi[j]
            e_i = unit_vectors[i]

            tensor[i,j] = assemble((dot(
                                        eps , dot(grad(xi_j), e_i)+ delta_ij
                                       )) *dx) / volume
            #tensor[i,j] = assemble((dot(grad(xi_j), e_i) + delta_ij)*dx) / volume
            #tensor[i,j] = assemble(1.*(dot(grad(xi_j), e_i) + delta_ij)*dx(1)) / volume + assemble(2.*(dot(grad(xi_j), e_i) + delta_ij)*dx(2)) / volume 
            
    return tensor

def mobility_tensor(xi_phi, geometry):
    
    mesh = geometry.get_mesh()
    domain_markers = geometry.get_materials()
    volume = geometry.get_mesh_volume()
    dx = Measure('dx', domain=mesh, subdomain_data=domain_markers)
    dim = mesh.geometric_dimension()
    eps = geometry.get_permittivity()
    
    if dim == 2:
        unit_vectors = [Constant([1., 0.]), Constant([0., 1.])]
    elif dim == 3:
        unit_vectors = [Constant([1., 0., 0.]), Constant([0., 1., 0.]), Constant([0., 0., 1.])]
    
    tensor = np.zeros([dim,dim])
    
    for i in range(dim):
        for j in range(dim):
            print('Computing mobility tensor component i,j = ', i,j)
            
            # compute delta_ij
            if i == j:
                delta_ij = Constant(1.)
            else:
                delta_ij = Constant(0.)

            # compute dx_k xi^l
            xi_j = xi_phi[j]
            e_i = unit_vectors[i]
            
            tensor[i,j] = assemble((dot(grad(xi_j), e_i) + delta_ij)*dx(1)) / volume
    return tensor

