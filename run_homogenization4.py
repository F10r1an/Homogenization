#!/usr/bin/python

import datetime
import pickle

from fenics import *
import numpy as np
import PeriodicBoundary as pb
import Homogenization as hm

#from GeometryBase import GeometryBase 
#gm = 'box_in__box_2d'
#gm = 'diagonal_box_in_box_2d'
#gm = 'straight_channel_2d'
gm = 'diagonal_channels_2d'
#gm = 'perturbed_channel_2d'
#gm = 'vertical_filaments_2d'
#gm = 'circle_in_a_box_2d'
#gm = 'four_boxes_in_a_box_2d'
#gm = 'crossing_channels_2d'

#gm = 'extended_box_3d'
#gm = 'crossing_filaments_3d'
#gm = 'box_in_a_box_3d'
#gm = 'actin_cytosceleton_3d'

if gm == 'box_in__box_2d':
    import GeometryBoxInABox2D as Geometry
elif gm == 'four_boxes_in_a_box_2d':
    import GeometryFourBoxesInABox2D as Geometry
elif gm == 'straight_channel_2d':
    import GeometryStraightChannel2D as Geometry
elif gm == 'perturbed_channel_2d':
    import GeometryPerturbedChannel2D as Geometry
elif gm == 'vertical_filaments_2d':
    import GeometryVerticalFilaments2D as Geometry
elif gm == 'diagonal_channels_2d':
    import GeometryDiagonalChannels2D as Geometry
elif gm == 'crossing_channels_2d':
    import GeometryCrossingChannels2D as Geometry
elif gm == 'diagonal_box_in_box_2d':
    import GeometryDiagonalBoxInABox2D as Geometry
elif gm == 'circle_in_a_box_2d':
    import GeometryCircleInABox2D as Geometry

# 3D meshes
elif gm == 'extended_box_3d':
    import GeometryExtendedBox3D as Geometry
elif gm == 'crossing_filaments_3d':
    import GeometryCrossingFilaments3D as Geometry
elif gm == 'box_in_a_box_3d':
    import GeometryBoxInABox3D as Geometry
elif gm == 'actin_cytosceleton_3d':
    import GeometryActinCytosceleton as Geometry

if __name__ == '__main__':
    
    common.cpp.log.set_log_level(1) 
    
    eps_l_list = [0.01, 0.1, 0.5, 2., 10., 100.]
    #eps_s_list = [1., 3., 9., 27., 81.]
    #eps_l = 2.
    eps_s = 1.
    #res_x = [20 * i for i in range(1,21) ]
    #res_list = [[9*i, 22*i] for i in range(1,21)]
    #res = 360
    #experiment = 'scale_independence_permittivity'
    #experiment = 'accuracy_vs_resolution'
    experiment = 'pourosity_vs_correction_tensors'
    
    all_results = []
    
    #for res in res_x:
    #for box_size in [1 * i for i in range(1,18)]:
    for eps_l in eps_l_list:
        #for eps_s in eps_s_list:
        #print('######   EPS_L EPS_S: ', eps_l, eps_s)
        for i in range(1,16):
            geom = Geometry.Geometry(eps_l, eps_s, resolution=32, ds = 32., sc=i*2.)
            
            pourosity = geom.get_liquid_mesh_volume()/geom.get_mesh_volume()

            xi_phi = hm.reference_cell_solver_potential(geom)
            
            xi_c = hm.reference_cell_solver_diffusion(geom)
            
            diffusion_tensor = hm.diffusion_tensor(xi_c, geom)
            
            permittivity_tensor = hm.permittivity_tensor(xi_phi, geom)
            
            
            
            now = datetime.datetime.now()
            date = now.isoformat()

            pourosity = geom.get_liquid_mesh_volume()/geom.get_mesh_volume()

            results_list = [gm, geom.eps_l, geom.eps_s, geom.res, geom.res_y, geom.res_z, diffusion_tensor, permittivity_tensor, pourosity, geom.domain_size, geom.solid_scale, date, experiment]
            
            all_results.append(results_list)

    pickle.dump(all_results, file=open('./../results/homogenization_' + gm + '.pcl', 'wb')) 


