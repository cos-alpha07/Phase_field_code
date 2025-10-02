# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 22:40:53 2025

@author: KP
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 15:55:09 2025

@author: Dell
"""

#%% imports
from file_func import *
from Plotting_file import *
from numpy import array, zeros, ix_, dot, asarray, empty, sum, append
from scipy.linalg import det, inv, norm
import time


#%% Pre-processing
#%%% reading mesh data
elem_mat = array([[1, 1, 2, 3, 4],
                  [2, 4, 3, 5, 0]], dtype=int)

node_mat = array([[1, 0, 0],
                  [2, 1, 0],
                  [3, 1, 1],
                  [4, 0, 1],
                  [5, 0.5, 2]], dtype=float)


elem_count, node_count = len(elem_mat), len(node_mat)

#%%% material props
nu = 0.3
E = 1e5
D = define_constitutive_stiffness_matrix(E, nu, 'PLANE_STRESS')


#%%% dof setting
# fix dofs, load_dofs, pres_dofs, free_dofs
nd_fix_u1 = array([1], dtype=int)
nd_fix_u2 = array([1, 2], dtype=int)

nd_load_u1 = array([], dtype=int)
nd_load_u2 = array([5], dtype=int)

bc_data = [nd_fix_u1, nd_fix_u2, nd_load_u1, nd_load_u2]

fill_crack = array([], dtype=int)
fill_set1 = array([], dtype=int)

fill_data = [fill_crack, fill_set1]

xfix_dof = [2*i-2 for i in nd_fix_u1]
yfix_dof = [2*i-1 for i in nd_fix_u2]
fix_dofs = xfix_dof + yfix_dof

xload_dof = [2*i-2 for i in nd_load_u1]
yload_dof = [2*i-1 for i in nd_load_u2]
load_dofs = xload_dof + yload_dof

total_disp_dofs = [i for i in range(2*node_count)]
pdofs = fix_dofs + load_dofs
fdofs = [i for i in total_disp_dofs if i not in pdofs]

Ug = zeros(2*len(node_mat))
dU = zeros(2*len(node_mat))


#%% mesh plot
mesh_display = False
if mesh_display:
    mdt = time.process_time()
    plot_mesh(elem_mat, node_mat)
    print(f'Mesh Plotted, time - {time.process_time() - mdt}')
    

#%% Solver
app_load = 0
stp = 0
SR_mat = zeros([1, 2]) 
TS_mat = empty([0, 3])

#%%% solver criteria 
    # - loading_type - load ctrl, disp-ctrl
    # - load_max, load_step, tolerance
    
load_type = 'disp'
umax = 1
num_steps = 100
u_step = 0.01
tol = 1e-4

#%%% NR begins
while app_load < umax:    # loop on steps
    dU[:] = 0.0
    app_load += u_step
    stp += 1
    print(f'Running load-step - {stp}')
    t1 = time.process_time()
    
    tempu = Ug.copy()

    err_ = 1
    itr = 0
    # loop for iterations
    while err_ > tol:
        itr += 1
        print(f'({stp} (load={app_load}))\t Iteration - {itr}', end=', ')
        
        # build Kmat, (Fint, Fext) -> Resd
        K_mat, F_int = assemble_KK_Fi(elem_mat, node_mat, tempu, D)
        
        # applying boundary conditions
        K_ff = K_mat[ix_(fdofs, fdofs)]
        K_fp = K_mat[ix_(fdofs, pdofs)]
        
        if itr == 1:
            dU[fdofs] = dot(-inv(K_ff), (dot(K_fp, dU[pdofs]) + F_int[fdofs]))
            dU[yload_dof] = app_load
            dU[fix_dofs] = 0
        else:
            dU[fdofs] = dot(-inv(K_ff), F_int[fdofs])
            dU[yload_dof] = app_load
            dU[fix_dofs] = 0

        # adding dU to tempu += dU
        tempu += dU
        err_ = norm(dU)/(norm(tempu) + 1e-8)
        
        # convergence check
        if check_convergence(err_, tol):
            step_time = time.process_time() - t1
            Ug = tempu.copy()
            
            # updating plot-data
            SR_mat = append(SR_mat, [[app_load, -sum(F_int[fix_dofs])]], axis=0)
            TS_mat = append(TS_mat, [[stp, itr, step_time]], axis=0)
            
            print(f'\n **> Solution converged at iteration-{itr}')
            break
        else:
            print(f'error={err_:5e}')
    
    # save_my_results():
        
    
    
    
    
    
    
    
    
    
    
    #plot_structural_response(SR_mat)
    #plot_time_statistics(TS_mat)

