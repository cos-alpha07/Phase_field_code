# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 22:40:53 2025

@author: KP
"""

#%% imports
from specimen_files.unibar_redXs import *

#%% material props
D = define_constitutive_stiffness_matrix(E, nu, 'PLANE_STRESS')

#%% [INITIALIZING] nodal state variables
Ug = zeros(2*len(node_data))
dU = zeros(2*len(node_data))
PHIg = zeros(len(node_data))
dPHI = zeros(len(node_data))

u1_df, u2_df, phi_df = initialize_save_files(file_save_directory, probname)

#%% mesh plot
mesh_display = False
if mesh_display:
    mdt = time.process_time()
    plot_mesh(elem_data, node_data)
    print(f'Mesh Plotted, time - {time.process_time() - mdt}')
    
#%% Solver
app_load = 0
stp = 0
SR_mat = zeros([1, 2]) 
TS_mat = empty([0, 3])

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
        K_mat, F_int = assemble_KK_Fi(elem_data, node_data, tempu, D)
        #K_disp, K_phi, F_disp, F_phi = assemble_KK_Fi(elem_data, node_data, tempu, tempd, D)
                
        # applying boundary conditions
        K_ff = K_mat[ix_(fdofs, fdofs)]
        K_fp = K_mat[ix_(fdofs, pdofs)]
        
        dU[load_dofs] = app_load
        
        if itr == 1:
            rhs = -(K_fp @ dU[pdofs] + F_int[fdofs])
            
        else:
            rhs = -F_int[fdofs]
            dU[pdofs] = 0 
           
        dU[fdofs] = solve(K_ff, rhs)
        
        # adding dU to tempu += dU
        tempu[load_dofs] = 0
        tempu += dU
        err_ = norm(dU)/(norm(tempu) + 1e-8)
        
        # convergence check
        if check_convergence(err_, tol):
            step_time = time.process_time() - t1
            Ug = tempu.copy()
            
            # updating plot-data
            SR_mat = append(SR_mat, [[app_load, -sum(F_int[fix_dofs])]], axis=0)
            TS_mat = append(TS_mat, [[stp, itr, step_time]], axis=0)
            
            print(f'  **> Solution converged at iteration-{itr}')
            break
        else:
            print(f'error={err_:5e}')
    
    # creating dataframes
    u1_df[f'Step-{stp}'] = Series(Ug[[2*i for i in range(node_count)]], index=range(node_count))
    u2_df[f'Step-{stp}'] = Series(Ug[[2*i + 1 for i in range(node_count)]], index=range(node_count))
    phi_df[f'Step-{stp}'] = Series(PHIg)
    save_my_results(u1_df, u2_df, phi_df, zeros(2), SR_mat, TS_mat, elem_data, node_data, zeros([2, 4]), save_loc=fr'{file_save_directory}\{probname}')
    
    
    


