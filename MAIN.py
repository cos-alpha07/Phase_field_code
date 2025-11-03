# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 22:40:53 2025

@author: KP
"""

#%% imports
from specimen_files.unibar_redXs import *

#%% material props
D = define_constitutive_stiffness_matrix(E, nu, 'PLANE_STRAIN')

#%% [INITIALIZING] nodal state variables
Ug = zeros(2*len(node_data))
dU = zeros(2*len(node_data))
Phig = zeros(len(node_data))
dPhi = zeros(len(node_data))
hist_param = zeros(len(gpts_crd))

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
    dPhi[:] = 0.0
    app_load += u_step
    
    stp += 1
    print(f'Running load-step - {stp}')
    t1 = time.process_time()
    
    tempu = Ug.copy()
    tempp = Phig.copy()
    
    err_ = 1
    itr = 0
    
    # loop for iterations
    while err_ > tol:
        itr += 1
        print(f'({stp}) load={app_load}\t Iteration - {itr}', end=', ')
        
        # build Kmat, (Fint, Fext) -> Resd
        K_mat, K_phi, F_int, F_phi, hist_param = assemble_forces_and_stiffness(elem_data, node_data, tempu, tempp, D,
                                                                   hist_param, E, f_t, G_c, l_c, thk)
        #K_disp, K_phi, F_disp, F_phi = assemble_KK_Fi(elem_data, node_data, tempu, tempp, D)
        
        # [DISP]
        
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
        
        # [PHI]
        dPhi = solve(K_phi, -F_phi)
        tempp += dPhi
        
        err_disp = norm(dU)/(norm(tempu) + 1e-8)
        err_phi = norm(dPhi)/(norm(tempp) + 1e-8)
        err_ = max(err_disp, err_phi)
        
        # convergence check
        if check_convergence(err_, tol):
            step_time = time.process_time() - t1
            Ug = tempu.copy()
            Phig = tempp.copy()
            
            # updating plot-data
            SR_mat = append(SR_mat, [[app_load, -sum(F_int[fix_dofs])]], axis=0)
            TS_mat = append(TS_mat, [[stp, itr, step_time]], axis=0)
            
            print(f'  **> Solution converged at iteration-{itr}')
            break
        else:
            print(f'err_disp={err_disp:5e}, err_phi={err_phi:5e}, error={err_:5e}')
    
    # creating dataframes
    u1_df[f'Step-{stp}'] = Series(Ug[[2*i for i in range(node_count)]], index=range(node_count))
    u2_df[f'Step-{stp}'] = Series(Ug[[2*i + 1 for i in range(node_count)]], index=range(node_count))
    phi_df[f'Step-{stp}'] = Series(Phig)
    save_my_results(u1_df, u2_df, phi_df, zeros(2), SR_mat, TS_mat, elem_data, node_data, zeros([2, 4]),
                    save_loc=fr'{file_save_directory}\{probname}')

sys.exit()

#%% Intermittent Plotting
plot_contour(plot_this=Ug[[2*i for i in range(node_count)]], plot_on='nodes', crds=node_data, text_TT='U1',
             stpNum=stp, fill_data=fill_data)

    
    
    


