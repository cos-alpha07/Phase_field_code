# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 11:35:58 2025

@author: Dell
"""

#%% imports
from file_func import *
from Plotting_file import *
from impLibrary import *

#%% Pre-processing
#%%% reading mesh data
root_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
input_fileDirectory = os.path.join(root_directory, 'inpData')

filename = os.path.abspath(input_fileDirectory, 'L_shaped_panel.inp')
file_save_directory = os.path.join(root_directory,'results')

el_data, nd_data, _, _ = get_gmsh_data(filename, dim=2)
plot_mesh(el_data, nd_data)

probname = 'L shaped panel'

el_data, nd_data, elem_count, node_count = get_gmsh_data(filename, dim=2)
plot_mesh(el_data, nd_data)
# gpt_crds
gpts_crd = get_gpts_crds(node_data, elem_data)
# make gpt_df 
gpts_df = DataFrame(gpts_crd)

#%%% material props
nu = 0.18
E = 2.585*1e4

#%%% dof setting
# fix dofs, load_dofs, pres_dofs, free_dofs
nd_fix_u1 = array([26, 27, 28, 29, 30, 31, 32, 33, 34], dtype=int)
nd_fix_u2 = array([26, 27, 28, 29, 30, 31, 32, 33, 34], dtype=int)

nd_load_u1 = array([], dtype=int)
nd_load_u2 = array([7], dtype=int)

fill_crack = array([], dtype=int)
fill_set1 = array([], dtype=int)
bc_data = [nd_fix_u1, nd_fix_u2, nd_load_u1, nd_load_u2]

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

#%%% solver criteria 
    # - loading_type - load ctrl, disp-ctrl
    # - load_max, load_step, tolerance
    
load_type = 'disp'    
umax = 0.006
num_steps = 60
u_step = 1e-4
tol = 1e-6
