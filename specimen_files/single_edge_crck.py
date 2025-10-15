# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 22:08:32 2025

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

filename = os.path.join(input_fileDirectory, 'single_edge_crck.inp')
file_save_directory = os.path.join(root_directory,'results')

probname = 'single edge crack'

el_data, nd_data, elem_count, node_count = get_gmsh_data(filename, dim=2)
plot_mesh(el_data, nd_data)
# gpt_crds
gpts_crd = get_gpts_crds(node_data, elem_data)
# make gpt_df 
gpts_df = DataFrame(gpts_crd)

#%%% material props
nu = 0.3
E = 2.1*1e5

#%%% dof setting
# fix dofs, load_dofs, pres_dofs, free_dofs
nd_fix_u1 = array([1, 2, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21], dtype=int)
nd_fix_u2 = array([1, 2, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21], dtype=int)

nd_load_u1 = array([], dtype=int)
nd_load_u2 = array([3, 4, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35], dtype=int)

fill_crack = array([5, 6, 7, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49], dtype=int)
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