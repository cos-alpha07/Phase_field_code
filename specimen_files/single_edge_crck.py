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

filename = os.path.join(input_fileDirectory, 'single_edge_crck_coarse.inp')
file_save_directory = os.path.join(root_directory,'results', 'single_edge_crck')
fixed_date = datetime.today().strftime('%d-%m-%Y')

elem_data, node_data, elem_count, node_count = get_gmsh_data(filename, dim=2)
# gpt_crds
gpts_crd = init_gpt_mat(elem_data)

#%%% material props
nu = 0.3
E = 2.1*1e5
G_c = 2.7
f_t = 2445.42
l_c = 0.03
thk = 1

#%%% dof setting
# fix dofs, load_dofs, pres_dofs, free_dofs
nd_fix_u1 = array([1, 2, 8, 9, 10, 11, 12, 13, 14, 15, 
16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 
26], dtype=int)
nd_fix_u2 = array([1, 2, 8, 9, 10, 11, 12, 13, 14, 15, 
16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 
26], dtype=int)

nd_load_u1 = array([], dtype=int)
nd_load_u2 = array([3, 4, 74, 75, 76, 77, 78, 79, 80, 81, 
82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 
92], dtype=int)

fill_crack = array([5, 6, 7, 102, 103, 104, 105, 106, 107, 108, 
109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 
119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 
129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 
139, 140, 141, 142, 143, 144, 145, 146, 147], dtype=int)

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
ustep = 5*1e-4
num_steps = int(umax/ustep)
tol = 1e-6

probname = fixed_date + f'__lc={l_c}'