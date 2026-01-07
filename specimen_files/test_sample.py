# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 12:54:26 2025

@author: Dell
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 16:51:54 2025

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

filename = os.path.join(input_fileDirectory, 'test_sample.inp')
file_save_directory = os.path.join(root_directory,'results', 'test_sample')
fixed_date = datetime.today().strftime('%d-%m-%Y')


elem_data, node_data, elem_count, node_count = get_gmsh_data(filename, dim=2)
# gpt_crds
gpts_crd = get_gpts_crds(node_data, elem_data)
# make gpt_df 
gpts_df = DataFrame(gpts_crd)

#%%% material props
nu = 0.2
E = 3e4       #MPa
G_c = 0.008    #N/mm
f_t = 3.0        #MPa
l_c = 5         # mm
thk = 1

#%%% dof setting
# fix dofs, load_dofs, pres_dofs, free_dofs
nd_fix_u1 = array([1], dtype=int)
nd_fix_u2 = array([2], dtype=int)

nd_load_u1 = array([], dtype=int)
nd_load_u2 = array([5], dtype=int)

fill_crack = array([], dtype=int)
fill_set1 = array([], dtype=int)
fill_set2 = array([], dtype=int)

bc_data = [nd_fix_u1, nd_fix_u2, nd_load_u1, nd_load_u2]

fill_data = [fill_crack, fill_set1, fill_set2]

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
umax = 0.02
ustep = 0.0002
num_steps = int(umax/ustep)
tol = 1e-5

probname = fixed_date + f'__lc={l_c}'