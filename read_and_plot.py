# %% [ IMPORTS ]
from specimen_files.unibar_redXs import fill_data, umax, bc_data, node_data
from specimen_files.unibar_redXs import file_save_directory, probname
from Plotting_file import *

# %% pickle reading from save_file_location
# for automated plotting get imports
# =============================================================================
# # manual plotting
# fsd = r'H:\kp\Phase_field_code\results'
# prob_name = 'unibar axial'
# =============================================================================

ref_loc = r'E:\Codes_GitSync\validation_data\unibar\redXS\actual_redXS.csv'
file_path = fr'{file_save_directory}\{probname}\data.pkl'

with open(file_path, 'rb') as f:
    loaded_variables = pickle.load(f)
    f.close()
globals().update(loaded_variables)


#%% PLOTTING
stp = -1

# mesh
# plot_mesh(el_data, node_data)

# displacement contour
plot_contour(plot_this=u1[f'{u1.columns[stp]}'], plot_on='nodes', crds=node_data, text_TT='U1 plot',
             stpNum=int(u1.columns[stp][-1]), fill_data=fill_data)


# damage contour


# Structural response / TSM
# plot_structural_response(srmat)
# plot_time_statistics(tsmat)

