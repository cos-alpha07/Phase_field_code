# %% [ IMPORTS ]
from specimen_files.unibar_redXs import fill_data, umax, bc_data, node_data
from Plotting_file import *

# %% pickle reading from save_file_location
# for automated plotting get imports
# manual plotting

# fsd = r'H:/kp/Phase_field_code/results/single_edge_crck'
fsd = r'H:/kp/Phase_field_code/results/unibar_tension'
prob_name = '24-11-2025__lc=5'

ref_loc = r'H:/kp/Phase_field_code/validationData/unibar_tension/2025_Pandey_adaptive_PFM.csv'
# ref_loc = r'H:/kp/Phase_field_code/validationData/single_edge_crck/2025_edge_crck_plate_lc=0.03.csv'
file_path = fr'{fsd}\{prob_name}\data.pkl'

with open(file_path, 'rb') as f:
    loaded_variables = pickle.load(f)
    f.close()
globals().update(loaded_variables)

#%% PLOTTING
stp = -1

# mesh
# plot_mesh(el_data, nd_data)

# displacement contour
plot_contour(plot_this=u1[f'{u1.columns[stp]}'], plot_on='nodes', crds=node_data, text_TT='U1 plot',
             stpNum=int(u1.columns[stp][5::]), fill_data=fill_data)

plot_contour(plot_this=u2[f'{u2.columns[stp]}'], plot_on='nodes', crds=node_data, text_TT='U2 plot',
             stpNum=int(u2.columns[stp][5::]), fill_data=fill_data)

# damage contour
plot_contour(plot_this=phi[f'{phi.columns[stp]}'], plot_on='nodes', crds=node_data, text_TT='Damage plot',
             stpNum=int(phi.columns[stp][5::]), fill_data=fill_data)

# Structural response / TSM
plot_structural_response(srmat, validate=True, ref_loc=ref_loc)
#plot_time_statistics(tsmat)
