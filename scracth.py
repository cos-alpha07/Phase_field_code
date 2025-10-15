# from specFiles.MODE1_2D import fill_set_data, load_max, Dim, bc_data
# from Plotting_file import plot_contour, plot_mesh
# import pickle, os

# # %% pickle reading from save_file_location
# fsd = r'E:\Codes_GitSync\results'
# prob_name = '26-09-2025__quasi-static, standard_UALCNP_VICH__UNIBAR_redXS_2D__lc=10__smth=1'

# ref_loc = r'E:\Codes_GitSync\validation_data\unibar\redXS\actual_redXS.csv'
# file_path1 = fr'{fsd}\{prob_name}\Mesh\Step-data.pkl'
# file_path2 = fr'{fsd}\{prob_name}\data.pkl'

# if os.path.exists(fr'{fsd}\prob_name\Mesh'):
#     with open(file_path1, 'rb') as f:
#         loaded_variables1 = pickle.load(f)
#         f.close()
#     globals().update(loaded_variables1)

# with open(file_path2, 'rb') as f:
#     loaded_variables2 = pickle.load(f)
#     f.close()
# globals().update(loaded_variables2)


# #%% PLOTTING
# # mesh
# # show_mesh(elem_df.values, node_df.values, Dim, show_bcs=True, bcs=bc_data)

# # displacement contour
# plot_contour(plot_this=u1[u1.columns[-1]].values, plot_on='nodes', plot_what='U1', dim=Dim, el_data=elem_df.values,
#              nd_data=node_df.values)
# plot_contour(plot_this=u2[u2.columns[-1]].values, plot_on='nodes', plot_what='U2', dim=Dim, el_data=elem_df.values,
#              nd_data=node_df.values)

# # damage contour
# plot_contour(plot_this=dmg[dmg.columns[-1]].values, plot_on='nodes', plot_what='PHI', dim=Dim, el_data=elem_df.values,
#              nd_data=node_df.values)

# # Structural response / TSM
# # kwargs for saving plots -> save_loc=fr'{fsd}/{prob_name}', sname

# plot(SR=SR_df, xunit='mm', yunit='N', compare_sr=False, ref_loc=ref_loc, annotate_data=True)

# plot(TSM=TSMat_df)
# plot(SH=stepHistory)

from numpy import asarray, array, append
from file_func import *



elem_mat = array([[1, 1, 2, 3, 4],
                  [2, 4, 3, 5, 0]], dtype=int)

node_mat = array([[1, 0, 0],
                  [2, 1, 0],
                  [3, 1, 1],
                  [4, 0, 1],
                  [5, 0.5, 2]], dtype=float)


def get_gpts_crds(nd_mat, el_conn):
    gpt_crds = empty([0, 2])
    for row in el_conn:
        cl_conn = [nd for nd in get_conn(row[0], el_conn) if nd != 0]
        coords =  asarray([get_node_coord(n, nd_mat) for n in cl_conn], dtype=float)
        n_nodes_elem = len(cl_conn)
        int_pts, loc, wts = get_el_type(n_nodes_elem)
        for i in range(int_pts):
            xi, eta  = loc[i, :]
            N, dN = get_shape_functions(n_nodes_elem, xi, eta)
            gpt_crds = append(gpt_crds, [N @ coords], axis=0)
    return gpt_crds
        
gpts = get_gpts_crds(node_mat, elem_mat)