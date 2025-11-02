# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 13:21:00 2025

@author: KP
"""

# %% [IMPORTS]
from impLibrary import *

# %% [FE functions]
def get_conn(elem_no, conn):
    res, = [row[1:] for row in conn if row[0] == elem_no]
    return res

def get_el_num(nodal_con, elcon):
    res, = [row[0] for row in elcon if set(row[1:]) == set(nodal_con)]
    return res

def get_node_coord(node_num, node_info):
    res, = [row[1:] for row in node_info if row[0] == node_num]
    return res

def get_node_no(node_coord, node_info):
    res, = [row[0] for row in node_info if set(row[1:]) == set(node_coord)]
    return res

def get_el_type(n_nodes):
    if n_nodes == 3:
        int_pts = 1
        loc =  array([[1/3, 1/3]])
        wts =  array([0.5])
    elif n_nodes == 4:
        int_pts = 4
        a = 1/sqrt(3)
        loc =  array([[-a, -a],
                        [a, -a],
                        [a, a],
                        [-a, a]])
        wts =  array([1, 1, 1, 1])
    else:
        raise ValueError("Unsupported element node count")   
        
    return int_pts, loc, wts

def get_shape_functions(n_nodes, xi, eta):
    if n_nodes == 3:
        N =  array([1.0 - xi - eta, xi, eta])
        dN =  array([[-1, -1],
              [1, 0],
              [0, 1]])   
    elif n_nodes == 4:
        N = 0.25* array([(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)])
        dN = 0.25* array([[-(1-eta), -(1-xi)],
             [(1-eta), -(1+xi)],
             [(1+eta), (1+xi)],
             [-(1+eta), (1-xi)]])
    else:
        raise ValueError("Unsupported element node count")      
    return N, dN
         
def comp_B_global(dN, n_nodes, coords):
    dN_dxi_deta = dN.T
    J = dN_dxi_deta @ coords
    det_J =   det(J)
    if det_J <= 0:    
        raise ValueError("Non-positive Jacobian determinant: detJ = %g" % det_J)
    inv_J  =   inv(J)
    B_phi = inv_J @ dN_dxi_deta     #2x4
    B =  zeros((3, 2*n_nodes))
    for i in range(n_nodes):
        dN_dx = B_phi[0,i]
        dN_dy = B_phi[1,i]
        B[0, 2*i] = dN_dx
        B[1, 2*i+1] = dN_dy
        B[2, 2*i] = dN_dy
        B[2, 2*i+1] = dN_dx 
    return B, B_phi, det_J

# def comp_B_phi(dN, n_nodes, coords):
#     dN_dxi_deta = dN.T 
#     J = dN_dxi_deta @ coords
#     inv_J  = inv(J)
#     dN_dx_dy = inv_J @ dN_dxi_deta
#     B_phi = zeros((2, n_nodes))
#     for i in range(n_nodes):
#         dN_dx = dN_dx_dy[]
        
def comp_stiffness_mat_and_forces(idx, D_mat, n_nodes, u_e, p_e, coords, hist, emod, f_t, G_c, l_c, thk):
    int_pts, loc, wts = get_el_type(n_nodes)
    K_e =  zeros((2*n_nodes, 2*n_nodes))
    K_p = zeros((n_nodes, n_nodes))
    f_int_e =  zeros(2*n_nodes)
    f_p = zeros(n_nodes)
    c_o = pi
    for i in range(int_pts):
        xi, eta  = loc[i, :]
        N, dN = get_shape_functions(n_nodes, xi, eta)
        # add function to get Bphi
        B, B_phi, det_J = comp_B_global(dN, n_nodes, coords)
        # add function to evaluate strain, stress, energy, phase-field value
        phi = N @ p_e.T
        phi_grad = B_phi @ p_e.T    #2x1
        # add function for evaluating alpha(phi), omega(phi) and there derivatives
        al_phi, al_phi_d, al_phi_dd, om_phi, om_phi_d, om_phi_dd = comp_phase_field_vars(phi, emod, f_t, G_c,
                                                                                         l_c, 'linear')
        strain = B @ u_e                  #
        stres = D_mat @ strain           #
        sigma_1 = eigvals(array([[stres[0], stres[2]], [stres[2], stres[1]]], dtype=float))[0]
        
        psi_0 = 0.5 * comp_macaulay(sigma_1, 'positive')**2/emod 
        
        hist_gpt = max(hist[idx], psi_0, 0.5*(f_t**2)/emod)
        
        # update history parameter in the global storage vector also
        
        # assembly for element
        dvol = det_J*wts[i]*thk
        K_e +=  om_phi*(B.T @ D_mat @ B)*dvol # update this
        
        K_p += ((al_phi_dd * G_c/(l_c*c_o)) + om_phi_dd * hist_gpt)*N[:, newaxis] @ N[newaxis, :]*dvol 
        K_p += 2*G_c*l_c/c_o * B_phi.T @ B_phi*dvol
        
        iK_p = inv(K_p)
        
        f_int_e += (B.T @ stres)*dvol
        
        f_p += (om_phi_d * hist_gpt + al_phi_d*G_c/(l_c*c_o))*N*dvol
        f_p += 2*G_c*l_c/c_o * B_phi.T @ phi_grad*dvol
        
        idx += 1
        
    return K_e, K_p, f_int_e, f_p, idx

def assemble_forces_and_stiffness(el_conn, node_info, a_global, p_global, D_mat, hist, emod, f_t, G_c, l_c, thk):
    max_node = int(amax(el_conn))
    n_dof = 2 * int(max_node)
    f_int = zeros(n_dof)
    f_phi = zeros(max_node)
    K = zeros((n_dof, n_dof))
    K_phi = zeros((max_node, max_node))
    idx = 0
    for row in el_conn:
        conn = get_conn(row[0], el_conn)
        cl_conn = [nd for nd in conn if nd != 0]
        dof_map = [dof for n in cl_conn for dof in (2*(n-1), 2*(n-1)+1)]
        dof_map_p = [n - 1 for n in cl_conn]
        u_e = a_global[dof_map]
        p_e = p_global[dof_map_p]
        coords =  asarray([get_node_coord(n, node_info) for n in cl_conn], dtype=float)
        n_nodes_elem = len(cl_conn)
        K_e, K_p, f_int_e, f_p, idx = comp_stiffness_mat_and_forces(idx, D_mat, n_nodes_elem,
                                                                    u_e, p_e, coords, hist, emod, 
                                                                    f_t, G_c, l_c, thk)
        for i, gi in enumerate(dof_map):
            f_int[gi] += f_int_e[i]
        for j, gj in enumerate(dof_map_p):
            f_phi[gj] += f_p[j]
            
        K[ix_(dof_map, dof_map)] += K_e
        K_phi[ix_(dof_map_p, dof_map_p)] += K_p
        
    return K, K_phi, f_int, f_phi

def define_constitutive_stiffness_matrix(emod, enu, load_state):
    # load_state - {Plane strain, Plane stress}
    D_mat = zeros((3, 3))
    
    if load_state.upper() in ('PLANE STRAIN', 'PLANE_STRAIN'):
         D_mat[0,0]=emod*(1.0-enu)/((1.0+enu)*(1.0-2.0*enu))
         D_mat[0,1]= D_mat[0,0]*enu/(1.0-enu)
         D_mat[1,0]= D_mat[0,1]
         D_mat[1,1]= D_mat[0,0]
         D_mat[2,2]= D_mat[0,0]*(0.5-enu)/(1.0-enu)
        
    elif load_state.upper() in ('PLANE STRESS', 'PLANE_STRESS'):
         D_mat[0,0]=emod/(1.0-enu*enu)
         D_mat[0,1]= D_mat[0,0]*enu
         D_mat[1,0]= D_mat[0,1]
         D_mat[1,1]= D_mat[0,0]
         D_mat[2,2]=emod/(2.0*(1+enu))
    else:
        raise ValueError('Incorrect load_state defined, exiting code')
        
    return D_mat

def check_convergence(error, tol):
    cnvg = False
    if error < tol:
        cnvg = True
    return cnvg

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

def comp_phase_field_vars(dmg, emod, f_t, G_c, l_c, soft_type):
    l_ch = emod*G_c/f_t**2
    a_1 = 4*l_ch/(pi*l_c)
    
    if soft_type.upper() == 'LINEAR':
        a_2, a_3, p = -0.5, 0, 2
    elif soft_type.upper() == 'CONCRETE':
        a_2, a_3, p = 1.3868, 0.9107, 2
        
    phi = symbols('phi')
    alpha_phi = 2*phi - phi**2
    alpha_phi_d = diff(alpha_phi, phi)
    alpha_phi_dd = diff(alpha_phi_d, phi)
    
    omega_phi = (1-phi)**p/((1-phi)**p + a_1*phi*(1 + a_2*phi + a_2*a_3*phi**2))
    omega_phi_d = diff(omega_phi, phi)
    omega_phi_dd = diff(omega_phi_d, phi)
    
    ret_res = [alpha_phi.evalf(subs={phi:dmg}), alpha_phi_d.evalf(subs={phi:dmg}), 
               alpha_phi_dd.evalf(subs={phi:dmg}), omega_phi.evalf(subs={phi:dmg}), 
               omega_phi_d.evalf(subs={phi:dmg}), omega_phi_dd.evalf(subs={phi:dmg})]

    return [float(_) for _ in ret_res]
    
def comp_macaulay(val, part):
    if part.upper() == 'POSITIVE':
        return (val + abs(val))/2
    elif part.upper() == 'NEGATIVE':
        return (abs(val) - val)/2
    else:
        return (val + abs(val))/2
        
    
    
#%% Utilty functions
def initialize_save_files(fsd, prob_name):
    if not os.path.exists(fr'{fsd}\{prob_name}'):
        os.makedirs(fr'{fsd}\{prob_name}')
    
    return DataFrame(), DataFrame(), DataFrame()


def save_my_results(u1, u2, phi, his_par, srmat, tsmat, el_data, nd_data, gpt_data, save_loc):   
    # dumping data-pickles
    dumps = dict(
        [
            ("u1", locals()["u1"]),
            ("u2", locals()["u2"]),
            ("srmat", locals()["srmat"]),
            ("tsmat", locals()["tsmat"]),
            ("nd_data", locals()["nd_data"]),
            ("el_data", locals()["el_data"]),
            ("gpt_data", locals()["gpt_data"]),
            ("phi", locals()["phi"]),
        ]
    )
    # Define the file path where the selected variables will be saved
    dump_file_path = os.path.join(save_loc, 'data.pkl')
    
    # Save the selected variables to a file
    with open(dump_file_path, "wb") as f:
        pickle.dump(dumps, f)
        
    print("data saved\n")

def get_gmsh_data(file_name, dim=2):
    nd_data = []
    el_data = []
    nd_idxs = []
    el_idxs = []

    with open(file_name, 'r') as file:
        line = file.readlines()
        file.close()

    # reading node data
    for i, l in enumerate(line):
        if l.startswith('*NODE'):
            nd_idxs.append(i+1)
        elif l.startswith('*ELEMENT'):
            nd_idxs.append(i-1)
            break

    for l in line[nd_idxs[0]:nd_idxs[1]]:
        nd_data.append(l.split(', '))
    nd_data = array(nd_data, dtype=float)

    # reading element data
    for i, l in enumerate(line):
        if l.startswith('*ELEMENT'):
            if 'Surface' in l:
                el_idxs.append(i+1)

    for l in line[el_idxs[0]:]:
        if l.startswith('*ELSET') or l.startswith('*NSET'):
            break
        elif l.startswith('*'):
            continue

        if len(l.split(', ')) == 4:
            el_data.append(l.split(', ') + ['0'])
        else:
            el_data.append(l.split(', '))

    el_data = array(el_data, dtype=int)
    el_data[:, 0] = arange(1, len(el_data)+1, 1, dtype=int)

    return el_data, nd_data[:, 0:dim+1], len(el_data), len(nd_data)