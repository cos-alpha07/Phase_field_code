# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 13:21:00 2025

@author: KP
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 17:03:19 2025

@author: Dell
"""
#%% [IMPORTS]
from math import sqrt
from numpy import array, zeros, ix_, dot, asarray, empty, max
from scipy.linalg import det, inv, norm
import time
import pickle


#%% [FE functions]
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
                        [-a, a],
                        [a, -a],
                        [a, a]])
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
    dN_dx_dy = inv_J @ dN_dxi_deta
    B =  zeros((3, 2*n_nodes))
    for i in range(n_nodes):
        dN_dx = dN_dx_dy[0,i]
        dN_dy = dN_dx_dy[1,i]
        B[0, 2*i] = dN_dx
        B[1, 2*i+1] = dN_dy
        B[2, 2*i] = dN_dy
        B[2, 2*i+1] = dN_dx 
    return B, det_J

def comp_stiffness_mat(D_mat, n_nodes, coords):
    int_pts, loc, wts = get_el_type(n_nodes)
    K_e =  zeros((2*n_nodes, 2*n_nodes))
    for i in range(int_pts):
        xi, eta  = loc[i, :]
        N, dN = get_shape_functions(n_nodes, xi, eta)
        B, det_J = comp_B_global(dN, n_nodes, coords)
        
        # assembly for element
        dvol = det_J*wts[i]
        K_e +=(B.T @ D_mat @ B)*dvol
        # Kglobal[ix_(dof_map, dof_map)] += (B.T @  D_mat @ B)*dvol
        # return Kglobal
    return K_e

def assemble_stiffnes_mat(el_conn, node_info, D_mat):
    max_node = int(max(el_conn))
    n_dof = 2 * int(max_node)
    K  =  zeros((n_dof, n_dof))
    for row in el_conn:
        cl_conn = [nd for nd in get_conn(row[0], el_conn) if nd != 0]
        coords =  asarray([get_node_coord(n, node_info) for n in cl_conn], dtype=float)
        n_nodes_elem = len(cl_conn)
        K_e = comp_stiffness_mat(D_mat, n_nodes_elem, coords)
        dof_map = []
        for node in cl_conn:
            dof_u = 2*(node-1)
            dof_v = 2*node-1
            dof_map.append(dof_u)
            dof_map.append(dof_v)
            
            K[ix_(dof_map), ix_(dof_map)] += K_e
    return K

def comp_fint(n_nodes, coords, u_e, D_mat):
    int_pts, loc, wts = get_el_type(n_nodes)
    f_int_e =  zeros(2*n_nodes)
    for i in range(int_pts):
        xi, eta = loc[i,0], loc[i,1]
        _, dN = get_shape_functions(n_nodes, xi, eta)
        B, det_J = comp_B_global(dN, n_nodes, coords)
        stres = D_mat @ (B @ u_e)
        f_int_e += (B.T @ stres)*det_J*wts[i]
    return f_int_e

def assemble_fint(el_conn, node_info, a_global, D_mat):
    max_node = int(max(el_conn))
    n_dof = 2 * int(max_node)
    f_int =  zeros(n_dof)
    for row in el_conn:
        conn = get_conn(row[0], el_conn)
        cl_conn = [nd for nd in conn if nd != 0]
        dof_map = [dof for n in cl_conn for dof in (2*(n-1), 2*(n-1)+1)]
        u_e = a_global[dof_map]
        coords =  asarray([get_node_coord(n, node_info) for n in cl_conn], dtype=float)
        n_nodes_elem = len(cl_conn)
        f_int_e = comp_fint(n_nodes_elem, coords, u_e, D_mat)
        for i, gi in enumerate(dof_map):
            f_int[gi] += f_int_e[i]         
    return f_int


def assemble_KK_Fi(el_data, nd_data, tempu, D_mat):
    K_mat = assemble_stiffnes_mat(el_data, nd_data, D_mat)
    F_int = assemble_fint(el_data, nd_data, tempu, D_mat)
    return K_mat, F_int


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
