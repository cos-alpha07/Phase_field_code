# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 12:25:12 2025

@author: Dell
"""

#%% [IMPORTS]
from file_func import *
from impLibrary import *
from matplotlib.patches import Polygon

#%% [Contour and Mesh Plots]
def plot_contour(**kwargs):
    # use kwargs for user input for 
    plot_this = kwargs.get('plot_this', None)   # vector value
    plot_on = kwargs.get('plot_on', None)   # character value
    crds = kwargs.get('crds', None)         # tensor value
    text_TT = kwargs.get('text_TT', None)     # character value
    stpNum = kwargs.get('stpNum', None)     # Integer value
    fill_data = kwargs.get('fill_data', None) # tensor data
    
    
    if plot_this is None:
        raise ValueError('No value to plot')
    else:
        
        # plot variable
        if plot_on is not None and plot_on.upper() in ('NODES', 'NODE'):
            print('Nodal plot')
        elif plot_on is not None and plot_on.upper() in ('GPTS', 'GPT', 'GAUSSPOINT'):
            print('Gaussian plot')
        
        # CYX system
        if crds is None:
            raise ValueError('coordinate system undefined')
        else:
            t = tri.Triangulation(crds[:,1], crds[:,2])
        
        # plot begins here
        fig, ax = plt.subplots()
        ax.axis('off')
    
        cntr = ax.tricontourf(t, plot_this, cmap='jet', levels=500,
                              vmin=plot_this.min(), vmax=plot_this.max())
        
        cbar = fig.colorbar(cntr, ax=ax, ticks=linspace(plot_this.min(), plot_this.max(), 5),
                            location='bottom', fraction=0.046, pad=0.06)
        cbar.ax.tick_params(direction='out', labelsize='5', width=0.2) # changes size of ticks only
        for spine in cbar.ax.spines.values():       # changes size of cbar edge width
            spine.set_linewidth(0.3)
            
        # changing the tick labels in scientific formatting
        cbar.ax.ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
        
        # Fill sets
        if fill_data is not None:
            crck_set = fill_data[0]
            if len(crck_set) != 0:
                crck_x = [get_node_coord(i, crds)[0] for i in crck_set]
                crck_y = [get_node_coord(i, crds)[1] for i in crck_set]
            
                ax.fill(crck_x, crck_y, fill=True, color='white', edgecolor='white', lw=0.1)
                
            for _ in range(len(fill_data) - 1):
                if len(fill_data[1+_]) != 0:
                    xcrd = [get_node_coord(i, crds)[0] for i in fill_data[1+_]]
                    ycrd = [get_node_coord(i, crds)[1] for i in fill_data[1+_]]
            
                    ax.fill(xcrd, ycrd, fill=True, color='white', edgecolor='white', lw=0.1, ls='-')
        
        # fill sets (EXTRA)
        # Lpanel
        # ax.fill([250, 250, 500, 500], [0, 250, 250, 0],fill=True, color='white', edgecolor='white', lw=2)
        
        plt.title(f'{text_TT} - Step-{stpNum}')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()
    
def plot_mesh(el_data, nd_data, **kwargs):
    show_nodes = kwargs.get('show_nodes', None)
    show_elems = kwargs.get('show_elems', None)
    show_crack = kwargs.get('show_crack', None)
    show_bcs = kwargs.get('show_bcs', None)
            
    list = []
    for row in el_data:
        conn = get_conn(row[0], el_data)
        cl_conn = [nd for nd in conn if nd != 0]        
        coords =  [get_node_coord(n, nd_data) for n in cl_conn]
        res =  [tuple(row) for row in coords]
        list.append(res)
    
    poly_collection = PolyCollection(list, facecolors='white', lw=0.1, edgecolors='blue', alpha=1)
    fig, ax = plt.subplots()
    ax.add_collection(poly_collection)

    # Set limits 
    ax.axis('off')
    ax.autoscale()
    ax.set_title(f'No. of Elements {len(el_data)}, No. of Nodes {len(nd_data)}')    
    
    if show_nodes:
        # display node no.
        for row in nd_data:
            x, y = row[1:]
            plt.text(x, y, str(int(row[0])), size=1)
    
    if show_elems:
        # show elem no.
        for elem in el_data:
            i, resx, resy = 0, 0, 0
            for nd in elem[1:]:
                if nd != 0:
                    x, y  = get_node_coord(nd, nd_data)
                    resx += x
                    resy += y
                    i += 1
            resx, resy = resx/i, resy/i
            plt.text(resx, resy, str(int(elem[0])), size=1)
    
    if show_crack:
        # display crack by using ax.fill again, just fill red color to display crack
        pass
    
    if show_bcs:
        
        bc_nodes = kwargs.get('bc_nodes', None)
        if bc_nodes is None:
            raise ValueError('BCs nodes (bc_nodes=) not specified in function call')
        else:
            # color codes for showing bcs
            # x fix -- square marker black, face color white
            for nd in bc_nodes[0]:
                x, y = get_node_coord(nd, nd_data)
                ax.plot(x, y, ls='', marker='>', ms=4, c='k')
            # y fix -- diamond marker back, face color white
            for nd in bc_nodes[1]:
                x, y = get_node_coord(nd, nd_data)
                ax.plot(x, y, ls='', marker='^', ms=4, c='k')
            # x load -- arrow annotation towards right
            for nd in bc_nodes[2]:
                x, y = get_node_coord(nd, nd_data)
                ax.plot(x, y, ls='', marker='>', ms=4, c='r')
            # y load -- arrow annotation towards up
            for nd in bc_nodes[3]:
                x, y = get_node_coord(nd, nd_data)
                ax.plot(x, y, ls='', marker='^', ms=4, c='r')
        pass
     
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()
  
#%% [Sci charts, lineplots]
def plot_structural_response(srmat, **kwargs):    
    # reading kwargs
    xunit = kwargs.get('xunit', None)
    yunit = kwargs.get('yunit', None)
    validate = kwargs.get('validate', False)
    ref_loc = kwargs.get('ref_loc', None)
    
    if validate:
        ref_data = read_csv(ref_loc)
    
    fig, ax = plt.subplots()
    #fig.suptitle(f'Step-{len(srmat)-1}')
    
    # gridlines
    ax.grid(which='major', linestyle='--', lw=0.1, c='royalblue')
    ax.grid(which='minor', linestyle=':', lw=0.01, c='royalblue')

    # axis ranges
    if validate:
        ax.set_xlim(0, round(max(ref_data['disp']), 2)*1)
        ax.set_ylim(0, round(max(ref_data['load']), 0)*1.2)
        xRange = round(max(ref_data['disp']), 2)
        yRange = round(max(ref_data['load']), 0)
    else:
        ax.set_xlim(0, max(srmat[:, 0]))
        ax.set_ylim(0, max(srmat[:, 1]))
    
    # axis ticks
    if validate:
        ax.set_xticks(linspace(0, round(max(ref_data['disp']), 2), 5, dtype=float))
        ax.set_yticks(linspace(0, round(max(ref_data['load']), 0), 5, dtype=float))
    else:   
        ax.set_xticks(linspace(0, max(srmat[:, 0])*1, 5, dtype=float))
        ax.set_yticks(linspace(0, max(srmat[:, 1])*1, 5, dtype=float))
    
    # axis labels
    if xunit is not None:
        ax.set_xlabel(f'Displacement [{xunit}]')
    else:
        ax.set_xlabel(f'Displacement [mm]')
        
    if yunit is not None:
        ax.set_ylabel(f'Load [{yunit}]')
    else:
        ax.set_ylabel(f'Load [N]')
        
    # plot
    if validate:
        ax.plot(ref_data['disp'], ref_data['load'], ls='--', c='red', lw=1, label='Reference data')
    
    ax.plot(srmat[:, 0], srmat[:, 1], ls='-', c='blue', lw=1, label='PF framework')
    ax.plot(srmat[-1, 0], srmat[-1, 1], ls='', marker='o', ms=5, mfc='None')
    
    # legend
    ax.legend()
    
    # tight_layout and show func call
    # ax.set_aspect(xRange/yRange)
    plt.title(f'Step-{len(srmat)-1}')
    plt.tight_layout()
    plt.show()
    
def plot_time_statistics(tsmat, **kwargs):
    # essentially this is a bar plot in format of (2, 1)
    fig, ax = plt.subplots(2, 1)
    
    # numerical iterations
    ax[0].bar(tsmat[:, 0], tsmat[:, 1], c='red')
    ax[0].set_ylabel("Iterations")
    ax[0].grid(which="major", linestyle="--", lw=0.1, c="royalblue")
    ax[0].grid(which="minor", linestyle=":", lw=0.01, c="royalblue")
    ax[0].set_xlabel(f'Step no. {tsmat[-1, 0]}')

    # process time
    ax[1].bar(tsmat[:, 0], tsmat[:, 2], c='blue')
    ax[1].set_ylabel("Time")
    ax[1].grid(which="major", linestyle="--", lw=0.1, c="royalblue")
    ax[1].grid(which="minor", linestyle=":", lw=0.01, c="royalblue")
    ax[1].set_xlabel(f'Step no. {tsmat[-1, 0]}')
    
    plt.tight_layout()
    plt.show()

