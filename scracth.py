from file_func import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from numpy import array, asarray, linspace

elem_mat = array([[1, 1, 2, 3, 4],
                  [2, 4, 3, 5, 0]], dtype=int)

node_mat = array([[1, 0, 0],
                  [2, 1, 0],
                  [3, 1, 1],
                  [4, 0, 1],
                  [5, 0.5, 2]], dtype=float)

crack_set = array([1, 2])
fill_set1 = array([2,4,5])
fill_set_data = [crack_set, fill_set1]

velocity = array([1, 2, -1, 3, 0.2])



def plot_contour(**kwargs ):
    # use kwargs for user input 
    plot_this = kwargs.get('plot_this', None)   # vector value
    plot_on = kwargs.get('plot_on', None)   # character value
    crds = kwargs.get('crds', None)         # tensor value
    text_TT = kwargs.get('title', None)     # character value
    stpNum = kwargs.get('stpNum', None)     # Integer value
    fill_data = kwargs.get('fill_data', None)   # tensor data
    
    
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
            raise ValueError('coordinate system undefiuned')
        else:
            t = tri.Triangulation(crds[:,1], crds[:,2])
        
        # plot begins here
        fig, ax = plt.subplots()
        ax.axis('off')
        
        # Fill sets
        if fill_data is not None:
            # filling crack
            crck_x = [get_node_coord(i, crds)[0] for i in fill_data[0]]
            crck_y = [get_node_coord(i, crds)[1] for i in fill_data[0]]
            
            ax.fill(crck_x, crck_y, fill=True, color='white', edgecolor='white', lw=2)
        
        cntr = ax.tricontourf(t, plot_this, cmap='jet', levels=500,
                              vmin=plot_this.min(), vmax=plot_this.max())
        
        fig.colorbar(cntr, ax=ax, 
                     ticks=linspace(plot_this.min(), plot_this.max(), 10), location='bottom')
        
        plt.title(f'{text_TT} - Step-{stpNum}')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()
    

plot_contour(plot_this=velocity, plot_on='nodes', crds=node_mat, title='Velocity', stpNum=1,
             fill_data=fill_set_data)
    
