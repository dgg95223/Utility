from turtle import pu
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np

def plot(x, y, plot_type='single', plot_mode='line', **setting):
    '''General basic plot setting for plot and subplots
    x, data for X Axis, 1-D array
    y, data for Y Axis, for single plot, 1-D array; for N subplots, N-D array
    
    '''
    default_key = ['title', 'xlabel', 'ylabel', 'plot_name']
    advanced_key = ['xlim', 'ylim']
    default_dict = {'title':'notitle', 'xlabel':'X Axis', 'ylabel':'Y Axis', 'plot_name':'noname.png'}

    input_keys = []

    for key in setting:
        input_keys.append(key)

    for key in input_keys:
        if key not in default_key:
            setting[key] = default_dict[key]

    if plot_type  == 'single':
        fig, ax = plt.subplots()
        if plot_mode == 'line':
            ax.plot(x,y)
        elif plot_mode == 'scatter':
            ax.scatter(x,y)
        else:
            raise NotImplemented
        ax.set_xlabel(setting['xlabel'])
        ax.set_ylabel(setting['ylabel'])
        ax.set_title(setting['title'])
        if 'xlim' in input_keys:
            ax.set_xlim(setting['xlim'])
        elif 'ylim' in input_keys:
            ax.set_ylim(setting['ylim'])
        
        fig.savefig(setting['plot_name'])

    elif plot_type == 'sub':
        plots_shape = setting['plots_shape']
        n_plot = len(y)
        n_col_max = plots_shape[1]
        if n_plot <= 3:
            n_row = 1
            n_col = n_plot
        elif n_plot > 3:
            n_row = n_plot // n_col_max + 1
            n_col = n_col_max
        
        fig, axs = plt.subplots(n_row, n_col, sharex='all', sharey='all')
        axs.set_title(setting['title'])
        for i in range(0, n_row):
            for j in range(0, n_col):
                if plot_mode == 'line':
                    axs[i,j].plot(x,y[i*n_col_max+j])
                elif plot_mode == 'scatter':
                    axs[i,j].scatter(x,y[i*n_col_max+j])                    
                else:
                    raise NotImplemented
                axs[i,j].set_xlabel(setting['xlabel'])
                axs[i,j].set_ylabel(setting['ylabel'])
                if 'xlim' in input_keys:
                    axs[i,j].set_xlim(setting['xlim'])
                elif 'ylim' in input_keys:
                    axs[i,j].set_ylim(setting['ylim'])
        
        fig.savefig(setting['plot_name'])
    else:
        raise NotImplemented

def printm(data):
    '''Print a matrix in a neat  way'''
    return NotImplemented

def printcm(data):
    '''Print a complex matrix in a neat way'''
    return NotImplemented

def read_txt(filename, usecol=0):
    return NotImplemented

def read_xyz(filename, index=None, output='regular'):
    '''
    index: '-1' refers to the last geometry
           'N' any integar larger than 0, refers to the N^th geometry
           '0' refers to all geometry
    output mode: 'regular' output atom number, atom symbols, a np array of coordinates
                 'pyscf' output atom number, atom symbols, a string include atom symbols and coordinates    
    '''
    if index is None:                                                                                           
        index = -1  # read the last geometry as default
    elif index > 0:
        index = index - 1
    elif index == 0:
        read_all = True

    with open(filename,'r') as xyz:
        molecules = xyz.readlines()
    
    for i in range(0, len(molecules)):
        if molecules[i] == '\n':
            molecules.pop(i)
    
    atoms_num = int(molecules[0])
    atom_symbol = np.loadtxt(filename, usecols=0, dtype='str', skiprows=2, max_rows=atoms_num+2)

    assert len(atom_symbol) == atoms_num, 'There is somthing wrong with the format of xyz file'
    
    if index == -1:
        # read the last geometry
        if output == 'regular':
            _geom = []
            for i in range(0, atoms_num):
                _geom_ = molecules[index * (atoms_num + 2) + 2 + i].split()
                _geom.append(_geom_[1:])
            _geom =np.array(_geom, dtype=np.float64)
        elif output == 'pycsf':
            _geom = ''.join(molecules[index * atoms_num:])
    elif index == 0:  ##################################################### to be done 2022/4/25
        # read all geometiers
        geoms = []
        geoms_num = len(np.loadtxt(filename, usecols=0, dtype='str')) // (atoms_num + 1)
        for i in range(0, geoms_num):
            _geom = ''.join(molecules[index * atoms_num:])
            geoms.append(_geom)
    elif molecules[(index + 1) * (atoms_num + 2) - 1][-1] == '\n':
        # read the specified geometry
        _geom = ''.join(molecules[index * (atoms_num + 2) + 2 :(index + 1) * (atoms_num + 2)])[:-1]
    else:
        _geom = ''.join(molecules[index * (atoms_num + 2) + 2 :(index + 1) * (atoms_num + 2)])


        
    
    return 

def read_json(filename):
    '''from json file to python dict'''
    import json
    with open(filename,'r') as file:
        json_setting = json.load(file)
    return json_setting

def write_json(json_dict, filename):
    '''from python dict to json file'''
    import json
    with open(filename,'w') as file:
        json.dump(json_dict, file, indent=4)