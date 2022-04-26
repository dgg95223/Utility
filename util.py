from turtle import pu
import matplotlib
from matplotlib.backend_bases import MouseEvent
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
           'N' any integar larger than 0, refers to the N^th geometry, '-' refers to count the geometry in reversed order
           '0' refers to all geometry
    output mode: 'regular' output atom number, atom symbols, a np.array of coordinates
                 'pyscf' output atom number, atom symbols, a string includes atom symbols and coordinates  
    
    Current version only support the geometries of the same molecule --2022/4/26
    '''
    with open(filename,'r') as xyz:
        molecules = xyz.readlines()
    
    # clear unnecessary empty rows
    reverse_i = list(range(0, len(molecules)))[::-1]
    for i in reverse_i:
        if molecules[i] == '\n':
            if (len(molecules[i-1]) > 10) or (len(molecules[i-1]) == 1):
                molecules.pop(i)

    # get the number of atoms in each geometry
    atoms_num = []
    ii = 0
    while ii < len(molecules) :
        atoms_num.append(int(molecules[ii]))
        ii += (2 + int(molecules[ii]))
        if ii == len(molecules):
            break

    # get the amount of geometries
    geoms_num = len(atoms_num)
    atom_symbol = []
    # get the symbol of atoms in each geometry
    _atom_symbol = np.loadtxt(filename, usecols=0, dtype='str')
    start = 1
    for i in range(0, geoms_num):    
        end = start + atoms_num[i]
        atom_symbol.append(_atom_symbol[start:end])
        start = end + 1

    if index is None:                                                                                           
        index = -1  # read the last geometry as default
    elif index > 0: # read the N^th geometry
        index = index - 1
    elif index <= -1:
        index = geoms_num + index 
    elif index == 0: # read all geometries
        pass
    
    if index == 0:
        # read all geometries
        geoms = []
        for i in range(0, geoms_num):
            if output == 'regular':
                _geom = []
                for j in range(0, atoms_num[i]):
                    _geom_ = molecules[ sum(atoms_num[:i+1]) +(i + 1) * 2 + j].split()
                    _geom.append(_geom_[1:])
                    print(_geom_)
                # print(_geom)
                _geom =np.array(_geom, dtype=np.float64)
            elif output == 'pyscf':
                _geom = ''.join(molecules[i * (atoms_num[i] + 2) + 2:(i + 1) * (atoms_num[i] + 2)])
            geoms.append(_geom)
    else: 
        # index == 'N' read the N^th geometry
        if output == 'regular':
            _geom = []
            for j in range(0, atoms_num[index][0]):
                _geom_ = molecules[index * (atoms_num[index] + 2) + 2 + j].split()
                _geom.append(_geom_[1:])
            _geom =np.array(_geom, dtype=np.float64)
        elif output == 'pyscf':
            _geom = ''.join(molecules[index * (atoms_num[index] + 2) + 2:(index + 1) * (atoms_num[index] + 2)])
        geoms = _geom
        atoms_num = atoms_num[index]
        atom_symbol =atom_symbol[index]
    
    return atoms_num, atom_symbol, geoms

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