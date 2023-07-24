# Aurthor: Jingheng Deng
# E-mail: deng.jing.heng223@hotmail.com
import numpy as np
import re

class Gaussian():
    '''
    This is a class of methods for parsing output files, '.log', of Gaussian package
    '''
    def __init__(self, file_path):
        self.file_path = file_path
        with open(file_path,'r') as f:
            self.out = f.readlines()

    def get_gaussian_opt_geom(self):
        igeom = []
        output = self.out    
        for ii, i in enumerate(output):
            if 'NAtoms' in i:
                n_atom = int(i.split()[1])
            if 'Standard orientation' in i:
                igeom.append(ii+5)
        ilast_start = igeom[-1]
        ilast_end   = ilast_start + n_atom
        last_geom = []
        atom_sym = []
        for i in range(ilast_start, ilast_end):
            last_geom.append(output[i].split()[3:])
            atom_sym.append(output[i].split()[1])
        geom     = last_geom
        return n_atom,atom_sym,np.array(geom,dtype=np.float64)
    
    def get_tddft_results(self):
        output = self.out
        ES_idx = []
        istate = []

        for ii,i in enumerate(output):
            if 'nstates' in i:
                i_ = re.split(' |\/|\=|\+|\(|\)|\,',i)
                idx = i_.index('nstates')
                self.nstates = i[idx+1]
            if 'Excited State' in i:
                ES_idx.append(ii)
                istate.append(i.split()[2].split(':')[0])
                
        self.ES = {}
        for ii,i in enumerate(ES_idx):
            es = {}
            es_info = output[i].split()
            es_pairs = []
            
            istart = 1
            idx = ES_idx[ii]+1
            while istart >= 1:          
                i_ = [j for j in re.split(' +|->|\n',output[idx]) if j]
                if (i_ != []) and (len(i_) ==3):
                    es_pairs.append(i_)
                    istart += 1
                    idx += 1
                else:
                    istart = 0

            es['es_pairs']     = es_pairs
            es['es_energy_ev'] = es_info[4]
            es['es_energy_nm'] = es_info[6]
            es['es_freq']      = es_info[8].split('=')[1]
            es['es_s2']        = es_info[9].split('=')[1]
            
            self.ES[istate[ii]] = es

        return self.ES

class QChem():# not finished
    '''
    This is a class of methods for parsing output files, '.out', of Q-Chem package
    '''
    def __init__(self, file_path):
        self.file_path = file_path
        with open(file_path,'r') as f:
            self.out = f.readlines()

    def get_opt_geom(self):
        igeom = []
        output = self.out    
        for ii, i in enumerate(output):
            if 'NAtoms' in i:
                n_atom = int(i.split()[1])
            if 'Standard orientation' in i:
                igeom.append(ii+5)
        ilast_start = igeom[-1]
        ilast_end   = ilast_start + n_atom
        last_geom = []
        atom_sym = []
        for i in range(ilast_start, ilast_end):
            last_geom.append(output[i].split()[3:])
            atom_sym.append(output[i].split()[1])
        geom     = last_geom
        return n_atom,atom_sym,np.array(geom,dtype=np.float64)
    
    def get_tddft_results(self):
        output = self.out
        ES_idx = []
        istate = []

        for ii,i in enumerate(output):
            if 'nstates' in i:
                i_ = re.split(' |\/|\=|\+|\(|\)|\,',i)
                idx = i_.index('nstates')
                self.nstates = i[idx+1]
            if 'Excited State' in i:
                ES_idx.append(ii)
                istate.append(i.split()[2].split(':')[0])
                
        self.ES = {}
        for ii,i in enumerate(ES_idx):
            es = {}
            es_info = output[i].split()
            es_pairs = []
            
            istart = 1
            idx = ES_idx[ii]+1
            while istart >= 1:          
                i_ = [j for j in re.split(' +|->|\n',output[idx]) if j]
                if (i_ != []) and (len(i_) ==3):
                    es_pairs.append(i_)
                    istart += 1
                    idx += 1
                else:
                    istart = 0

            es['es_pairs']     = es_pairs
            es['es_energy_ev'] = es_info[4]
            es['es_energy_nm'] = es_info[6]
            es['es_freq']      = es_info[8].split('=')[1]
            es['es_s2']        = es_info[9].split('=')[1]
            
            self.ES[istate[ii]] = es

        return self.ES

    def get_cdft_energy(self):
        return 
    def get_cdftci_coupling(self):
        return
