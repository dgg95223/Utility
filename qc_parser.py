# Aurthor: Jingheng Deng
# E-mail: deng.jing.heng223@hotmail.com
import numpy as np
import re
import subprocess as sp

au2ev = 27.2114

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
            es['es_osci_str']  = es_info[8].split('=')[1]
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
        self.err_scf_converge = None
        self.err_cis_converge = None
        self.err_img_root_rpa = None
        self.have_error = None

    def get_opt_geom(self): 
        igeom = []
        output = self.out    
        for ii, i in enumerate(output):
            if 'NAtoms' in i:
                n_atom = int(output[ii+1].split()[0])
            if '**  OPTIMIZATION CONVERGED  **' in i:
                igeom.append(ii+5)

        ilast_start = igeom[-1]
        ilast_end   = ilast_start + n_atom
        last_geom = []
        atom_sym = []
        for i in range(ilast_start, ilast_end):
            last_geom.append(output[i].split()[2:])
            atom_sym.append(output[i].split()[1])
        geom     = last_geom
        return n_atom,atom_sym,np.array(geom,dtype=np.float64)

    def diagnolizer(self): # not finished
        output = self.out
        print('125',self.file_path[0:-4])
        for i in output:
            if 'Thank you' in i:
                self.have_error = False
            if 'SCF fail' in i:
                self.err_scf_converge = True
                print('SCF fail to converge.')
            if 'CIS/TDDFT calculation failed to' in i:
                self.err_cis_converge = True
                print('CIS/TDDFT fail to converge.')
            if 'Imaginary RPA root detected' in i: # generally if the S0 can converge easily while imaginary roots detected in direct TDDFT, the singlet/triplet MOs are unstable
                self.err_img_root_rpa = True
                print('Imaginary RPA root detected.')
            if 'Coordinates do not transform' in i: # don't really know the reason
                self.err_standard_ori = True
                print('Unkonwn error: coordinates do not transform.')
                
    def fixer(self):
        tag = self.file_path[0:-4]
        with open(tag+'.inp','r') as f:
            inp = f.readlines()
        job_ctrl = {'MAX_SCF_CYCLES':200,'MAX_CIS_CYCLES':200,'internal_stability':1,}
        status = [self.err_scf_converge, self.err_cis_converge, self.err_img_root_rpa]
        inp_ = ''.join(inp)
        for idx,key in enumerate(job_ctrl.keys()):
            if status[idx] is True:
                if key.lower() in inp_.lower():
                    sp.run(["sed -i '/ %s/c\\%s %d' %s.inp"%(key,key,job_ctrl[key],tag)],shell=True)
                else:
                    sp.run(["sed -i '/rem/a\\%s %d' %s.inp"%(key,job_ctrl[key],tag)],shell=True)
            status[idx] = False
            
        sp.run('. /opt/q-chem/qcenv.sh',shell=True)
        sp.run('qchem -nt 1 %s.inp %s.out'%(tag,tag),shell=True)
    
    def get_tddft_results(self):
        output = self.out
        ES_idx = []
        istate = []
        amp_state = 0   # range: 0,1,2
        pbht = 5

        itd   = None
        itda  = None
        not_converge = False
        set_tri  = 0     # range: 0,1
        set_sin  = 0 
        do_tri   = 0
        do_sin   = 0
        ptss_idx = [] #  f do state specific pcm correction with perturbation

        iosci = -1
        do_u  = 0 # unrestricted calculation or not, defaultly no.
        imulti = 2
        show_multi = 0
        ies_pair = 3

        for ii,i in enumerate(output):
            if 'cis_triplets' in i.lower():
                set_tri   = 1
                i_ = [j for j in re.split(' |\n',i) if j]
                if (i_[-1] == '1') or (i_[-1].lower() == 'true'):
                    do_tri     = 1
            if 'cis_singlets' in i.lower():
                set_sin   = 1 
                i_ = [j for j in re.split(' |\n',i) if j]
                if (i_[-1] == '1') or (i_[-1].lower() == 'true'):
                    do_sin     = 1
            if 'cis_n_roots' in i.lower():
                i_ = [j for j in re.split(' |\n',i) if j]
                nstate = int(i_[-1])              
            if 'Excited state' in i:
                ES_idx.append(ii)
                istate.append(i.split()[2].split(':')[0])
            if 'unrestricted' in i.lower():
                i_ = [j for j in re.split(' |\n',i) if j]
                if (i_[-1] == '1') or (i_[-1].lower() == 'true'):
                    pbht = 5
                    do_u = 1
                    ies_pair = 4
            if 'pbht_analysis' in i.lower():
                i_ = [j for j in re.split(' |\n',i) if j]
                if (i_[-1] == '1') or (i_[-1].lower() == 'true'):
                    if do_u == 0:
                        pbht = 6 # read pbht
                    else:
                        pbht = 5 # pbht doesn't show up when 'unrestricted' is true
            if 'Total 1st-order (ptSS+ptLR)  exc. energy' in i: # Excitation energy with state specific correction using perturbation
                ptss_idx.append(ii)                    
            if 'TDDFT/TDA Excitation Energies' in i:
                itda  = ii
            if 'TDDFT Excitation Energies' in i:
                itd  = ii
            if 'NRoots was altered' in i:
                nstate = int([j for j in re.split(' |-->|\n',i) if j][-1])
            if 'CIS/TDDFT calculation failed to' in i:
                not_converge = True            
            if 'Multiplicity' in i:
                show_multi = 1

        if not_converge is True:
            print('The cis/tddft job is not converged')

        if (set_sin == 0) and (set_tri == 0):
            print('Defaultly, both singlet and triplet CIS states are calculated, so both states will be read.')
            do_sin = 1
            do_tri = 1
        elif (set_sin == 0) and (set_tri == 1):
            do_sin = 1
        elif (set_sin == 1) and (set_tri == 0):
            do_tri = 1

        amp_state = do_sin + do_tri
        self.nstates = nstate * amp_state        

        if itda is not None:
            ES_tda_idx = [i for i in ES_idx if i < itd]
            ES_tda_ptss_idx = ptss_idx[0:self.nstates]
            assert len(ES_tda_idx) == self.nstates, 'Something wrong with the output file'
        ES_idx = [i for i in ES_idx if i > itd]  # defaultly only read tddft results
        ES_ptss_idx = ptss_idx[self.nstates:]
        assert len(ES_idx) == self.nstates, 'Something wrong with the output file'
            
        self.ES = {}
        for ii,i in enumerate(ES_idx): # read tddft
            es = {}
            es_info = [j.split('\n') for j in output[i:i+pbht+1]]
            es_pairs = []
            
            istart = 1
            
            idx = i+pbht-do_u
            while istart >= 1:          
                i_ = [j for j in re.split(' +|X\:|-->|D|V|\(|\)|amplitude|\=|\n',output[idx]) if j]
                if (i_ != []) and (len(i_) ==ies_pair):
                    if do_u == 1:
                        i_ = i_[0:3]
                    es_pairs.append(i_)
                    istart += 1
                    idx += 1
                else:
                    istart = 0

            es['es_pairs']              = np.float64(es_pairs)
            es['es_energy_ev']          = np.float64(es_info[0][0].split()[-1])                   # eV
            es['es_energy_ptss_ev']     = np.float64(output[ES_ptss_idx[ii]].split()[-2])         # ptLR + ptSS, eV 
            es['es_tot_energy']         = np.float64(es_info[1][0].split()[-2]) * au2ev           # eV
            if show_multi == 1:
                es['es_multi']         = es_info[imulti][0].split()[-1]
            es['es_osci_str']      = np.float64([i for i in es_info[imulti + 2][0].split() if i][iosci])
            if pbht == 6:
                es['es_pbht']      = np.float64([i for i in es_info[imulti + 3][0].split() if i][-1])
            
            self.ES[istate[ii]] = es

        self.ES['nstates'] = self.nstates

        return self.ES
    
    def get_dft_results(self):
        output = self.out
        self.GS = {}
        have_sym = False
        for ii,i in enumerate(output):
            if 'beta electrons' in i:
                n_elec = np.array([j for j in re.split(' +|There are|alpha and|beta electrons|\n',i) if j ],dtype=np.int64)
            if 'Total energy in the final basis set' in i:
                energy = np.float64(i.split()[-1]) * au2ev
            if 'Alpha MOs' in i:
                iMO_occ_start = ii + 2
            if '-- Virtual --' in i:
                iMO_occ_end   = ii
                iMO_vir_start = ii + 1
            if 'basis functions' in i: 
                n_MO = np.array([j for j in re.split(' |\n',i) if j ][-3],dtype=np.int64)
            if 'Orbital Energies (a.u.) and Symmetries' in i:
                have_sym = True
            if 'Number of orthogonalized atomic orbitals' in i: # The format change when the symmetry information is print out
                n_MO = np.array([j for j in re.split(' |\n',i) if j ][-1],dtype=np.int64)
        n_MO_occ = n_elec[0]
        n_MO_occ_ = n_MO_occ // 8
        if n_MO_occ_ * 8 < n_MO_occ:
            n_MO_occ_ += 1
        n_MO_vir = n_MO - n_MO_occ
        n_MO_vir_ = n_MO_vir // 8
        if n_MO_vir_ * 8 < n_MO_vir:
            n_MO_vir_ += 1
        iMO_vir_end = iMO_vir_start + n_MO_vir_

        MOs_occ = []
        MOs_vir = []
        if have_sym is not True:
            for i in range(iMO_occ_start,iMO_occ_end):            
                MOs_occ.append(np.array(output[i].split(),dtype=np.float64))
            for i in range(iMO_vir_start,iMO_vir_end):
                MOs_vir.append(np.array(output[i].split(),dtype=np.float64))
        else:
            for i in range(0, n_MO_occ_):
                i_ = iMO_occ_start + 2*i
                MOs_occ.append(np.array(output[i_].split(),dtype=np.float64))
            for i in range(0 ,n_MO_vir_):
                i_ = iMO_vir_start + 2*i
                MOs_vir.append(np.array(output[i_].split(),dtype=np.float64))
        MOs_occ = np.concatenate(MOs_occ,dtype=np.float64) * au2ev
        MOs_vir = np.concatenate(MOs_vir,dtype=np.float64) * au2ev

        e_homo = MOs_occ[-1]
        e_lumo = MOs_vir[0]

        self.GS['gs_energy'] = energy   # eV
        self.GS['n_elecs']   = n_elec
        self.GS['E_mo_occ']  = MOs_occ  # eV
        self.GS['E_mo_vir']  = MOs_vir  # eV
        self.GS['E_homo']    = e_homo
        self.GS['E_lumo']    = e_lumo

        return self.GS
    
    def get_dm(self,dm_a=True,dm_b=False):
        output = self.out

        dm_a_tags = []
        dm_b_tags = []
        for ii,i in enumerate(output):
            if 'basis functions' in i:
                n_MO = np.array([j for j in re.split(' |\n',i) if j ][-3],dtype=np.int64)
            if (dm_a is True) and ('Final Alpha density matrix' in i):
                dm_a_tags.append(ii + 1)
            if (dm_b is True) and ('Final Beta density matrix' in i):
                dm_b_tags.append(ii + 1)

        n_col = 6
        n_block = n_MO // 6
        if n_block * n_col < n_MO:
            n_block += 1        
        
        if (dm_a is True) and (dm_b is False):
            n_dm = len(dm_a_tags)
            dm_shape = (n_dm,n_MO,n_MO)
            dms = np.zeros(dm_shape)
            for idm in range(0,n_dm):
                icol = 0
                for i in range(0,n_block):
                    istart = dm_a_tags[idm] + i * (n_MO + 1)
                    n_col = len(output[istart+1].split())-1
                    for j in range(1,n_col+1):
                        dms[idm,:,icol+j-1] = np.loadtxt(output[istart+1:istart+n_MO+1],dtype=np.float64)[:,j]
                    icol += n_col
        elif (dm_a is True) and (dm_b is True):
            n_dm = len(dm_a_tags)
            dm_shape = (n_dm,2,n_MO,n_MO)
            dms = np.zeros(dm_shape)
            for idm in range(0,n_dm):
                icol_a = 0
                icol_b = 0
                for i in range(0,n_block):
                    istart_a = dm_a_tags[idm] + i * (n_MO + 1)
                    istart_b = dm_b_tags[idm] + i * (n_MO + 1)
                    n_col_a = len(output[istart_a+1].split())-1
                    n_col_b = len(output[istart_b+1].split())-1
                    for j in range(1,n_col_a+1):
                        dms[idm,0,:,icol_a+j-1] = np.loadtxt(output[istart_a+1:istart_a+n_MO+1],dtype=np.float64)[:,j]
                    for j in range(1,n_col_b+1):
                        dms[idm,1,:,icol_b+j-1] = np.loadtxt(output[istart_b+1:istart_b+n_MO+1],dtype=np.float64)[:,j]
                    icol_a += n_col_a
                    icol_b += n_col_b
        
        return dms        
        
            
    def get_cdft_energy(self):
        return 
    def get_cdftci_coupling(self):
        '''
        Read CDFT-CI coupling from single output file, return the coupling value, str
        '''
        output = self.out
        success = False
        for ii,i in enumerate(output):
            if i == ' CDFT-CI Hamiltotnian matrix in orthogonalized basis\n':
                coupling = output[ii+2].split()[-1]
                success = True
        if not success:
            ec = '%15s'%'NaN'
        else:
            ec =  '%10.6f'%(np.float64(coupling)*au2ev)
        return ec 
    def get_binding_energy(self):
        '''
        Read binding energy and BSSE errors from BSSE jobs, dict
        '''
        return
if __name__ == '__main__':
    dft_ = '/home/jingheng/MO_DGNN/data/qm9/b3lyp_631+gd/input/dft/qm9_dft_111111.out'
    tddft_ = '/home/jingheng/MO_DGNN/data/qm9/b3lyp_631+gd/input/tddft/qm9_tddft_111111.out'

    dft = QChem(dft_).get_dft_results()
    print(dft)

    tddft = QChem(tddft_).get_tddft_results()
    print(tddft)
    
