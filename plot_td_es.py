import numpy as np
import re,glob
import matplotlib.pyplot as plt

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

file_paths = glob.glob('./*log')
for ifile,file_path in enumerate(file_paths):
    Output = Gaussian(file_path)
    results = Output.get_tddft_results()
    
    istate = np.array(list(results.keys()))
    es_pairs_amp = [np.array(results[i]['es_pairs'],dtype=np.float64)[:,2] for i in istate]
    es_pairs_idx = [np.array(results[i]['es_pairs'])[:,:2] for i in istate]
    
    # print(es_pairs_amp)
    es_enes = np.array([results[i]['es_energy_ev'] for i in results.keys()],dtype=np.float64)
    es_freqs = np.array([results[i]['es_freq'] for i in results.keys()],dtype=np.float64)
    imax_freq = np.argmax(es_freqs)
    imax_pairs_amp = [np.argmax(i) for i in es_pairs_amp]
    # print(imax_pairs_amp)

    max_es_pairs_amp = es_pairs_amp[imax_freq][imax_pairs_amp[imax_freq]]
    max_es_pairs_idx = es_pairs_idx[imax_freq][imax_pairs_amp[imax_freq]]
    # print(max_es_pairs_idx)

    # adding trendline with polyfit
    
    z = np.polyfit(es_enes, es_freqs,5)
    p = np.poly1d(z)
    x = np.arange(min(es_enes),max(es_enes),0.01)
    
    fig,ax = plt.subplots(figsize=(12,8))
    plt.rcParams.update({'font.size':18})
    ax.stem(es_enes,es_freqs,basefmt='black')
    ax.set_ylim(0,round(max(es_freqs),1)+0.1)
    ax.set_xlim(min(es_enes),max(es_enes))
    ax.text(min(es_enes)+0.03,round(max(es_freqs),1)+0.075,'Max oscillation strength: %5.3f'%es_freqs[imax_freq])
    ax.annotate('State %d: %s:  %5.3f'%(imax_freq+1,'->'.join(max_es_pairs_idx),max_es_pairs_amp),\
                (es_enes[imax_freq]-0.35,max(es_freqs)+0.02))
    ax.plot(x, p(x)*3, "r")
    ax.set_xlabel('Excitation energy',fontsize=18)
    ax.set_ylabel('Oscillation strength',fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.savefig('%s.png'%file_path.split('.')[1][1:])
