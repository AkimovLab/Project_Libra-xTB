import os
import sys
import time
import glob
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from liblibra_core import *
from libra_py import units, data_stat

#plt.figure(num=None, figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
fig, ax = plt.subplots(num=None, figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
folders = ['Si1009H412','Si501H228','Si329H172','Si265H140','Si123H100','Si59H60']
isteps = [0,2000,2000,2000,2000,2000]
# 1.0 eV above LUMO
active_spaces = [349, 182, 133, 116, 60, 32]
#active_spaces = [76,40,35,35,20,13]
labels = ['Si$_{1009}$H$_{412}$','Si$_{501}$H$_{228}$','Si$_{329}$H$_{172}$',
          'Si$_{265}$H$_{140}$','Si$_{123}$H$_{100}$','Si$_{59}$H$_{60}$']

data = []
for c, folder in enumerate(folders):
    print(folder)
    a = list(range(active_spaces[c]))
    nac = []
    for k in range(1,3999):
        print(k)
        hvib_ave = sp.load_npz(F'../../4_nacs/{folder}/res-electron-only/Hvib_sd_{k+isteps[c]}_im.npz')[a,:][:,a]
        hvib_ave_dense = hvib_ave.todense().real
        for i in range(hvib_ave.shape[0]):
            for j in range(hvib_ave.shape[0]):
                if j != i:
                    nac_ij = np.abs(hvib_ave_dense[i,j])* 1000.0 * units.au2ev
                    x_mb = MATRIX(1,1)
                    x_mb.set(0, 0, nac_ij )
                    nac.append( x_mb )
    bin_supp, dens, cum = data_stat.cmat_distrib( nac, 0, 0, 0, 0, 50, 0.1)
    plt.plot( bin_supp, dens, label=labels[c], linewidth=3.5 )
plt.xlim(0,4)
# plt.ylim(0.0,0.03)
#plt.legend(fontsize=10)
plt.xlabel('|NAC|, meV', fontsize=10)
plt.ylabel('PD, 1/meV', fontsize=10)
plt.tight_layout()
plt.savefig('nac_dist_1_3.jpg', dpi=600)

plt.xlim(4,10)
plt.ylim(0.0,0.08)
ax.get_legend().remove()
plt.tight_layout()
plt.savefig('nac_dist_2_3.jpg', dpi=600)

plt.xlim(10,50)
plt.ylim(0.0,0.02)
plt.tight_layout() 
plt.savefig('nac_dist_3_3.jpg', dpi=600)






