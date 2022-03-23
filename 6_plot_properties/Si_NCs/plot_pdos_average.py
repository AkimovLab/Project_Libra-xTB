import os
import sys
import time
import math
import glob
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
from libra_py import units

def gaussian_function(a, mu, sigma, num_points, x_min, x_max):

    pre_fact = (a/sigma)/(np.sqrt(2*np.pi))
    x = np.linspace(x_min, x_max, num_points)
    x_input = np.array((-1/2)/(np.square(sigma))*np.square(x-mu))
    gaussian_fun = pre_fact*np.exp(x_input)

    return x, gaussian_fun

def gaussian_function_vector(a_vec, mu_vec, sigma, num_points, x_min, x_max):
    
    for i in range(len(a_vec)):
        if i==0:
            sum_vec = np.zeros(num_points)
        energy_grid, conv_vec = gaussian_function(a_vec[i], mu_vec[i], sigma, num_points, x_min, x_max)
        sum_vec += conv_vec
    return energy_grid, sum_vec


# In[5]:


# %matplotlib notebook
folders = ['Si1009H412','Si501H228','Si329H172','Si265H140','Si123H100','Si59H60']
labels = ['Si$_{1009}$H$_{412}$','Si$_{501}$H$_{228}$','Si$_{329}$H$_{172}$',
          'Si$_{265}$H$_{140}$','Si$_{123}$H$_{100}$','Si$_{59}$H$_{60}$']
params = {}
params['sigma'] = 0.05
params['energy_shift'] = 3.0
params['atoms_list'] = ['Si', 'H']
# for i in [500]:
params['compute_average'] = True
params['path_to_pdos'] = F{'../../3_overlaps/{folder}/step2/all_pdosfiles'}
# plot_pdos(params)
# def plot_pdos(params):
sigma = params['sigma']
shift = params['energy_shift']
atoms = params['atoms_list']
path_to_pdos = params['path_to_pdos']
compute_average = params['compute_average']
for c1, folder in enumerate(folders):    
    # plot_title = params['plot_title']
    plt.figure(num=None, figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
    pdos_all = []
    for kind in [1,2]:
        if compute_average:
            pdos_files = sorted(glob.glob(path_to_pdos+F'/*k{kind}*pdos'))
            for c, pdos_file in enumerate(pdos_files):
                if c==0:
                    # initialize pdos average
                    tmp = np.loadtxt(pdos_file)
                    pdos_average = np.zeros(tmp.shape)
                pdos_average += np.loadtxt(pdos_file)
            pdos = pdos_average/(c+1)
        else:        
            pdos_files = sorted(glob.glob(path_to_pdos+F'/*k{kind}*pdos'))
            pdos = np.loadtxt(pdos_files[0])
        pdos_all.append(pdos)
    
    energy = pdos_all[0][:,1] * units.au2ev
    e_min = energy[0] - shift
    e_max = energy[-1] + shift
    # The number of points must be more than the number of states
    npoints = pdos_all[0].shape[0] + 2000
    orbitals = ['s','p','d','f']
    #                 S        P           D            F
    orbitals_cols = [[3], range(4,7), range(7,12), range(12,19)]
    homo_level = np.max(np.where(pdos_all[0][:,2]==2.0))
    homo_energy = energy[homo_level]
    total_pdos = np.zeros((npoints))
    for i in range(len(pdos_all)):
        for c, orbital_cols in enumerate(orbitals_cols):
            try:
                sum_pdos = np.sum(pdos_all[i][:,orbital_cols],axis=1)
                energy_grid, pdos_convolved = gaussian_function_vector(sum_pdos, energy, sigma,
                                                                               npoints, e_min, e_max)
                total_pdos += pdos_convolved
                pdos_label = atoms[i]+F', {orbitals[c]}'
                plt.plot(energy_grid-homo_energy, pdos_convolved, label=pdos_label)
                print('Done with sum', i, orbital_cols)
            except:
                pass
    plt.plot(energy_grid-homo_energy, total_pdos, label='Total', color='black')
    # plt.ylim(0,700)
    plt.xlim(-5,5)
    plt.xlabel('Energy, eV', fontsize=12)
    plt.ylabel('DOS, 1/eV', fontsize=12)
    plt.xticks([-5,-2.5,0.0,2.5,5.0],fontsize=12)
    plt.yticks(fontsize=12)
    plt.title(F'{labels[c1]}', fontsize=12)
    # plt.legend(fontsize=6, loc='best') #, bbox_to_anchor=(0.15, 0.5, 0.5, 0.15), frameon=True)#upper center')
    plt.tight_layout()
    plt.savefig(F'{folder}-ave-pdos.jpg', dpi=600)

