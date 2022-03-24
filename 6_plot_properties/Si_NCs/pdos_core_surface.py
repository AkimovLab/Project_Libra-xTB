import os
import numpy as np
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




names = ['Si59H60', 'Si123H100', 'Si265H140', 'Si329H172', 'Si501H228', 'Si1009H412']
titles = ['Si$_{59}$H$_{60}$', 'Si$_{123}$H$_{100}$', 'Si$_{265}$H$_{140}$',
          'Si$_{329}$H$_{172}$', 'Si$_{501}$H$_{228}$', 'Si$_{1009}$H$_{412}$']
npoints = 5000
shift = 2.0 # eV
sigma = 0.1 # eV
labels = ['Surface', 'Core']
colors = ['blue','red']
for c1, name in enumerate(names):
    title = titles[c1]
    plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
    for c, i in enumerate([1,2]):
        pdos = np.loadtxt(F'{name}/Si_QD_xtb-si_qd-list{i}-1.pdos')
        pdos[:,1] *= units.au2ev
        e_min = np.min(pdos[:,1])-shift
        e_max = np.max(pdos[:,1])+shift
        homo_level = np.max(np.where(pdos[:,2]==2.0))
        homo_energy = pdos[:,1][homo_level]
        sum_pdos = np.sum(pdos[:,3::],axis=1)
        energy_grid, pdos_convolved = gaussian_function_vector(sum_pdos, pdos[:,1], sigma, npoints, e_min, e_max)
        plt.plot(energy_grid-homo_energy, pdos_convolved, label=labels[c], color=colors[c])
    plt.xlim(-5,5)
    if c1==0:
        plt.legend(fontsize=12)
    plt.ylabel('DOS, 1/eV',fontsize=12)
    plt.xlabel('Energy, eV',fontsize=12)
    plt.xticks([-5,-2.5,0,2.5,5],fontsize=12)
    plt.yticks(fontsize=12)
    plt.title(title,fontsize=12)
    plt.tight_layout()
    plt.savefig(F'{name}.jpg',dpi=300)
    plt.close()

