import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

dirs = ['../5_namd/2x2_C3N4']
titles = ['C$_3$N$_4$']
R2s = [0.1,0.2,0.3,0.4,0.5,0.6]
nsteps = 2000
for ct, dir1 in enumerate(dirs):
    iconds = range(0,1000,2)
    ist = 0
    icond_avgs = []
    print('Reading data of dir', dir1)
    for icond in iconds:
        files = glob.glob(F'{dir1}/namd_res_{icond}/*/SH_pop.txt')
        c = 0
        if len(files)>0:
            for file in files:
                x = np.loadtxt(file)
                if x.shape[0]==nsteps:
                    c += 1
                    if c==1:
                        x_ave = np.zeros(x.shape)
                    x_ave += x
            x_ave /= c
            if np.min(x_ave[:,0])<1:
                icond_avgs.append(1-x_ave)
    icond_avgs = np.array(icond_avgs)
    min_val = 100
    for R2 in R2s:
        plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
        
        def linear_func(t, tau):
            return 1-t/tau
        def exp_func(t, tau):
            return np.exp(-t/tau)
        taus = []
        x_ave_1 = np.zeros(x_ave.shape)
        c1 = 0
        for i in range(len(icond_avgs)):
            tmp = np.min(icond_avgs[i][:,ist])
            #if tmp<1:
            md_time = np.arange(0,len(icond_avgs[i][:,ist]))
            #popt, pcov = curve_fit( linear_func, md_time, icond_avgs[i][:,ist],bounds=([0.0],[40000000]) )
            popt, pcov = curve_fit( exp_func, md_time, icond_avgs[i][:,ist],bounds=([0.0],[40000000]) )
            #print(popt)
            #tau = 1/(1-tmp)*2/1000
            tau = popt[0]
            #residuals  = icond_avgs[i][:,ist] - linear_func(md_time, *popt)
            residuals  = icond_avgs[i][:,ist] - exp_func(md_time, *popt)
            ss_res     = np.sum(residuals**2)
            ss_tot     = np.sum((icond_avgs[i][:,ist] - np.mean(icond_avgs[i][:,ist]))**2)
            r_squared  = 1.0 - (ss_res / ss_tot)
            if r_squared>R2:
                c1 += 1
                taus.append(tau)
                plt.plot(md_time, icond_avgs[i][:,ist], color='gray')
                min_val = min([min_val, np.min(icond_avgs[i][:,ist])])
                print('R2:',r_squared,', tau:', tau/1000000,'ns')
                x_ave_1 += icond_avgs[i]
        x_ave_1 /= c1
        taus = np.array(taus)
        Z = 1.96
        N = len(taus)
        s = np.std(taus)
        tau_ave = np.average(taus)
        error_bar = Z*s/np.sqrt(N)
        print('For ', titles[ct])
        print('Average timescale:',tau_ave/1000000,'+-',error_bar/1000000, 'ns', 'R2:',R2, 'Fitted ', N, 'samples')
        print(np.average(icond_avgs,axis=0)[:,ist])
        plt.plot(x_ave_1[:,ist],color='red',lw=2)
        plt.ylabel('Population', fontsize=10)
        plt.xlabel('Time, fs', fontsize=10)
        plt.text(-20,min_val,'%0.2f$\pm$%0.2f ns'%(tau_ave/1000000, error_bar/1000000), fontsize=14, fontweight='bold')
        plt.title(titles[ct], fontsize=10)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.tight_layout()
        plt.savefig(F'c3n4_e_h_recomb_{titles[ct]}_{R2}_.jpg', dpi=600)
        plt.close()
