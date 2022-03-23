import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt


plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
iconds = range(0,1000,2)
icond_avgs = []
for icond in iconds:
    files = glob.glob(F'../../5_namd/2x2_C3N4/namd_res_{icond}/*/SH_pop.txt')
    c = 0
    if len(files)>0:
        for file in files:
            x = np.loadtxt(file)
            #print(x.shape)
            if x.shape[0]==2000:
                c += 1
                if c==1:
                    x_ave = np.zeros(x.shape)
                x_ave += x
        x_ave /= c
        if np.min(x_ave[:,1])<1:
            #plt.plot(x_ave[:,1],color='gray')
            icond_avgs.append(x_ave)
        #print(c)
icond_avgs = np.array(icond_avgs)


def linear_func(t, tau):
    return 1-t/tau


taus = []
for i in range(len(icond_avgs)):
    tmp = np.min(icond_avgs[i][:,1])
    #if tmp<1:
    md_time = np.arange(0,len(icond_avgs[i][:,1]))
    popt, pcov = curve_fit( linear_func, md_time, icond_avgs[i][:,1],bounds=([0.0],[40000000]) )
    #print(popt)
    #tau = 1/(1-tmp)*2/1000
    tau = popt[0]
    residuals  = icond_avgs[i][:,1] - linear_func(md_time, *popt)
    ss_res     = np.sum(residuals**2)
    ss_tot     = np.sum((icond_avgs[i][:,1] - np.mean(icond_avgs[i][:,1]))**2)
    r_squared  = 1.0 - (ss_res / ss_tot)
    if r_squared>0.5:
        taus.append(tau)
        plt.plot(md_time, icond_avgs[i][:,1], color='gray')
        plt.plot(md_time, linear_func(md_time, *popt), color='blue',ls='--')
        print('R2:',r_squared,', tau:', tau/1000000,'ns')
taus = np.array(taus)
Z = 1.96
N = len(taus)
print(N)
s = np.std(taus)
tau_ave = np.average(taus)
error_bar = Z*s/np.sqrt(N)
print(tau_ave/1000000,'+-',error_bar/1000000, 'ns')
plt.plot(np.average(icond_avgs,axis=0)[:,1], color='red',lw=2)
plt.ylabel('Population', fontsize=10)
plt.xlabel('Time, fs', fontsize=10)
plt.title('C$_3$N$_4$ - e-h recombination', fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()
plt.savefig('../c3n4_e_h_recomb.jpg', dpi=600)

