import os
import sys
import time
import glob
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from libra_py import units


# # Find all the folders for a system and decoherence scheme
# We first need to find all the folders that contain the NAD results. In order to do that, we use Linux `find` command and then read them through Python.

# In[130]:

folder_names = ['Si59H60', 'Si123H100', 'Si265H140', 'Si329H172', 'Si501H228', 'Si1009H412']
titles = ['Si$_{59}$H$_{60}$', 'Si$_{123}$H$_{100}$', 'Si$_{265}$H$_{140}$', 'Si$_{329}$H$_{172}$', 'Si$_{501}$H$_{228}$', 'Si$_{1009}$H$_{412}$'] 
schemes = ['FSSH','mSDM', 'IDA']
isteps = [2500,2500,2500,2500,2500,1]
fsteps = [4500,4500,4500,4500,4500,2001]
dt = 0.5 # fs
# This will be used to find all of the indices within this range above the first excited state.
energy_window = 0.25#0.12 # eV
"""
About the method:
    - 1: It will generate the hot "excess" energy decay vs time data and fit.
    - 2: The same as 1 but it will generate the excess energy using another procedure.
    - 3: Fitting the recovery population dynamics vs time.
"""
# Fit energy vs time data
method = 1 
for c, folder_name in enumerate(folder_names):
    for scheme in schemes:
        title = titles[c]
        os.system(F'find ../../../5_namd/{folder_name} -type d -name "*{scheme}_*" > out.log')
        #os.system(F'find {folder_name}/step4-0.4ev -type d -name "*{scheme}_*" > out.log')
        if scheme=='IDA':
            scheme = 'ID-A'
        file = open('out.log','r')
        lines = file.readlines()
        file.close()
        os.system('rm out.log')
        folders = []
        for i in range(len(lines)):
            folders.append(lines[i].split()[0])
        print('# Found folders:',len(folders))
        print(folders)
        
         # # Reading the energies from Step3
        # Next, we need the energies of the basis we used to perform the NAD which in here is the electron-only excitation basis. Since the energies should be in order, we have to specify the `istep` and `fstep` for each system. The energies are then subtracted from the ground state energy so that we have only the excitation energies.
        
        # In[131]:
        
        
        istep = isteps[c]
        fstep = fsteps[c]
        nstates = np.loadtxt(folders[0]+'/SH_pop.txt').shape[1]
        a = range(nstates)
        energies = []
        for step in range(istep, fstep):
            e_tmp = np.diag( sp.load_npz(F'../../../4_nacs/{folder_name}/res-electron-only/Hvib_sd_{step}_re.npz').todense().real[a,:][:,a] )
            energies.append(e_tmp)
        energies = np.array(energies)
        for i in range(1,energies.shape[1]):
            energies[:,i] -= energies[:,0]
        energies[:,0] -= energies[:,0]
        energies *= units.au2ev
        
        
        # Then, we compute the average energies of each state over the MD trajectory. This will be helpful for fitting the functions to compute the NAD timescales with good accuracies. Infact, what we want to do is to compute the recovery dynamics timescales for a specific energy window which in here is chosen to be `0.1 eV`.
        
        # In[132]:
        
        
        average_energies = np.average(energies,axis=0)
        recovery_states_indices = np.where(np.logical_and(average_energies-average_energies[1]<energy_window,
                                                          average_energies-average_energies[1]>=0))[0]
        if len(recovery_states_indices)<3:
            recovery_states_indices = [1,2,3]
        print('The population of these states will be used to compute the recovery dynamics timescales:',
              recovery_states_indices)
        
        
        # # Where is the population?
        # An intuitive way to answer this question is to plot the population matrix, which shows the population of each state in time. The following plot shows the map of this matrix. Sometimes, the population is distributed over many states and the map does not clearly show where the population is? Therefore, it would be intuitive to plot the population with respect to excited state energies not with respect to excite states indices. Another way is to show that where the first, second, third, and fourth maximum population are at each step. This can be helpful since sometimes some states might have the same population and both are maximum. Here is how we plot it for our systems.
        # 
        # Note that, the colors are indicative of the first, second, third, ... maximum population. 
        
        # In[128]:
        
        
        #get_ipython().run_line_magic('matplotlib', 'notebook')
        fontsize=12
        plt.figure(figsize=(3.21*1.5, 2.41*1.5), dpi=600, edgecolor='black', frameon=True)
        
        colors = ['red','blue','green']
        for i in range(len(folders)):
            folder = folders[i]
            pop_mat = np.loadtxt(folder+'/SH_pop.txt')
            if i==0:
                # The md_time is for plotting and will also be used for fitting the data.
                md_time = np.arange(0,pop_mat.shape[0]*dt, dt)
                # This line of code is not really needed but you can check 
                # which state was the initial state
                #istate_index = np.where(pop_mat[0,:]==1.0)[0][0]
            # We want to plot the map with respect to energies so we need meshgrid
            # This energy mesh is from the first excited state to the maximum excited state
            [X, Y] = np.meshgrid(md_time, np.linspace(np.min(energies[:,1::]),
                                                      np.max(energies[:,:]),nstates-1))
            # The contour energy levels to be plotted (populations) from 0 to 1
            levels = np.arange(0,1.01,0.01)
            plt.contourf(X, Y, pop_mat[:,1::].T, levels=levels, vmax=1.0, cmap='hot')
            # How many maximum to be considered is obtained from the colors
            for k in range(1,len(colors)+1):
                ene_all_step = []
                for j in range(pop_mat.shape[0]):
                    # sort and find the maximum population in the pop_mat
                    sorted_args = pop_mat[j,np.argsort(pop_mat[j,:])]
                    max_index = np.where(pop_mat[j,:]==sorted_args[-k])[0][0]
                    # The excitation energy at which the kth maximum population is
                    ene_tmp = energies[j,max_index]
                    ene_all_step.append(ene_tmp)
                # Plot the max population
                plt.plot(md_time, ene_all_step, '.', markersize=0.1, color=colors[k-1], linewidth=0.1)
        
                
        cb = plt.colorbar()
        cb.ax.tick_params(labelsize=fontsize)
        plt.xlim(0,1000)
        plt.ylim(np.min(energies[:,1::]), np.max(energies[:,:]))
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlabel('Time, fs', fontsize=fontsize)
        plt.ylabel('Excitation energy, eV', fontsize=fontsize)
        plt.title(F'{title} {scheme}', fontsize=fontsize)
        plt.tight_layout()
        plt.savefig(F'{folder_name}_{scheme}_pop_mat.jpg', dpi=600)
        plt.close()
        

        # Plot energy vs time
        fontsize = 10
        plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
        for i in range(energies.shape[1]):
            plt.plot(md_time, energies[:,i], linewidth=1)
        
        plt.xlim(0,100)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.xlabel('Time, fs', fontsize=fontsize)
        plt.ylabel('Excitation energy, eV', fontsize=fontsize)
        plt.title(F'{title} e$^-$', fontsize=fontsize)
        plt.tight_layout()
        plt.savefig(F'{folder_name}_energy_time.jpg', dpi=600)
        plt.close()


        # # Fitting the data
        # Finally, we want to fit the data to Gaussian-exponential function of the form:
        # 
        # $$f(t, A, \tau_1, \tau_2, E_0) = A\exp(-\frac{t}{\tau_1}) + (E_0-A)\exp(-(\frac{t}{\tau_2})^2)$$
        # 
        # Then, the overall timescale is obtained from:
        # 
        # $$\tau = \frac{A}{E_0}\tau_1 + \frac{(E_0-A)}{E_0}\tau_2$$
        # 
        # The error bars are computed with respect to the following formula:
        # 
        # $$\epsilon = Z\frac{s}{\sqrt{N}}$$
        # 
        # where $Z$ is the confidence interval and $s$ is the standard deviation of the data, and $N$ is the number of samples. We use a $Z$ value of $1.96$ for $95\%$ confidence interval.
        
        # In[129]:
        
        
        #get_ipython().run_line_magic('matplotlib', 'notebook')
        
        def gaussian_exponential(t, a, tau1, tau2, E0):
            return a*np.exp( -( t/tau1 ) ) + (E0-a)*np.exp( -np.power(( t/tau2 ),2) )

        #def gaussian_exponential(t, a, tau1, tau2):
        #    return a*np.exp( -( t/tau1 ) ) + (1-a)*np.exp( -np.power(( t/tau2 ),2) )
        
        plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
        
        tau_vals = []
        a_vals = []
        tau1_vals = []
        tau2_vals = []
        E0_vals = []
        e2 = (energies.T-energies[:,1]).T
        counter = 0
        print(F'--------------------------{folder_name}----{scheme}---------------------')
        for i in range(len(folders)):
            folder = folders[i]
            pop_mat = np.loadtxt(folder+'/SH_pop.txt')
            if method==1:
                data_to_fit = np.sum(np.multiply(pop_mat[:,1::], energies[:,1::]), axis=1)
                data_to_fit -= min(data_to_fit)
            if method==2:
                data_to_fit = np.sum(np.multiply(pop_mat[:,1::], e2[:,1::]), axis=1)
            if method==3:
                data_to_fit = 1-np.sum(pop_mat[:,recovery_states_indices],axis=1)
            #pop_recovery = 1-np.sum(pop_mat[:,1:15],axis=1)
            #plt.plot(md_time, pop_recovery, color='gray')
            #popt, pcov = curve_fit( gaussian_exponential, md_time, pop_recovery, 
            #                       bounds=([0.0, 0.0, 0.0, 0.0],[np.inf, np.inf, np.inf, np.inf]))
            #popt, pcov = curve_fit( gaussian_exponential, md_time, data_to_fit,
            #                       bounds=([0.0, 0.0, 0.0],[1, np.inf, np.inf]))
            if method==1 or method==2:
                #popt, pcov = curve_fit( gaussian_exponential, md_time, data_to_fit, 
                #                       bounds=([0.0, 0.0, 0.0, data_to_fit[0]-0.01],[data_to_fit[0],
                #                                                                      np.inf, np.inf,data_to_fit[0]+0.01]))
                popt, pcov = curve_fit( gaussian_exponential, md_time, data_to_fit,
                                       bounds=([0.0, 0.0, 0.0, data_to_fit[0]-0.01],[data_to_fit[0],
                                                                                      np.inf, np.inf,data_to_fit[0]+0.01]))

            if method==3:
                popt, pcov = curve_fit( gaussian_exponential, md_time, data_to_fit,
                                       bounds=([0.0, 0.0, 0.0, 1.00000],[1, np.inf, np.inf, 1.00001]))

            a, tau1, tau2, E0 = popt
            #E0 = 1
            #a, tau1, tau2= popt
            tau = a*tau1/E0+(E0-a)*tau2/E0
            # Computing the R-squared
            residuals  = data_to_fit - gaussian_exponential(md_time, *popt)
            ss_res     = np.sum(residuals**2)
            ss_tot     = np.sum((data_to_fit - np.mean(data_to_fit))**2)
            r_squared  = 1.0 - (ss_res / ss_tot)
            # plot and use the data only if the R^2>0.85
            if counter==0:
                average_fit = np.zeros(gaussian_exponential(md_time, *popt).shape)
            if r_squared>0.9 and tau1<5000:
                average_fit += gaussian_exponential(md_time, *popt)
                plt.plot(md_time, data_to_fit, color='gray')
                print(a, tau1, tau2, E0, r_squared)
                tau_vals.append(tau)
                a_vals.append(a)
                tau1_vals.append(tau1)
                tau2_vals.append(tau2)
                E0_vals.append(E0)
                print(tau)
                counter += 1
        
        
        tau_vals = np.array(tau_vals)
        a_vals = np.array(a_vals)
        tau1_vals = np.array(tau1_vals)
        tau2_vals = np.array(tau2_vals)
        E0_vals = np.array(E0_vals)
        # The confidence interval
        Z = 1.96
        N = counter
        s = np.std(tau_vals)
        tau_ave = np.average(tau_vals)
        error_bar = Z*s/np.sqrt(N)
        s1 = np.std(tau1_vals)
        tau_ave1 = np.average(tau1_vals)
        error_bar1 = Z*s1/np.sqrt(N)
        s2 = np.std(tau2_vals)
        tau_ave2 = np.average(tau2_vals)
        error_bar2 = Z*s2/np.sqrt(N)
        sa = np.std(a_vals)
        a_ave = np.average(a_vals)
        error_bar_a = Z*sa/np.sqrt(N)
        se0 = np.std(E0_vals)
        e0_ave = np.average(E0_vals)
        error_bar_e0 = Z*se0/np.sqrt(N)


        if counter>0:
            plt.plot(md_time, average_fit/counter, color='red')
        
        print(F'The overall timescale for Gaussian-exponential function: {tau_ave}+-{error_bar}, tau1: {tau_ave1}+-{error_bar1}, tau2: {tau_ave2}+-{error_bar2}, a: {a_ave}+-{error_bar_a}, E0: {e0_ave}+-{error_bar_e0}, number of samples: {N}')
        
        plt.xlabel('Time, fs', fontsize=10)
        if method==1 or method==2:
            plt.ylabel('Excess energy, eV', fontsize=10)
        if method==3:
            plt.ylabel('Population', fontsize=10)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.title(F'{title} {scheme}')
        plt.xlim(0,1000)
        #plt.ylim(0,1)
        plt.tight_layout()
        plt.savefig(F'{folder_name}_{scheme}_fit_gaussian_exponential.jpg', dpi=600)
        plt.close()
        
        # In[72]:
        
        
        # np.sum(pop[:,0],axis=1)


        def stretched_compressed(t, tau, beta):
            return np.exp( -np.power(( t/tau ),beta) )
        
        plt.figure(figsize=(3.21, 2.41), dpi=600, edgecolor='black', frameon=True)
        
        tau_vals = []
        beta_vals = []
        e2 = (energies.T-energies[:,1]).T
        counter = 0
        for i in range(len(folders)):
            folder = folders[i]
            pop_mat = np.loadtxt(folder+'/SH_pop.txt')
            if method==1:
                data_to_fit = np.sum(np.multiply(pop_mat[:,1::], energies[:,1::]), axis=1)
                data_to_fit -= min(data_to_fit)
            if method==2:
                data_to_fit = np.sum(np.multiply(pop_mat[:,1::], e2[:,1::]), axis=1)
            if method==3:
                data_to_fit = 1-np.sum(pop_mat[:,recovery_states_indices],axis=1)
            #pop_recovery = 1-np.sum(pop_mat[:,recovery_states_indices],axis=1)
            #pop_recovery = 1-np.sum(pop_mat[:,1:15],axis=1)
            #plt.plot(md_time, data_to_fit, color='gray')
            popt, pcov = curve_fit( stretched_compressed, md_time, data_to_fit, 
                                   bounds=([0.0, 0.0],[np.inf, np.inf]))
            
            tau, beta = popt
            # Computing the R-squared
            residuals  = data_to_fit - stretched_compressed(md_time, *popt)
            ss_res     = np.sum(residuals**2)
            ss_tot     = np.sum((data_to_fit - np.mean(data_to_fit))**2)
            r_squared  = 1.0 - (ss_res / ss_tot)
            # plot and use the data only if the R^2>0.8
            if r_squared>0.9 and tau<50000:
                plt.plot(md_time, data_to_fit, color='gray')

                if counter==0:
                    average_fit = stretched_compressed(md_time, *popt)
                else:
                    average_fit += stretched_compressed(md_time, *popt)
                print(tau, beta, r_squared)
                tau_vals.append(tau)
                beta_vals.append(beta)
                counter += 1
        
        
        tau_vals = np.array(tau_vals)
        beta_vals = np.array(beta_vals)
        
        # The confidence interval
        Z = 1.96
        N = counter
        s = np.std(tau_vals)
        tau_ave = np.average(tau_vals)
        beta_ave = np.average(beta_vals)
        s_b = np.std(beta_vals)
        error_bar = Z*s/np.sqrt(N)
        error_bar_b = Z*s_b/np.sqrt(N)
        if counter>0:
            plt.plot(md_time, stretched_compressed(md_time, tau_ave, beta_ave), color='red')
            plt.plot(md_time, stretched_compressed(md_time, tau_ave+s, beta_ave), color='red', ls='dashed')
            plt.plot(md_time, stretched_compressed(md_time, tau_ave-s, beta_ave), color='red', ls='dashed')
        
        print(F'The overall timescale for stretched-compressed function: {tau_ave}+-{error_bar}, average beta value: {beta_ave}+-{error_bar_b}, number 0f samples: {N}')
        
        plt.xlabel('Time, fs', fontsize=10)
        if method==1 or method==2:
            plt.ylabel('Excess energy, eV', fontsize=10)
        if method==3:
            plt.ylabel('Population', fontsize=10)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.title(F'{title} {scheme}')
        plt.xlim(0,1000)
        #plt.ylim(0,1)
        #plt.yticks([0,0.25,0.5,0.75,1.0])
        plt.tight_layout()
        plt.savefig(F'{folder_name}_{scheme}_fit_stretched_compressed.jpg', dpi=600)
        plt.close()

