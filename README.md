# Project_Libra-xTB

This repository contains the files to perform NA-MD simulations in the extended tight-binding (xTB) framework for Si nanocrystals (NCs) of different sizes, a 2D 
graphitic carbon nitride monolayer, and a titanium-based MIL-125-NH2 MOF structure.

More details on how to run the code are brought in [Libra CompChemCyberTraining](https://github.com/compchem-cybertraining/Tutorials_Libra/tree/master/6_dynamics/2_nbra_workflows) repository. Detailed explanations about installation and running the CP2K inputs can be found in [here](https://github.com/compchem-cybertraining/Tutorials_CP2K).


The prepared items are the main iputs that we used for our recent project, "Nonadiabatic Molecular Dynamics with Extended Density Functional Tight-Binding: Application to Nanocrystals and Periodic Solids". Below are the details about how one can use and modify these inputs to generate the results we have obtained.

## 1. Molecular dynamics

This folder contains the CP2K inputs for molecular dynamics along with the trajectory `xyz` files. To uncompress the tar files that contain the trajectories that have been split into multiple parts, you can do this in these folders:

```
cat *.tar.bz2.part* > file.tar.gz.joined
tar xvf file.tar.gz.joined
```

For the the C3N4 monolayer and MIL-125-NH2, the timestep is 1.0 fs and for Si NCs, the trajectories are available for time steps of 0.5 and 0.25 fs. These files are used for computing the MO overlaps and time-overlaps `3_overlaps` folder.

## 2. MO overlaps

The detialed explanation ondifferent inputs of this section are brought in [here](https://github.com/compchem-cybertraining/Tutorials_Libra/tree/master/6_dynamics/2_nbra_workflows/7_step2_cp2k). To compute the MO overlaps for the Si NCs with smaller time step you need to specify the correct path in the `run_template.py` to their corresponding trajectories in `2_molecular dynamics/Si_0.25fs_dt` folder.

## 3. Nonadiabatic couplings

The computation of the NACs requires the path to MO overlaps and time-overlaps. Detailed explanation about the parameters used in the `step3.py` files are available in [thislink](https://github.com/compchem-cybertraining/Tutorials_Libra/tree/master/6_dynamics/2_nbra_workflows/8_step3_cp2k). 

The NACs for the Si NCs are in the electron-only excitation basis. Also, you can compute the NACs between KS states as well. 

For MIL-125-NH2 and C3N4 monolayer, we first obtain the corrected MO overlaps obtained by running `step3.run_step3_ks_nacs_libint(params_ks)`, then, we use them, which are stored in folder `res-ks`, to compute the SD overlaps in the mixed electron and hole excitation basis using `step3.run_step3_sd_nacs_libint(params_sd)`. The reason we did not use this approach for the Si NCs is that the state-tracking algorithm for a large number of states would be prohibitively expensive but it is used for the periodic and smaller systems with low number of states.

## 4. Nonadiabatic molecular dynamics

The inputs prepare here are adapted to the most recent version of Libra 5.2.1. Here, we set the most important parameters that you can use to perform different types of calculations, including longer trajectories by repeating the Hamiltonian matrices. The important parameters in the `run_namd.py` files are as follows:

`path_to_npz_files`: The full path to the NAC folder. 

`istep`: The index of the initial step for reading the files.

`fstep`: The index of the final step for reading the files.

`igeo`: The initial geometry to start the dynamics from the read files from `istep` to `fstep`. For example, if you set `igeo = 10`, it  will start the dynamics from the `istep+10`th step. 

`nfiles`: The number of files that were read. This value is `fstep-istep`.

`active_space`: The active space means that which states you want to involve in the dynamics. For example, to include only the ground state and the first excited state, you can set it to `range(0,2)`. You can include more states in the dynamics as desired.

`user_nsteps`: This parameter is dependent on how long you want to do the dynamics. If this value is larger than `nfiles`, it will use a repeating Hamiltonian approach. So, you can set it to an arbitrary positive integer value. For example, you can run 100 ps dynamics using the 2 ps Hamiltonian files.

`isNBRA`: This keyword, in the `params_nbra`, is used only for FSSH, ID-A, and mSDM type of calculations. If it is set to 1, the Hamiltonian related properties are only computed for one of the trajectories. This is good for nanoscale structures and for dynamics that involve a high number of states and trajectories. For DISH scheme, this value should be set to `0`. Currently, large number of states with a large number of trajectories cannot be used for this type of calcualtions. 

Other parameters are fully explained in details in [this link](https://github.com/compchem-cybertraining/Tutorials_Libra/tree/master/6_dynamics/2_nbra_workflows/9_step4_cp2k).

After setting up your input, you can submit the jobs through slurm for multiple geometries using `sh exe.sh`. If you do not have access to multiple nodes, you can uncomment `# for icond in range(igeo,fgeo):` line and set your geometries as desired (don't forget to uncomment `for icond in [igeo]:`!).

The results will be stored for each initial geometry in `namd_res_*` folders for each initial state, scheme, and batch.


Note that you can apply changes to the Hamiltonian files when reading them. For example, in order to perform the electron-hole recombination dynamics in the spin adapted configuration basis (SAC), you need to multiply the NAC values by `sqrt(2)` which can be done using:

```
hvib_im[:,0] *= np.sqrt(2)
hvib_im[0,:] *= np.sqrt(2)
```

or in Si NCs, you can compute the dynamics by allowing only adjacent states transitions using:

```
# First make a new matrix with the same dimension as the original matrix, hvib_im1
hvib_im = np.zeros(hvib_im1.shape)
for k1 in range(hvib_im1.shape[0]):
    try:
        # only appending the first off-diagonals in the matrix
        hvib_im[k1,k1]    = hvib_im1[k1,k1]
        hvib_im[k1,k1+1]  = hvib_im1[k1,k1+1]
        hvib_im[k1+1,k1]  = hvib_im1[k1+1,k1]
    except:
        pass
# Append hvib_im
hvib_im_MATRIX = data_conv.nparray2MATRIX(hvib_im)

```

## 5. Visualizing

The Python files in folder `6_plot_properties` are a set of diverse files that are used for analyzing and plotting the properties of the generated data from previous steps. We highly recommend to run the `.py` files in `6_plot_properties` folder in Jupyter notebooks. The reason is that the user can better understand the procedure of how the data are found, analyzed, and plotted. The `namd_results.py` files are for fitting the NA-MD results and `prop.py` files are for plotting the electronic structure poperties (such as PDOS and energy vs time plots) and NAC maps in different basis. The user can modify these scripts as desired. One of the most important libraries that we use is `glob` which can be used to find files with specific names. For example, finding the imaginary part of the Hamiltonians, one can use `glob.glob('res-ks/*im*')` and use the found files for plotting the average NAC maps. 


Other files are used for analyzing the PDOS and movement of the core and surface atoms in Si NCs. The coordinates of the core and surface atoms are provided in `6_plot_properties/Si_NCs` folder.



