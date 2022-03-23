import os
import sys
import time
from liblibra_core import *
from libra_py.workflows.nbra import step3


params = {
          'lowest_orbital': 296-100, 'highest_orbital': 296+100, 'num_occ_states': 100, 'num_unocc_states': 100,
          'use_multiprocessing': True,'nprocs': 12, 'time_step': 1.0, 'es_software': 'cp2k',
          'path_to_npz_files': os.getcwd()+'/../../3_overlaps/Si123H100/res',
          'logfile_directory': os.getcwd()+'/../../3_overlaps/Si123H100/all_logfiles',
          'path_to_save_ks_Hvibs': os.getcwd()+'/res-ks',
          'path_to_save_sd_Hvibs': os.getcwd()+'/res-sd',
          'start_time': 2000, 'finish_time': 5999,
          'apply_phase_correction': True, 'apply_orthonormalization': True,
          'do_state_reordering': 2, 'state_reordering_alpha':0
         }

#### For KS states - Applying correction to KS overlaps and computing the NACs in KS space
step3.run_step3_ks_nacs_libint(params)

#### For excited states - Computing the excited states SDs and their overlaps and NACs
# First, for electron-only basis
params = {
          'lowest_orbital': 296-100, 'highest_orbital': 296+100, 'num_occ_states': 1, 'num_unocc_states': 100,
          'isUKS': 0, 'number_of_states': 10, 'tolerance': 0.01, 'verbosity': 0, 
          'use_multiprocessing': True, 'nprocs': 12, 
          'is_many_body': False, 'time_step': 1.0, 'es_software': 'cp2k',
          'path_to_npz_files': os.getcwd()+'/../../3_overlaps/Si123H100/res',
          'logfile_directory': os.getcwd()+'/../../3_overlaps/Si123H100/all_logfiles',
          'path_to_save_sd_Hvibs': os.getcwd()+'/res-mixed-basis',
          'outdir': os.getcwd()+'/res-mixed-basis',
          'start_time': 2000, 'finish_time': 5999, 'sorting_type': 'energy',
          'apply_phase_correction': True, 'apply_orthonormalization': True,
          'do_state_reordering': 2, 'state_reordering_alpha':0
         }


step3.run_step3_sd_nacs_libint(params)


