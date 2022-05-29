import os
import sys
import time
from liblibra_core import *
from libra_py.workflows.nbra import step3


params_ks = {
          'lowest_orbital': 512-20, 'highest_orbital': 512+21, 'num_occ_states': 10, 'num_unocc_states': 10,
          'use_multiprocessing': True,'nprocs': 12, 'time_step': 1.0, 'es_software': 'cp2k',
          'path_to_npz_files': os.getcwd()+'/../../3_overlaps/2x2_C3N4/res',
          'logfile_directory': os.getcwd()+'/../../3_overlaps/2x2_C3N4/all_logfiles',
          'path_to_save_ks_Hvibs': os.getcwd()+'/res-ks',
          'path_to_save_sd_Hvibs': os.getcwd()+'/res-ks',
          'start_time': 900, 'finish_time': 3999,
          'apply_phase_correction': True, 'apply_orthonormalization': True,
          'do_state_reordering': 2, 'state_reordering_alpha':0
         }

#### For KS states - Applying correction to KS overlaps and computing the NACs in KS space
step3.run_step3_ks_nacs_libint(params_ks)

#### For excited states - Computing the excited states SDs and their overlaps and NACs
# First, for mixed basis sorted by their energy
params_sd = {
          'lowest_orbital': 512-params_ks['num_occ_states']+1, 'highest_orbital': 512+params_ks['num_occ_states'], 'num_occ_states': 5, 'num_unocc_states': 5,
          'isUKS': 0, 'number_of_states': 10, 'tolerance': 0.01, 'verbosity': 0, 
          'use_multiprocessing': True, 'nprocs': 12, 
          'is_many_body': False, 'time_step': 1.0, 'es_software': 'cp2k',
          'path_to_npz_files': os.getcwd()+'/res-ks',
          'logfile_directory': os.getcwd()+'/../../3_overlaps/2x2_C3N4/all_logfiles',
          'path_to_save_sd_Hvibs': os.getcwd()+'/res-mixed-basis-energy',
          'outdir': os.getcwd()+'/res-mixed-basis-energy',
          'start_time': 900, 'finish_time': 3999, 'sorting_type': 'energy',
          'apply_phase_correction': True, 'apply_orthonormalization': True,
          'do_state_reordering': 2, 'state_reordering_alpha':0
         }


step3.run_step3_sd_nacs_libint(params_sd)

#### For excited states - Computing the excited states SDs and their overlaps and NACs
# First, for mixed basis sorted by their identity
params_sd = {
          'lowest_orbital': 512-params_ks['num_occ_states']+1, 'highest_orbital': 512+params_ks['num_occ_states'], 'num_occ_states': 5, 'num_unocc_states': 5,
          'isUKS': 0, 'number_of_states': 10, 'tolerance': 0.01, 'verbosity': 0,
          'use_multiprocessing': True, 'nprocs': 12,
          'is_many_body': False, 'time_step': 1.0, 'es_software': 'cp2k',
          'path_to_npz_files': os.getcwd()+'/res-ks',
          'logfile_directory': os.getcwd()+'/../../3_overlaps/2x2_C3N4/all_logfiles',
          'path_to_save_sd_Hvibs': os.getcwd()+'/res-mixed-basis-identity',
          'outdir': os.getcwd()+'/res-mixed-basis-identity',
          'start_time': 900, 'finish_time': 3999, 'sorting_type': 'identity',
          'apply_phase_correction': True, 'apply_orthonormalization': True,
          'do_state_reordering': 2, 'state_reordering_alpha':0
         }


step3.run_step3_sd_nacs_libint(params_sd)

