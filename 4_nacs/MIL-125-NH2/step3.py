import os
import sys
import time
from liblibra_core import *
from libra_py.workflows.nbra import step3


params_ks = {
          'lowest_orbital': 504-100, 'highest_orbital': 504+100, 'num_occ_states': 5, 'num_unocc_states': 5,
          'isUKS': 0, 'number_of_states': 10, 'tolerance': 0.01, 'verbosity': 0,
          'use_multiprocessing': True, 'nprocs': 12,
          'is_many_body': False, 'time_step': 1.0, 'es_software': 'cp2k',
          'path_to_npz_files': os.getcwd()+'/../../3_overlaps/MIL-125-NH2/res',
          'logfile_directory': os.getcwd()+'/../../3_overlaps/MIL-125-NH2/all_logfiles',
          'path_to_save_sd_Hvibs': os.getcwd()+'/res-ks',
          'outdir': os.getcwd()+'/res-ks',
          'start_time': 900, 'finish_time': 3997, 'sorting_type': 'energy',
          'apply_phase_correction': True, 'apply_orthonormalization': True,
          'do_state_reordering': 2, 'state_reordering_alpha':0
         }
#### For KS states - Applying correction to KS overlaps and computing the NACs in KS space
step3.run_step3_ks_nacs_libint(params_ks)

#### For excited states - Computing the excited states SDs and their overlaps and NACs
# Mixed basis sorted based on their energy
params_sd = {
          'lowest_orbital': 504-params_ks['num_occ_states']+1, 'highest_orbital': 504+params_ks['num_unocc_states'], 'num_occ_states': 5, 'num_unocc_states': 5,
          'isUKS': 0, 'number_of_states': 10, 'tolerance': 0.01, 'verbosity': 0,
          'use_multiprocessing': True, 'nprocs': 12,
          'is_many_body': False, 'time_step': 1.0, 'es_software': 'cp2k',
          #'path_to_npz_files': os.getcwd()+'/../step2/res-ks',
          'path_to_npz_files': os.getcwd()+'/res-ks',
          'logfile_directory': os.getcwd()+'/../../3_overlaps/MIL-125-NH2/all_logfiles',
          'path_to_save_sd_Hvibs': os.getcwd()+'/res-mixed-basis-energy',
          'outdir': os.getcwd()+'/res-mixed-basis-energy',
          'start_time': 900, 'finish_time': 3997, 'sorting_type': 'energy',
          'apply_phase_correction': True, 'apply_orthonormalization': True,
          'do_state_reordering': 2, 'state_reordering_alpha':0
         }


step3.run_step3_sd_nacs_libint(params_sd)

# Mixed basis sorted based on their identity
params_sd = {
          'lowest_orbital': 504-params_ks['num_occ_states']+1, 'highest_orbital': 504+params_ks['num_unocc_states'], 'num_occ_states': 5, 'num_unocc_states': 5,
          'isUKS': 0, 'number_of_states': 10, 'tolerance': 0.01, 'verbosity': 0,
          'use_multiprocessing': True, 'nprocs': 12,
          'is_many_body': False, 'time_step': 1.0, 'es_software': 'cp2k',
          #'path_to_npz_files': os.getcwd()+'/../step2/res-ks',
          'path_to_npz_files': os.getcwd()+'/res-ks',
          'logfile_directory': os.getcwd()+'/../../3_overlaps/MIL-125-NH2/all_logfiles',
          'path_to_save_sd_Hvibs': os.getcwd()+'/res-mixed-basis-identity',
          'outdir': os.getcwd()+'/res-mixed-basis-identity',
          'start_time': 900, 'finish_time': 3997, 'sorting_type': 'identity',
          'apply_phase_correction': True, 'apply_orthonormalization': True,
          'do_state_reordering': 2, 'state_reordering_alpha':0
         }


step3.run_step3_sd_nacs_libint(params_sd)


