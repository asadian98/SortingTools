from kilosort import run_kilosort, DEFAULT_SETTINGS
from kilosort.io import load_probe

import numpy as np

import sys

def runKS4(path_dat_file, save_path, path_probe_file, n_chan_bin, datFile):
    
    settings = DEFAULT_SETTINGS

    # General setting
    settings['save_preprocessed_copy'] = False
    settings['clear_cache'] = True
    settings['do_CAR'] = True
    settings['invert_sign'] = False
    settings['data_dir'] = path_dat_file
    settings['filename'] = path_dat_file + datFile
    settings['data_file_path'] = path_dat_file + datFile
    settings['results_dir'] = save_path
    settings['data_dtype'] = 'int16'
    settings['fs'] = 30000.0
    settings['nblocks'] = 0 # 0 means no drift correction
    settings['Th_universal'] = 9.0
    settings['Th_learned'] = 8.0
    settings['n_chan_bin'] = int(n_chan_bin)

    p = load_probe(path_probe_file)

    ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate, kept_spikes = \
        run_kilosort(settings=settings, probe=p)


if __name__ == "__main__":

    path_dat_file = sys.argv[1]  # First argument from the command line
    save_path = sys.argv[2]
    path_probe_file = sys.argv[3]
    n_chan_bin = sys.argv[4]
    datFile = sys.argv[5]
    runKS4(path_dat_file, save_path, path_probe_file, n_chan_bin, datFile)





