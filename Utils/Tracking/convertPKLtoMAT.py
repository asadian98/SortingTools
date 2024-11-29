# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pickle
import scipy

def convertPKLtoMAT(pklDirectory, ses_name): 
    # The path to your pickle file
    pickle_file_path = pklDirectory + ses_name + '_timestamps' + '.pkl'
    # Open the file in binary read mode
    with open(pickle_file_path, 'rb') as file:
        # Load the data from the file
        frames_timestamps = pickle.load(file)
    
    # Save the dictionary to a .mat file
    scipy.io.savemat(pklDirectory + ses_name + '_timestamps' + '.mat', {'frames_timestamps': frames_timestamps['DEV_000F315D15FA']})
    
    