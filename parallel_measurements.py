#
#   for running the measurements parallelly
#   

import os
import sys
sys.path.append('./quQST')
import measurements


import pickle
from Utility import Gen_Rho, state_measure
import time


def label_to_measurement_dict(Nk, StateName, num_shots, meas_path, file):        
    
    ID_label = int(file[10:])
    lab_file = '{}/{}'.format(meas_path, file)                 

    print('\n -------------- to LOAD label file = {}  -> ID = {}'.format(lab_file, ID_label))

    print('           Nk    = {}'.format(Nk))
    print('       num_shots = {}'.format(num_shots))
    print('       StateName = {}'.format(StateName))
    print('       meas_path = {}'.format(meas_path))
    print('         file    = {}'.format(file))


    with open(lab_file, 'rb') as f:
        labels = pickle.load(f)

    # ------------------------- #
    #   Do measurement_dict     #
    # ------------------------- #
    tm1 = time.time()

    input_S, target_density_matrix, rho = Gen_Rho(Nk, StateName)

    measurement_dict_part, data_dict_list_part, backend = \
        state_measure(num_shots, labels, input_S)

    tm2 = time.time()

    print('             state_measure for {}  -->  time = {}\n'.format(file, tm2-tm1))
    print('\n      -------------  get result (count_dict_list -> saveInOne) from measurement_dict  --------------------- \n')    

    measurement_store = measurements.MeasurementStore(measurement_dict_part)
    measurement_store.saveInOne(meas_path, ID_label)            #  only save count_dict


if __name__ == "__main__":


    Nk        = int(sys.argv[1])
    num_shots = int(sys.argv[2])
    StateName = sys.argv[3]
    meas_path = sys.argv[4]
    file      = sys.argv[5]


    label_to_measurement_dict(Nk, StateName, num_shots, meas_path, file)


 





