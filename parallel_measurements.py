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
    #print(" sys.argv content = {}".format(sys.argv))
    #print('     argv length  = {}'.format(len(sys.argv)))
    #print("      sys.argv[0] = {}\n".format(sys.argv[0]))

    print('           Nk    = {}'.format(Nk))
    print('       num_shots = {}'.format(num_shots))
    print('       StateName = {}'.format(StateName))
    print('       meas_path = {}'.format(meas_path))
    print('         file    = {}'.format(file))


    with open(lab_file, 'rb') as f:
        labels = pickle.load(f)
    #print('\n        labels   = {} \n'.format(labels))

    # ------------------------- #
    #   Do measurement_dict     #
    # ------------------------- #
    tm1 = time.time()

    input_S, target_density_matrix, rho = Gen_Rho(Nk, StateName)

    measurement_dict_part, data_dict_list_part, backend = \
        state_measure(num_shots, labels, input_S)
    #print('      measurement_dict = {}'.format(measurement_dict_part))

    #state_measure_save(Dir, params_setup, label_list, 
    #        measurement_dict_part, data_dict_list_part, backend,
    #        input_S, target_density_matrix, rho, ID_label)

    tm2 = time.time()

    print('             state_measure for {}  -->  time = {}\n'.format(file, tm2-tm1))
    print('\n      -------------  get result (count_dict_list -> saveInOne) from measurement_dict  --------------------- \n')    

    measurement_store = measurements.MeasurementStore(measurement_dict_part)
    measurement_store.saveInOne(meas_path, ID_label)            #  only save count_dict


if __name__ == "__main__":

    #sys.argv = ['parallel_measurements.py']
    #sys.argv.append(3) 
    #sys.argv.append(500)
    #sys.argv.append('GHZ')
    #sys.argv.append('data/GHZ-3/GHZ-3_m10_s1/GHZ-3_m10_s1_shot500_v1_Measure')
    #sys.argv.append('ml_labels_2')


    Nk        = int(sys.argv[1])
    num_shots = int(sys.argv[2])
    StateName = sys.argv[3]
    meas_path = sys.argv[4]
    file      = sys.argv[5]

    #print('   type(Nk)         = {}'.format(type(Nk)))
    #print('   type(num_shots)  = {}'.format(type(num_shots)))
    #print('   type(StateName)  = {}'.format(type(StateName)))
    #print('   type(meas_path)  = {}'.format(type(meas_path)))


    label_to_measurement_dict(Nk, StateName, num_shots, meas_path, file)


 





