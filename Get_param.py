
import time
import os
import shutil

import sys
sys.path.append('./quQST')
import measurements
import projectors

from importlib import reload
reload(measurements)


def Do_measure(measurement_store_path, measurement_dict, mea_method, parallel=1, Name=None):

    print('\n  -------------  get result (saveInOne & Save_measurement_list) from measurement_dict  --------------------- \n')    
    print('        -->  Name = {}'.format(Name))

    # ----------------------------------------------------- #
    #   saving measurement results into pickle              #
    #   i.e. storing  self.count_dict for each label        #
    # ----------------------------------------------------- #
    #print(' --------------   measurements  ------------------- ')

    tm1 = time.time()

    measurement_store = measurements.MeasurementStore(measurement_dict)

    #print('  *****  Do_measure -->  measurement_store = {}'.format(measurement_store))
    #measurement_store_path = '{}_measurement'.format(Dir)

    #if os.path.exists(measurement_store_path):
    #    shutil.rmtree(measurement_store_path)
        #os.unlink(measurement_store_path)
        
    if not os.path.exists(measurement_store_path):
        os.makedirs(measurement_store_path)

    if mea_method == 0:         #  original method: to save each Pauli measurement separately
        measurement_store.save(measurement_store_path)
    elif mea_method == 1:       #  the new  method: to save ALL Pauli in ONE file
        tt1 = time.time()

        #measurement_store.saveInOne(measurement_store_path)                 #  only save count_dict
        measurement_store.saveInOne(measurement_store_path, Name)            #  only save count_dict
        tt2 = time.time()

        #print('measurement_store.label_list = {}'.format(measurement_store.label_list))
        #print('measurement_store.data_dict  = {}'.format(measurement_store.data_dict))
        
        #measurement_list = measurement_store.Save_measurement_list(measurement_store_path)     #  save calculated measurement_list
        #measurement_list = measurement_store.Save_measurement_list(measurement_store_path, 1)     #  save calculated measurement_list

        labs, measurement_list = measurement_store.Save_measurement_list(measurement_store_path, parallel, Name)     #  save calculated measurement_list

        #print('  ****  Do_measure -->  measurement_store = {}'.format(measurement_store))


        tt3 = time.time()

        print('  *****   measurement_store.saveInOne              -->  time = {}'.format(tt2-tt1))
        print('  *****   measurement_store.Save_measurement_list  -->  time = {}'.format(tt3-tt2))


    tm2 = time.time()
    print('\n          The running time for measurement = {}'.format(tm2 - tm1))

    if mea_method == 0:
        return
    elif mea_method == 1:        
        #return measurement_store_path
        return labs, measurement_list

def Do_projector(s_label, projector_store_path, Pj_method):

    ## ---------------------------------------- ##
    ##  Projectors: Pauli string operators      ##
    ## ---------------------------------------- ##
    #print(' -------------  projectors  ----------------------')

    tp1 = time.time()

    projector_store = projectors.ProjectorStore(s_label)

    #projector_store_path = '{}_projectors'.format(Dir)
    #if os.path.exists(projector_store_path):
    #    shutil.rmtree(projector_store_path)
    
    if not os.path.exists(projector_store_path):
        os.makedirs(projector_store_path)

    if Pj_method == 0:         #  original method: saving each projectors separately
        projector_store.populate(projector_store_path)
        num_cpus     = 1
        saveP_bulk   = 0
        Partition_Pj = 0

    elif Pj_method == 1:       #  new method:  saving into one file
        num_cpus, saveP_bulk, Partition_Pj = projector_store.mpPool_map(projector_store_path)

    tp2 = time.time()
    print('\n *****        num_cpus = {}, saved # bulk Pj = {}        ****'.format(num_cpus, saveP_bulk))
    print('\n *****   The running time for projector = {}  **** \n'.format(tp2 - tp1))

    #return projector_store_path
    return num_cpus, saveP_bulk, Partition_Pj

def Initial_params(Dir, params_setup, target_density_matrix):
#def Initial_params(Dir, params_qiskit, input_S):

    #print('\n ------------------   to get parameters    ------------------- ')

    #measurement_store_path = '{}_measurement'.format(Dir)
    #projector_store_path = '{}_projectors'.format(Dir)

    measurement_store_path = params_setup['measurement_store_path']
    projector_store_path = params_setup['projector_store_path']

    if not os.path.exists(measurement_store_path):
        print(measurement_store_path)
        print('\n **** measurement_store_path  NOT  exist  ****  \n')
        return 
    
    if not os.path.exists(projector_store_path):
        print('\n **** projector_store_path  NOT  exist  ****  \n')
        return 

    ##  State reconstruction using MiFGD in qutomo  ##
    ##
    Noise = params_setup['Noise']

    if type(Noise) == int:          ## the usual shot measurement

        params_dict = {'measurement_store_path': measurement_store_path,
                   'projector_store_path': projector_store_path,
                   'StateName': params_setup['StateName'],
                   'version': params_setup['version'],
                   'Pj_method': params_setup['Pj_method'],
                   'mea_method': params_setup['mea_method'],   
                   'measure_method': params_setup['measure_method'],                
                   #'num_iterations': 1000,
                   'num_iterations': 150,
                   'n': params_setup['n'],
                   'num_labels': params_setup['num_labels'],
                   'Nr': params_setup['Nr'],
                   #'backend': params_setup['backend'],
                   'num_shots': params_setup['num_shots'],
                   'convergence_check_period': 1,
                   #'target_state': input_S.get_state_vector(),
                   #'target_DM': input_S.get_state_matrix()
                   'target_DM': target_density_matrix}
        
    elif type(Noise) == list:      ##  add in the noise model

        params_dict = {'measurement_store_path': measurement_store_path,
                   'projector_store_path': projector_store_path,
                   'StateName': params_setup['StateName'],
                   'version': params_setup['version'],
                   'Pj_method': params_setup['Pj_method'],
                   'mea_method': params_setup['mea_method'],
                   'measure_method': params_setup['measure_method'],
                   #'num_iterations': 1000,
                   'num_iterations': 250,
                   'n': params_setup['n'],
                   'num_labels': params_setup['num_labels'],
                   'Nr': params_setup['Nr'],
                   #'backend': params_setup['backend'],
                   #'num_shots': params_setup['num_shots'],
                   'Noise': Noise,
                   'convergence_check_period': 1,
                   #'target_state': input_S.get_state_vector(),
                   #'target_DM': input_S.get_state_matrix()
                   'target_DM': target_density_matrix}

    return params_dict

# ----------------- #
#   below old       #
# ----------------- #

def Initial_params_v0(m_store, p_store, m, mea, version, measurement_dict, params_qiskit, input_S):

    #m = params_qiskit['num_labels']
    #mea = params_qiskit['num_shots']
    print('\n ------------------   to get parameters    ------------------- ')

    # ----------------------------------------------------- #
    #   saving measurement results into pickle              #
    #   i.e. storing  self.count_dict for each label        #
    # ----------------------------------------------------- #
    #print(' --------------   measurements  ------------------- ')

    tm1 = time.time()
    #measurement_store_path = '{}-m{}-shot{}-v{}-MiFGD-mu{:.3}'.format(m_store, m, mea, version, mu)
    measurement_store_path = '{}-m{}-shot{}-v{}'.format(m_store, m, mea, version)
    if not os.path.exists(measurement_store_path):
        os.makedirs(measurement_store_path)

    measurement_store = measurements.MeasurementStore(measurement_dict)
    measurement_store.save(measurement_store_path)
    #print(dir(measurement_store))

    tm2 = time.time()
    print('\n          The running time for measurement = {}'.format(tm2 - tm1))

    ## ---------------------------------------- ##
    ##  Projectors: Pauli string operators      ##
    ## ---------------------------------------- ##
    #print(' -------------  projectors  ----------------------')
    #projector_store_path = '{}-m{}-shot{}-v{}-MiFGD-mu{:.3}'.format(p_store, m, mea, version, mu)
    projector_store_path = '{}-m{}-shot{}-v{}'.format(p_store, m, mea, version)

    # ---------------------------------- #
    #  worker   --> defined in methods   #
    #   to get param_dict                #
    # ---------------------------------- #

    ##  State reconstruction using MiFGD in qutomo  ##
    ##
    params_dict = {'measurement_store_path': measurement_store_path,
                   'projector_store_path': projector_store_path,
                   #'num_iterations': 1000,
                   'num_iterations': 1,
                   'n': params_qiskit['n'],
                   'num_labels': params_qiskit['num_labels'],
                   'backend': params_qiskit['backend'],
                   'num_shots': params_qiskit['num_shots'],
                   'convergence_check_period': 1,
                   'target_state': input_S.get_state_vector(),
                   'target_DM': input_S.get_state_matrix()}

    return params_dict



