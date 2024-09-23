
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
    """ Given the shot measurement results and calculate the coefficients for each Pauli operator

    Args:
        measurement_store_path (str): path to store the measurements 
        measurement_dict (_type_): 
        mea_method (int): method to store / load the measurements
        parallel (int, optional): To do the parallel measurements or not. Defaults to 1.
        Name (str, optional): name of the target state. Defaults to None.

    Returns:
        if mea_method == 0:
            return None
        elif mea_method == 1:        
            return labs, measurement_list
            where labs             = list of sampled Pauli operators
                  measurement_list = list of calculated coefficients for each Pauli operator
    """

    print('\n  -------------  get result (saveInOne & Save_measurement_list) from measurement_dict  --------------------- \n')    
    print('        -->  Name = {}'.format(Name))

    # ----------------------------------------------------- #
    #   saving measurement results into pickle              #
    #   i.e. storing  self.count_dict for each label        #
    # ----------------------------------------------------- #
    #print(' --------------   measurements  ------------------- ')

    tm1 = time.time()

    measurement_store = measurements.MeasurementStore(measurement_dict)
        
    if not os.path.exists(measurement_store_path):
        os.makedirs(measurement_store_path)

    if mea_method == 0:         #  original method: to save each Pauli measurement separately
        measurement_store.save(measurement_store_path)
    elif mea_method == 1:       #  the new  method: to save ALL Pauli in ONE file
        tt1 = time.time()

        measurement_store.saveInOne(measurement_store_path, Name)            #  only save count_dict
        tt2 = time.time()

        labs, measurement_list = measurement_store.Save_measurement_list(measurement_store_path, parallel, Name)     #  save calculated measurement_list

        tt3 = time.time()

        print('  *****   measurement_store.saveInOne              -->  time = {}'.format(tt2-tt1))
        print('  *****   measurement_store.Save_measurement_list  -->  time = {}'.format(tt3-tt2))


    tm2 = time.time()
    print('\n          The running time for measurement = {}'.format(tm2 - tm1))

    if mea_method == 0:
        return
    elif mea_method == 1:        
        return labs, measurement_list

def Do_projector(s_label, projector_store_path, Pj_method):
    """ to create Pauli projectors according to given s_label

    Args:
        s_label (list): list of sampled Pauli operator labels
        projector_store_path (str): path to store the sampled Pauli projectors
        Pj_method (int): specify the method to save | load  projectors

    Returns:
        _type_: _description_
        int: (num_cpus) number of processors for producing projectors parallelly
         saveP_bulk, Partition_Pj
    """

    ## ---------------------------------------- ##
    ##  Projectors: Pauli string operators      ##
    ## ---------------------------------------- ##
    #print(' -------------  projectors  ----------------------')

    tp1 = time.time()

    projector_store = projectors.ProjectorStore(s_label)
    
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

    return num_cpus, saveP_bulk, Partition_Pj

def Initial_params(params_setup, target_density_matrix, input_S):
    """ To prepare dictionary of parameters needed for the tomography optimization

    Args:
        params_setup (dict): dictionary of parameters
        target_density_matrix (ndarray): the target density matrix
        input_S (class instance): the target state generated from the circuit

    Returns:
        dict: (params_dict) the dictionary of parameters needed for the optimization        
    """

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

    if Noise == -1:                 ## the usual shot measurement

        params_dict = {'measurement_store_path': measurement_store_path,
                   'projector_store_path': projector_store_path,
                   'DirRho': params_setup['DirRho'],
                   'StateName': params_setup['StateName'],
                   'Proj version': params_setup['Proj version'],
                   'measure version': params_setup['measure version'],
                   'Pj_method': params_setup['Pj_method'],
                   'mea_method': params_setup['mea_method'],   
                   'measure_method': params_setup['measure_method'],                
                   'num_iterations': 150,
                   'n': params_setup['n'],
                   'num_labels': params_setup['num_labels'],
                   'Nr': params_setup['Nr'],
                   'num_shots': params_setup['num_shots'],
                   'Noise': 'shot',        # must change corresponding methods_GetParam.py
                   'convergence_check_period': 1,
                   'target_DM': target_density_matrix}

    elif Noise == 0:               ##  Exact coef:  no noise
        params_dict = {'measurement_store_path': measurement_store_path,
                   'projector_store_path': projector_store_path,
                   'DirRho': params_setup['DirRho'],
                   'StateName': params_setup['StateName'],
                   'Proj version': params_setup['Proj version'],
                   'measure version': params_setup['measure version'],
                   'Pj_method': params_setup['Pj_method'],
                   'mea_method': params_setup['mea_method'],
                   'measure_method': params_setup['measure_method'],
                   'num_iterations': 200,
                   'n': params_setup['n'],
                   'num_labels': params_setup['num_labels'],
                   'Nr': params_setup['Nr'],
                   'Noise': 'Exact',
                   'convergence_check_period': 1,
                   'target_DM': target_density_matrix}

    if input_S:
        print(' ********************************************************* \n')
        print('   there exists pure state   ')
        params_dict['target_state'] = input_S.get_state_vector()            

    if params_setup['StateName'] != 'KapRnd':
        params_dict['Dir2m'] = params_setup['Dir2m']

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
    measurement_store_path = '{}-m{}-shot{}-v{}'.format(m_store, m, mea, version)
    if not os.path.exists(measurement_store_path):
        os.makedirs(measurement_store_path)

    measurement_store = measurements.MeasurementStore(measurement_dict)
    measurement_store.save(measurement_store_path)

    tm2 = time.time()
    print('\n          The running time for measurement = {}'.format(tm2 - tm1))

    ## ---------------------------------------- ##
    ##  Projectors: Pauli string operators      ##
    ## ---------------------------------------- ##
    projector_store_path = '{}-m{}-shot{}-v{}'.format(p_store, m, mea, version)

    # ---------------------------------- #
    #  worker   --> defined in methods   #
    #   to get param_dict                #
    # ---------------------------------- #

    ##  State reconstruction using MiFGD in qutomo  ##
    ##
    params_dict = {'measurement_store_path': measurement_store_path,
                   'projector_store_path': projector_store_path,
                   'num_iterations': 1,
                   'n': params_qiskit['n'],
                   'num_labels': params_qiskit['num_labels'],
                   'backend': params_qiskit['backend'],
                   'num_shots': params_qiskit['num_shots'],
                   'convergence_check_period': 1,
                   'target_state': input_S.get_state_vector(),
                   'target_DM': input_S.get_state_matrix()}

    return params_dict



