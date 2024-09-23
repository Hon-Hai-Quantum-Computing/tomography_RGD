
import time
import sys
sys.path.append('./quQST')
import measurements
import projectors

import numpy as np
import os
import pickle
import multiprocessing

import Get_param

from Utility import Gen_Rho, Gen_randM, state_measure, state_measure_save, state_S_coef, split_list
from Utility import Gen_Wst


def create_measurement_bash(Nk, mea, StateName, meas_path, ml_lab_files, Fname):
    """ to crate a bash script (run in Linux) for parallel processing measurements for different Pauli labels

    Args:
        Nk (int): number of qubits for the state
        mea (int): number of sampled Pauli operators
        StateName (str): name of the state 
        meas_path (str): path to store the measurements 
        ml_lab_files (list): list of file names of parallel measurement results  
        Fname (str): the bash script file name to run parallel measurements 
                = 'Run_measure.sh' defined in input_measurement/Run_measurement_dict
    """

    os.system('rm {}'.format(Fname))

    f2 = open(Fname, 'w')
    f2.write('#!/bin/bash\n\n')
    f2.write('echo " Do the measurement by each part"\n\n')
    f2.write('Nk={}\n'.format(Nk))
    f2.write('mea={}\n'.format(mea))
    f2.write('StateName={}\n'.format(StateName))
    f2.write('meas_path={}\n'.format(meas_path))
    f2.write('Rec=$meas_path"/Rec_qiskit_meas_dict_"\n\n')
    
    f2.write('echo "   Nk       = " $Nk\n')
    f2.write('echo "   mea      = " $mea\n')
    f2.write('echo " StateName  = " $StateName\n')
    f2.write('echo " meas_path  = " $meas_path\n')
    f2.write('echo " Rec Header = " $Rec\n\n')

    cmd = 'python parallel_measurements.py $Nk $mea $StateName $meas_path'
    print('  cmd   = {}'.format(cmd))

    for ii, file in enumerate(ml_lab_files):
        ID_file = int(file[10:])
        print(' file = {}  --> ID = {}\n'.format(file, ID_file))

        f2.write('taskset -c {} {} {} > $Rec"{}.txt" &\n'.format(ii, cmd, file, ID_file))
    f2.write('\n')
    f2.write('wait\n\n')
    f2.write('echo "All the measurement_dict DONE"')
    f2.close()

    print('  Run parallel measurement file = {} is updated\n'.format(Fname))


def Run_measurement_dict(Nk, mea, StateName, meas_path, Run_method=1):
    """ to load each separate part of sample Pauli labels
        and then do the actual measurements separately | in parallel

    Args:
        Nk (int): number of qubits of the state
        mea (int): number of shot measurements
        StateName (str): name of the target state
        meas_path (str): path to store the measurement results 
        Run_method (int, optional): method to do the parallel measurement. Defaults to 1.

    Returns:
        float: time of computation
    """

    tm1 = time.time()            

    ml_dir_file = os.listdir(meas_path)
    ml_lab_files = [xx for xx in ml_dir_file if xx.startswith("ml_labels_")]
    
    print('  ml_lab_files = {}\n'.format(ml_lab_files))
    print('    Run_method = {}\n'.format(Run_method))

    if Run_method == 0:
        method_Name = 'sequential qiskit measure each label part'

        for file in ml_lab_files:
            print('  file = {}\n'.format(file))

            params = '{} {} {} {} {} '.format(Nk, mea, StateName, meas_path, file)
            Rec    = '{}/Rec_qiskit_meas_dict_{}.txt'.format(meas_path, file[10:])
            
            cmd    = 'python parallel_measurements.py {} > {} &'.format(params, Rec)

            print('cmd = {}'.format(cmd))
            os.system(cmd)

        os.system('wait')

    elif Run_method == 1:
        method_Name = 'parallel CPU measurement_dict (Run.measure.sh)'

        Fname = 'Run_measure.sh'

        create_measurement_bash(Nk, mea, StateName, meas_path, ml_lab_files, Fname)

        os.system('chmod +x {}'.format(Fname))
        os.system('./{}'.format(Fname))


    tm2 = time.time()

    print('  ---------- qiskit measurement_dict DONE time = {} by  {}'.format(tm2-tm1, method_Name))
    print('\n  ---------- DONE of measurement_dict for ALL label_part  ----------  \n')
    
    return tm2-tm1

def combine_data_dict_Calc_Mlist(num_split, meas_path, Run_method=2):
    """  to read out all stored results from the measurement data dictionary (meas_path)
        & combine them into a list of measurement results 
                    corresponding to coefficients of Pauli operators

    Args:
        num_split (int): number of stored files representing the measurement results, where
            each file corresponds to separate Pauli labels
        meas_path (str): directory path of storing measurements
        Run_method (int, optional): control parameter for combining results. Defaults to 2.

    Returns:
        tuple: (measurement_outALL) = (label_list, measurement_list)
                passed from measurements.MeasurementStore.Save_measurement_by_data_dict function, where
            label_list       = list of Pauli operator labels
            measurement_list = list of coefficients corresponding to the label_list

    """
    # --------------------------------------------------------------- #
    #      combiner data_dict_each  & calc measurement_list           #
    # --------------------------------------------------------------- #

    td1 = time.time()

    print('   #####   start loading each data_dict from each label_part   ####\n')

    label_file = [xx for xx in os.listdir(meas_path) if xx.startswith('ml_labels_')]

    if len(label_file) == num_split:
        print('          num_split  = {}'.format(num_split))
    else:
        print('          ERROR:  the ml_labels_xx  file is not consistent with num_split')
        return


    if Run_method == 1:

        labels_ALL = []
        data_dict_ALL = {}
        for ID_label in range(num_split):

            labels_each, data_dict_each = \
                measurements.MeasurementStore.load_labels_data(meas_path, ID_label)
            labels_ALL    += labels_each

            data_dict_ALL.update(data_dict_each)
            
        if len(labels_ALL) ==  len(data_dict_ALL):
            print('\n   #####     num of data_dict_ALL = {}      #####'.format(len(data_dict_ALL)))
        else:
            print('   ******   ERROR  in combining  data_dict_ALL  \n')
            return
        td2 = time.time()
        print('   #####     DONE of loading and combining data_dict_ALL  -->  time = {}\n'.format(td2-td1))

        method = 0
        measurement_outALL = measurements.MeasurementStore.Save_measurement_by_data_dict( meas_path, \
                sorted(labels_ALL), data_dict_ALL, method, Name='ALL')

    elif Run_method == 2:

        measurement_combine = []
        labels_ALL = []
        for ID_label in range(num_split):
            print('\n +++++++++++        start loading   {}-th data_dict  (to calc measurement_list)    ++++++++++++ '.format(ID_label))

            labels_each, data_dict_each = \
                measurements.MeasurementStore.load_labels_data(meas_path, ID_label)
            labels_ALL    += labels_each

            method = 1
            measurement_part = measurements.MeasurementStore.Save_measurement_by_data_dict( meas_path, \
                labels_each, data_dict_each, method, Name=ID_label, ToSave=0)


            if labels_each != measurement_part[0]:
                print('  the measurement_part is NOT consistent with labels')
                break
            print(' ===========        Done calc  {}-th  measurement_list    =========== '.format(ID_label))

            del labels_each
            del data_dict_each

            measurement_combine  += measurement_part[1]


        ID_sort = sorted(range(len(labels_ALL)), key=lambda k: labels_ALL[k])
        
        sort_labels = sorted(labels_ALL)
        measurement_sort   = [measurement_combine[xx] for xx in ID_sort]

        del labels_ALL
        del measurement_combine

        measurement_outALL = (sort_labels, measurement_sort)
        print('\n    >>>>>>>>        measurement_list is combined and sorted         <<<<< \n')

        measurements.MeasurementStore.Save_measurement_list_file(meas_path, sort_labels, measurement_sort, 'ALL')


    return measurement_outALL

def split_label_4_measure(label_list):
    """ to split the sampled Pauli label list into several chunks for later parallel processing
        according to the number of multiproccessors, 
        which will be set to 3 chunks if number of all labels is less than the number of multiprocessors

    Args:
        label_list (list): list of all labels of sampled Pauli matrices

    Returns:
        _type_: _description_
        int : (num_cpus)   number of processors used for parellel generation of measurement results
        list: (label_part) list of labels sperated by parellel generation
    """

    num_cpus   = multiprocessing.cpu_count()
    #num_cpus2  = num_cpus - 2
    if num_cpus > len(label_list):
        num_cpus = 3
    else:
        num_cpus = min(48, num_cpus)

    print(' num_cpus = {},  labels size ={}'.format(num_cpus, len(label_list)))

    label_part = split_list(label_list, num_cpus)

    return num_cpus, label_part


def parallel_calc_measure_dict(Nk, m, mea, StateName, label_list, meas_path, Run_measure = 1):
    """ to call the parallel computation of measurements

    Args:
        Nk (int): number of qubits of the state
        m (int): number of sampled Pauli operators 
        mea (int): number of shot measurements 
        StateName (str): name of the target state
        label_list (list): list of labels for sampled Pauli operators 
        meas_path (str): directory path to store the measurement results
        Run_measure (int, optional): control parameter for parallel computation. 
                Defaults to 1 (running bash script in the Linux system).

    Returns:
        float: (t_parallel) time to do parallel measurements
        int  : (num_cpus)   number of parallel cpu
        list : (label_part) list of separate Pauli label list for parallel computation
    """

    num_cpus, label_part = split_label_4_measure(label_list)

    # ------------------------------------- #
    #   save label_part in meas_path        #
    # ------------------------------------- #
    if not os.path.exists(meas_path):
        os.makedirs(meas_path)

    for ii, labs in enumerate(label_part):
        lab_file = '{}/ml_labels_{}'.format(meas_path, ii)

        with open(lab_file, 'wb') as f:
            pickle.dump(labs, f)
        print('       saving {}-th label_part DONE'.format(ii))

    print('\n *****************   saving label_part files DONE  *********** \n\n')

    # --------------------------------------------- #
    #   loading label_part & calc measurement_dict  #
    # --------------------------------------------- #

    t_parallel = Run_measurement_dict(Nk, mea, StateName, meas_path, Run_measure)

    return t_parallel, num_cpus, label_part


def compare_direct_and_exact(Nk, m, mea, StateName, label_list, mea_method, Pj_method, \
                meas_path, proj_path):
    """ Do the measurement and convert 

    Args:
        Nk (int): number of qubits
        m (): _description_
        mea (int): number of shot measurements
        StateName (str): state name of target state
        label_list (list): list of sampled Pauli operators
        mea_method (int): the method to save | load measurement_dict (count_dict) 
        Pj_method (int): the method to save | load  projectors
        proj_path (str): the path to store the generated Pauli operators 
    """
    print('  *****   the original whole label_list measurements for comparison   ***** \n')
    input_S, target_density_matrix, rho = Gen_Rho(Nk, StateName)

    measurement_dict, data_dict_list, backend = \
        state_measure(mea, label_list, input_S)

    measurement_output = Get_param.Do_measure(meas_path, measurement_dict, mea_method)
    print('measurement_output = {}'.format(measurement_output))

    # --------------------------------------------- #
    #   to find projectors --> cacl exact ym        #
    # --------------------------------------------- #
    if Pj_method == 0:         
        projector_dict_list = projectors.ProjectorStore.load(proj_path, label_list)
    elif Pj_method == 1:    
        projector_dict_list = projectors.ProjectorStore.load_PoolMap(proj_path, label_list)

    projector_list = [projector_dict_list[label] for label in sorted(label_list)]


def Pure_state_measure(params_setup, label_list, measure_method=3, Run_measure=1):
    """ generation of measurement data for the pure state

    Args:
        params_setup (dict): dictionary of parameters
        label_list (list): list of label names representing the sampled Pauli operators 
        measure_method (int, optional): the method of numerical processing the measurement. Defaults to 3.
                    if = 1: direct processing all the label_list
                    if = 3: multiple cpu processors to parallelly producing the measurement data
        Run_measure (int, optional): control parameter for parallel computation. Defaults to 1.
                    (default 1 for running bash script in Linux)
                    set to 0 if the bash script cannot be run in the (Windows) system
    """

    Nk         = params_setup['n']
    m          = params_setup['num_labels']
    mea        = params_setup['num_shots']
    StateName  = params_setup['StateName']

    meas_path  = params_setup['measurement_store_path']

    mea_method = params_setup['mea_method']


    if len(label_list) >  10000:
        measure_method = 3             #    change to parallel computation

    if measure_method == 1:            #    direct calc
        tt1 = time.time()

        input_S, target_density_matrix, rho = Gen_Rho(Nk, StateName)

        measurement_dict, data_dict_list, backend = \
            state_measure(mea, label_list, input_S)

        parallel = 0
        if len(label_list) > 5000:
            parallel = 1
        Get_param.Do_measure(meas_path, measurement_dict, mea_method, parallel)

        # --------------------------------------------------- #
        #   save measurement_dict  in  xxx_qiskit.dat         #
        # --------------------------------------------------- #

        state_measure_save(meas_path, params_setup, label_list, measurement_dict, 
                input_S, target_density_matrix, rho)

        tt2 = time.time()
        t_parallel = tt2 - tt1
        num_cpus   = 1


    elif measure_method == 3:          #  each CPU save result

        t_parallel, num_cpus, label_part = parallel_calc_measure_dict(Nk, m, mea, StateName, label_list, meas_path, Run_measure)

        measurement_outALL = combine_data_dict_Calc_Mlist(num_cpus, meas_path)   #  combine data_dict  &  calc measurement_list

    # ------------------------------------------------------------------- #
    #   (comparison) direct calc measurement_dict ->  measurement_list    #
    #                   &    projectors  -->  exact   measurement_list    #
    # ------------------------------------------------------------------- #

    return t_parallel, num_cpus


def Get_measurement_by_labels(params_setup, label_list, T_rec):
    """to obtain the measurement results (the coefficients of each Pauli operator)
    from the given label_list corresponding to the sampled Pauli operators 

    Args:
        params_setup (dict): dictionary storing all parameters
        label_list (list): list of sampled Pauli operator labels
        T_rec (dict): dictionary storing the time for generating 
                    Pauli operators and their corresponding measurements

    Returns:
        np.array: (target_density_matrix) the matrix representation of the density matrix of interest
        circuit : (input_S) the circuit to generate the pure state vector of interest (eg) GHZ, Had state
                 will be set to None if the state vector is not generated from the circuit
                        (eg) the random density matrix 

        dict    : (T_rec) dictionary storing the time for generating 
                    Pauli operators and their corresponding measurement results
        int     : (num_cpus) number of processors for parallel calculation of measurement results

    """
                              
    tm0 = time.time()

    Nk         = params_setup['n']
    StateName  = params_setup['StateName']
    DirRho     = params_setup['DirRho']
    meas_path  = params_setup['measurement_store_path']
    mea_method = params_setup['mea_method']

    Gen_New_Meas = params_setup['Gen New Measure']   # control parameter for obtaining measurements
    measure_method = params_setup['measure_method'] 

    if Gen_New_Meas == 1:        # create new measurement
        print('  ********   starting to DO measurement (qiskit or qutip)  ************** \n')

        if StateName == 'rand' or StateName == 'KapRnd':

            target_density_matrix, rho = Gen_randM(params_setup, DirRho)
            input_S = None

            s_label, yProj, zN, yProj_Exact = \
                state_S_coef(params_setup, target_density_matrix, rho)
            
            num_cpus = 1

        elif StateName == 'GHZ' or StateName == 'Had':
            Run_measure = 0  # control parameter for parallel computation (default 1 for running bash script in Linux)
                             # set to 0 if the bash script cannot be run in the (Windows) system
            t_parallel, num_cpus = Pure_state_measure(params_setup, label_list, measure_method, Run_measure)
            T_rec['parallel_measure'] = t_parallel

        elif StateName == 'quGH' or StateName == 'quHa' or StateName == 'quWst':

            Wst, target_density_matrix = Gen_Wst(Nk, DirRho)

            s_label, yProj, zN, yProj_Exact = \
                state_S_coef(params_setup, target_density_matrix)

            input_S = None
            num_cpus = 1

    elif Gen_New_Meas == 0:           #  to load measurement_list from  file

        if StateName == 'rand' or StateName == 'KapRnd':

            version  = params_setup['version']
            Noise    = params_setup['Noise']

            ver_meas = version[1]
            zModel    = Noise

            Fname1 = '{}/zN{}_v{}_Noise'.format(meas_path, zModel, ver_meas)
            Fname2 = '{}/zN{}_v{}_measurements'.format(meas_path, zModel, ver_meas)
            print("  Fname1 = {}".format(Fname1))
            print("  Fname2 = {}".format(Fname2))

            with open(Fname1, 'rb') as f:
                yProj_Exact, zN, yProj, Noise, labels, params_setup, target_density_matrix, rho = pickle.load(f)
            
            input_S = None

            print('  labels[:10]       = {}'.format(labels[:10]))
            print('  yProj[:10]        = {}'.format(yProj[:10]))
            print('  yProj_Exact[:10]  = {}'.format(yProj_Exact[:10]))
            print('  ZN[:10]           = {}'.format(zN[:10]))

        else:
            if mea_method == 1:         #  measurement_list.pickle   or ml_ALL.pickle
                if measure_method == 1:
                    measurement_output = measurements.MeasurementStore.Load_measurement_list(meas_path)

                elif measure_method == 3:
                    measurement_output = measurements.MeasurementStore.Load_measurement_list(meas_path, 'ALL')

                print('len(measurement_output)    = {}'.format(len(measurement_output)))    
                print('len(measurement_output[0]) = {}'.format(len(measurement_output[0])))    
                print('len(measurement_output[1]) = {}'.format(len(measurement_output[1])))    


        num_cpus = 1


    elif Gen_New_Meas == 2:     #  

        if measure_method == 3:
            num_cpus, label_part = split_label_4_measure(label_list)

            print('   ********    start combining data_dict  [ Ncpu = {} ] -->  to do measurement_list    *****\n'.format(num_cpus))

            measurement_outALL = combine_data_dict_Calc_Mlist(num_cpus, meas_path)   #  combine data_dict  &  calc measurement_list

    elif Gen_New_Meas == -1:           #  to load measurement_dict from  xxx_qiskit.dat
        # ---------------------------------------------------------- #
        #   to download measurement data from qiskit recorded data   #
        #           measurement                                      #
        #       &   projectors                                       #
        #   --> not necessary (already got from the above)           #
        # ---------------------------------------------------------- #
        num_cpus = 1

        Fname = '{}_qiskit.dat'.format(Dir)

        print(' to load qiskit measurement data from {}'.format(Fname))
        with open(Fname, 'rb') as f:
            measurement_dict, data_dict_list, s_label, params_setup, backend, \
            target_density_matrix, input_S, rho = pickle.load(f)

    tm1 = time.time()
    dt = tm1 - tm0
    print('\n  <<<<<<<<               Do_measure time             = {}           >>>>>>\n'.format(dt))
    T_rec['measurement'] =  dt
    T_rec['Ncpu_meas']   = num_cpus

    if StateName != 'rand':
        try:
            target_density_matrix.shape
            print('   exists')
        except:
            input_S, target_density_matrix, rho = Gen_Rho(Nk, StateName)  # for returning
            print('  target_density_matrix is created for returning')


    return target_density_matrix, input_S, T_rec, num_cpus
