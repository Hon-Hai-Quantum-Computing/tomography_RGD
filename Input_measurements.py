
import time
import sys
sys.path.append('./quQST')
import measurements
import projectors
from methodsRGD import Amea

import numpy as np
import os
import pickle
import multiprocessing

import Get_param

from Utility import Gen_Rho, Gen_randM, state_measure, state_measure_save, state_S_coef, split_list

#from Utility import State_Naming, Gen_Rho, Gen_randM, state_measure, \
#                    state_measure_save, state_measure_wID, state_S_coef, split_list, mp_state_measure


from importlib import reload
reload(measurements)


def create_measurement_bash(Nk, mea, StateName, meas_path, ml_lab_files, Fname):
        #print('  Fname = {}'.format(Fname))

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

    tm1 = time.time()            

    ml_dir_file = os.listdir(meas_path)
    ml_lab_files = [xx for xx in ml_dir_file if xx.startswith("ml_labels_")]
    
    print('  ml_lab_files = {}\n'.format(ml_lab_files))
    print('    Run_method = {}\n'.format(Run_method))

    if Run_method == 0:
        method_Name = 'sequential qiskit measure each label part'

        for file in ml_lab_files:
            print('  file = {}\n'.format(file))

            #  This seems not working for called file, i.e. not in the main
            #sys.argv = ['parallel_measurements.py']
            #sys.argv.append(Nk) 
            #sys.argv.append(mea)
            #sys.argv.append(StateName)
            #sys.argv.append(meas_path)
            #sys.argv.append(file)
            #exec(open('parallel_measurements.py').read())

            params = '{} {} {} {} {} '.format(Nk, mea, StateName, meas_path, file)
            Rec    = '{}/Rec_qiskit_meas_dict_{}.txt'.format(meas_path, file[10:])
            
            cmd    = 'python parallel_measurements.py {} > {}'.format(params, Rec)
            
            print('cmd = {}'.format(cmd))
            os.system(cmd)

        #os.system('wait')

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
    # --------------------------------------------------------------- #
    #      combiner data_dict_each  & calc measurement_list           #
    # --------------------------------------------------------------- #

    td1 = time.time()

    print('   #####   start loading each data_dict from each label_part   ####\n')

    label_file = [xx for xx in os.listdir(meas_path) if xx.startswith('ml_labels_')]
    #print('label_file = {}'.format(label_file))

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

            #for label in labels_each:
            #    data_dict_ALL[label] = data_dict_each[label]
            data_dict_ALL.update(data_dict_each)
            
        #print('   combined -->  data_dict_ALL  = {}'.format(data_dict_ALL))

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

            #print(' labels_each = {}'.format(labels_each))
            #print(' measurement_part = {}\n'.format(measurement_part))

        ID_sort = sorted(range(len(labels_ALL)), key=lambda k: labels_ALL[k])
        
        sort_labels = sorted(labels_ALL)
        measurement_sort   = [measurement_combine[xx] for xx in ID_sort]

        del labels_ALL
        del measurement_combine

        measurement_outALL = (sort_labels, measurement_sort)
        print('\n    >>>>>>>>        measurement_list is combined and sorted         <<<<< \n')

        measurements.MeasurementStore.Save_measurement_list_file(meas_path, sort_labels, measurement_sort, 'ALL')

        #print('  labels_ALL = {}'.format(labels_ALL))
        #print('  sorted(labels_ALL) = {}'.format(sorted(labels_ALL)))
        #print('  ID_sort = {}'.format(ID_sort))
        #print('  measurement_combine = {}'.format(measurement_combine))
        #print('  measurement_sort    = {}'.format(measurement_sort))



    #print('  measurement_outALL  = {}'.format(measurement_outALL))
    #print('  measurement_outALL[0] = {}'.format(measurement_outALL[0]))
    #print('  measurement_outALL[1] = {}'.format(measurement_outALL[1]))

    return measurement_outALL

def split_label_4_measure(label_list):

    num_cpus   = multiprocessing.cpu_count()
    #num_cpus2  = num_cpus - 2
    if num_cpus > len(label_list):
        num_cpus = 3
    else:
        num_cpus = min(48, num_cpus)

    print(' num_cpus = {},  labels size ={}'.format(num_cpus, len(label_list)))

    label_part = split_list(label_list, num_cpus)
    #print('  label_part = {}'.format(label_part))

    return num_cpus, label_part


def parallel_calc_measure_dict(Nk, m, mea, StateName, label_list, meas_path):

    num_cpus, label_part = split_label_4_measure(label_list)

    # ------------------------------------- #
    #   save label_part in meas_path        #
    # ------------------------------------- #
    if not os.path.exists(meas_path):
        os.makedirs(meas_path)

    for ii, labs in enumerate(label_part):
        lab_file = '{}/ml_labels_{}'.format(meas_path, ii)
        #print(' ------------------------------------- ')
        #print(' ii = {}, labs = {}, lab_file ={}'.format(ii, labs, lab_file))

        with open(lab_file, 'wb') as f:
            pickle.dump(labs, f)
        print('       saving {}-th label_part DONE'.format(ii))

    print('\n *****************   saving label_part files DONE  *********** \n\n')

    # --------------------------------------------- #
    #   loading label_part & calc measurement_dict  #
    # --------------------------------------------- #
    #cmd = 'ls {}/ml_labels*'.format(meas_path)
    #os.system(cmd)

    Run_measure = 1
    t_parallel = Run_measurement_dict(Nk, mea, StateName, meas_path, Run_measure)

    #print(' measurement_outALL = {}'.format(measurement_outALL))

    return t_parallel, num_cpus, label_part


def compare_direct_and_exact(Nk, m, mea, StateName, label_list, mea_method, Pj_method, \
                meas_path, proj_path):
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

    ym = Amea(projector_list, target_density_matrix, m, 2**Nk) * np.sqrt(m/2**Nk)
    print(' ym = {}'.format(ym))

def Pure_state_measure(params_setup, label_list, measure_method=3):

    Nk         = params_setup['n']
    m          = params_setup['num_labels']
    mea        = params_setup['num_shots']
    StateName  = params_setup['StateName']
    Dir        = params_setup['Dir']
    #Dir0       = params_setup['Dir0']
    #Nr         = params_setup['Nr']
    #proj_path  = params_setup['projector_store_path']
    #Pj_method  = params_setup['Pj_method']

    meas_path  = params_setup['measurement_store_path']

    mea_method = params_setup['mea_method']


    if len(label_list) >  10000:
        measure_method = 3

    if measure_method == 1:            #    direct calc
        tt1 = time.time()

        input_S, target_density_matrix, rho = Gen_Rho(Nk, StateName)

        measurement_dict, data_dict_list, backend = \
            state_measure(mea, label_list, input_S)
        #print('measurement_dict = {} \n'.format(measurement_dict))

        parallel = 0
        if len(label_list) > 5000:
            parallel = 1
        Get_param.Do_measure(meas_path, measurement_dict, mea_method, parallel)

        # --------------------------------------------------- #
        #   save measurement_dict  in  xxx_qiskit.dat         #
        # --------------------------------------------------- #
        state_measure_save(Dir, params_setup, label_list, 
                measurement_dict, data_dict_list, backend,
                input_S, target_density_matrix, rho)

        tt2 = time.time()
        t_parallel = tt2 - tt1
        num_cpus   = 1

    #elif measure_method == 2:   #   use parallel CPU  (not OK since subprocess call multiprocessing)
    #    measurement_dict2, data_dict_list2, backend = mp_state_measure(params_setup, label_list)

    elif measure_method == 3:          #  each CPU save result

        t_parallel, num_cpus, label_part = parallel_calc_measure_dict(Nk, m, mea, StateName, label_list, meas_path)

        measurement_outALL = combine_data_dict_Calc_Mlist(num_cpus, meas_path)   #  combine data_dict  &  calc measurement_list

    # ------------------------------------------------------------------- #
    #   (comparison) direct calc measurement_dict ->  measurement_list    #
    #                   &    projectors  -->  exact   measurement_list    #
    # ------------------------------------------------------------------- #

    #compare_direct_and_exact(Nk, m, mea, StateName, label_list, mea_method, Pj_method, \
    #                    meas_path, proj_path)

    return t_parallel, num_cpus


def Get_measurement_by_labels(params_setup, label_list, New_Pj_shot, StVer, \
            T_rec, measure_method = 3):
                              
    tm0 = time.time()

    Nk         = params_setup['n']
    m          = params_setup['num_labels']
    mea        = params_setup['num_shots']
    StateName  = params_setup['StateName']
    Nr         = params_setup['Nr']
    Dir0       = params_setup['Dir0']
    Dir        = params_setup['Dir']
    proj_path  = params_setup['projector_store_path']
    meas_path  = params_setup['measurement_store_path']
    Pj_method  = params_setup['Pj_method']
    mea_method = params_setup['mea_method']


    if New_Pj_shot[1] == 1:           #  do qiskit -> do shot measurement
        print('  ********   starting to DO measurement (qiskit or qutip)  ************** \n')

        if StateName == 'rand':

            target_density_matrix, rho = Gen_randM(Nk, StateName, Nr, StVer, Dir0)
            input_S = None
            #print('   Dir0   = {}'.format(Dir0))
            #print('   StVer  = {}'.format(StVer))
            #print('    rho   = {}'.format(rho))

            s_label, yProj, zN, yProj_Exact = \
                state_S_coef(Dir, params_setup, target_density_matrix, rho)
            
            num_cpus = 1

        else:
            t_parallel, num_cpus = Pure_state_measure(params_setup, label_list, measure_method)
            T_rec['parallel_measure'] = t_parallel


    elif New_Pj_shot[1] == 2:           #  to load measurement_dict from  xxx_qiskit.dat

        if measure_method == 3:
            num_cpus, label_part = split_label_4_measure(label_list)

            print('   ********    start combining data_dict  [ Ncpu = {} ] -->  to do measurement_list    *****\n'.format(num_cpus))

            measurement_outALL = combine_data_dict_Calc_Mlist(num_cpus, meas_path)   #  combine data_dict  &  calc measurement_list
            #print('  -------   measurement_outputALL  = {}\n'.format(measurement_outALL))


    elif New_Pj_shot[1] == 0:           #  to load measurement_list from  file

        if StateName == 'rand':

            version  = params_setup['version']
            Noise    = params_setup['Noise']

            ver_meas = version[1]
            Model    = Noise[0]

            Fname1 = '{}/zN{}_v{}_Noise'.format(meas_path, Model, ver_meas)
            Fname2 = '{}/zN{}_v{}_measurements'.format(meas_path, Model, ver_meas)
            print("  Fname1 = {}".format(Fname1))
            print("  Fname2 = {}".format(Fname2))

            with open(Fname1, 'rb') as f:
                yProj_Exact, zN, yProj, Noise, labels, params_setup, target_density_matrix, rho = pickle.load(f)

            #with open(Fname2, 'rb') as f:
            #    labels, yProj, zN = pickle.load(f)
            
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

                #print('measurement_output         = {}'.format(measurement_output))    
                print('len(measurement_output)    = {}'.format(len(measurement_output)))    
                print('len(measurement_output[0]) = {}'.format(len(measurement_output[0])))    
                print('len(measurement_output[1]) = {}'.format(len(measurement_output[1])))    


        num_cpus = 1


    elif New_Pj_shot[1] == -1:           #  to load measurement_dict from  xxx_qiskit.dat
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
            #In_Nk, In_m, In_shot, input_S, measurement_dict, data_list, s_label, params_qiskit, \
            #In_Nk, In_m, In_shot, input_S, measurement_dict, data_dict_list, s_label, params_qiskit, \
            #target_density_matrix, rho = pickle.load(f)
            measurement_dict, data_dict_list, s_label, params_setup, backend, \
            target_density_matrix, input_S, rho = pickle.load(f)

    tm1 = time.time()
    dt = tm1 - tm0
    print('\n  <<<<<<<<               Do_measure time             = {}           >>>>>>\n'.format(dt))
    T_rec['measurement'] =  dt

    #print(" sys.argv content = {}".format(sys.argv))
    #print('     argv length  = {}'.format(len(sys.argv)))
    #print("      sys.argv[0] = {}".format(sys.argv[0]))

    if StateName != 'rand':
        try:
            target_density_matrix.shape
            print('   exists')
        except:
            input_S, target_density_matrix, rho = Gen_Rho(Nk, StateName)  # for returning
            print('  target_density_matrix is created for returning')


    return target_density_matrix, input_S, T_rec, num_cpus

