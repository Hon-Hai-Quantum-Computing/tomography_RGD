
import time
import numpy as np

import sys
sys.path.append('./quQST')
import projectors

import Get_param

from Gen_localPauli import Gen_Site_2s_Pauli

def compare_proj_dict(label_list, Proj_dict1, Proj_dict2):
    """to compare the projector matrices of all labels 
        i.e. compare source 1 and source 2 of instances of ProjectorStore
        (class ProjectorStore is a ditionary storing 
                        all projectors from class Projector)
        
    Args:
        label_list (list): list of all sampled projector labels
        Proj_dict1 (class ProjectorStore): source 1 of instance of ProjectorStore
        Proj_dict2 (class ProjectorStore): source 2 of instance of ProjectorStore
    """

    diff = []
    for label in label_list:
        xx = Proj_dict1[label].matrix - Proj_dict2[label].matrix

        diff += [len(xx.data)]

    if np.sum(diff) > 0:
        print('   ERROR in producing and combining Projectors')
    else:
        print('   CHECK OK for combining proejectors')


def Check_combine_projectors(proj_path, label_list, Pj_combine):
    """wrapper to compare/check the input Pj_combine (all projectors) the same as 
        all the stored proejectors contained in the pro_path (projector directory)

    Args:
        proj_path (str): directory path that stores the proejectors
        label_list (list): list of label names of all proejectors
        Pj_combine (class): an instance of class ProjectorStore that represents
                            all generated projectors

    Returns:
        int : (num_cpus) number of processors for paralle generation of projectors
        int : (saveP_bulk) number of chucks that projectors are stored in the directory,
              i.e. all the proejctors in stored in several chunks (files), where the number of projects
                 stored in each chunk (file) is limited to some maximum number.  
    """
    # --------------------------------------------- #
    #   for comparsion with Pj_method = 0  method   #
    # --------------------------------------------- #

    Pj_method = 0
    num_cpus, saveP_bulk = Get_param.Do_projector(label_list, proj_path, Pj_method)

    #labels = projectors.ProjectorStore.load_labels(proj_path)

    #projector_dict = projectors.ProjectorStore.load(proj_path, label_list)
    projector_dict = projectors.ProjectorStore.load(proj_path)

    compare_proj_dict(label_list, projector_dict, Pj_combine)
    
    return num_cpus, saveP_bulk

def combine_proj_SaveInOne(proj_path, If_Return=0):
    """already produces several projectors(sampling martices) 
       if produced by parallel processes, then need to combine them together

    Args:
        proj_path (str): the directory storing the projectors
        If_Return (int, optional): to return a value for later comparison. Defaults to 0.

    Returns:
        int:  (bulk_Pj)  number of projector chunks (each chunk contains some number of projectors)  
        list: (label_sorted) list of sorted labels of projectors
        int:  (num_cpus) number of cpu processors used to parallellly produce the projectors
        class: (Pj_combine) an instance of class ProjectorStore (all projectors)
        
        if If_Return == 0:
            return bulk_Pj, label_sorted, num_cpus
        elif If_Return == 1:
            return bulk_Pj, label_sorted, Pj_combine, num_cpus

    """

    print('\n   --------   start loading ALL Pj_list_xx.pickle to have Pj_combine  ------ \n')
    # --------------------------------------------- #
    #   combine bulk projectors  &  save            #
    # --------------------------------------------- #
    tc1 = time.time()

    method_combine = 0
    label_sorted, Pj_combine, bulk_Pj, num_cpus = projectors.ProjectorStore.combine_proj_bulk(proj_path, method_combine)
    print('      using bulk_Pj = {} --> to get Pj_combine\n'.format(bulk_Pj))

    tc2 = time.time()
    print('      DONE of combining ALL Pj_list_xx.pickle  into Pj_combine     -->   Time = {}\n'.format(tc2-tc1))

    projectors.ProjectorStore.Save_label_Pj_dict(Pj_combine, label_sorted, proj_path)
    tc3 = time.time()
    print('      Saving file for Pj_combine & sorted(label_combine)           -->   Time = {}\n'.format(tc3-tc2))

    if If_Return == 0:
        return bulk_Pj, label_sorted, num_cpus
    elif If_Return == 1:
        return bulk_Pj, label_sorted, Pj_combine, num_cpus

def Get_label_list_by_Projector(params_setup, New_Pj):
    """the wrapper function for producing/loading the projectors from
        projectors/ProjectorStore

    Args:
        params_setup (_type_): _description_
        New_Pj (int): basically = params_setup['Gen New Proj sampling']
                 = -1   -> specify fixed Pauli list, i.e. onsite and nearest Pauli operator
                 =  0   -> load existing sampling (will check the existence of the directory)
                 =  1   -> generate new sampling, i.e. to create new projectors/Pauli matrices 
    Returns:
        list: (label_list) list of Pauli operator labels
        dict: (T_rec) dictionary recording time for generating all Pauli operators (matrices)
                T_rec['Nk'] = Nk = qubit number
                T_rec['m']  = m  = number of sampled Pauli matrices
                T_rec['proj'] = time for generating all Pauli matrices
        int : (num_cpus) number of processors used for parallel generation of Pauli matrices
        int : (saveP_bulk) number of projector chunks (each chunk contains some number of projectors)  
                i.e. all projectors are generated and stored in some number of chunks,
                     where each chunk contains some number of projectors
        int : (Partition_Pj) = 0: no partition in generating matrices; =1: generate matrices by more than one chunk

    """

    Nk         = params_setup['n']
    m          = params_setup['num_labels']
    proj_path  = params_setup['projector_store_path']
    Pj_method  = params_setup['Pj_method']
    #Obtain_Pj  = params_setup['Obtain_Pj']

    # ------------------------------------------------------------------ #
    #           projectors  ==> generate / loading  label_list           #
    # ------------------------------------------------------------------ #

    tp0 = time.time()

    saveP_bulk   = 0              #  the default saved bulk of Pj_list is 0
    Partition_Pj = 0              #  default:  no partition in Projectors

    if New_Pj == -1:            #  from specification of Pauli list
        tp0 = time.time()

        print(' -------  Generate labels from single & nearest sites    ------------')

        ALL_sitePauli, ALL_2sPauli, symb_P9, numB2s = Gen_Site_2s_Pauli(Nk)
        label_list = ALL_sitePauli + ALL_2sPauli
        print(' len(label_list) = {}'.format(len(label_list)))

        print('  symb_P9       = {}'.format(symb_P9))
        print('  ALL_sitePauli = {}'.format(ALL_sitePauli))
        print('  ALL_2sPauli   = {}'.format(ALL_2sPauli))
        print("     numB2s     = {}".format(numB2s))

        # --------------------------------- #
        #   generate & store projectors     #
        # --------------------------------- #
        num_cpus, saveP_bulk, Partition_Pj = Get_param.Do_projector(label_list, proj_path, Pj_method)


    elif New_Pj == 1:       #  generate new projectors
        tp0 = time.time()

        print(' -------------     Generate new labels    ------------')
        
        #  -----  from random sampling of Pauli list  ----------- #
        label_list = projectors.generate_random_label_list(m, Nk)

        # --------------------------------- #
        #   generate & store projectors     #
        # --------------------------------- #
        num_cpus, saveP_bulk, Partition_Pj = Get_param.Do_projector(label_list, proj_path, Pj_method)

    elif New_Pj == 2:     #  loading bulk projectors and combine  -->  save in One File

        saveP_bulk, label_list, num_cpus = combine_proj_SaveInOne(proj_path)

    elif New_Pj == 21:     #  loading bulk projectors and combine & compare

        saveP_bulk, label_list, Pj_combine, num_cpus = combine_proj_SaveInOne(proj_path, 1)
        
        proj_combine_dict = projectors.ProjectorStore.load_PoolMap(proj_path, label_list)

        num_cpus, saveP_bulk = Check_combine_projectors(proj_path, label_list, proj_combine_dict)
        num_cpus, saveP_bulk = Check_combine_projectors(proj_path, label_list, Pj_combine)
        compare_proj_dict(label_list, proj_combine_dict, Pj_combine)


    elif New_Pj == 0:     #  loading projector from proj_path
        # ------------------------- #
        #   loading projectors      #
        # ------------------------- #
        if Pj_method == 0:      #  original method
            label_list = [fname.split('.')[0] for fname in os.listdir(proj_path)]
        elif Pj_method == 1:    #  new method
            label_list = projectors.ProjectorStore.load_labels_from_file(proj_path)
        #print('label_list = {}'.format(label_list))

        #proj_combine_dict = projectors.ProjectorStore.load_PoolMap(proj_path, label_list)
        num_cpus = 1

    tp1 = time.time()

    print('    label_list[:10]   = {}\n'.format(label_list[:10]))
    print('    num of projectors = len(label_list) = {}'.format(len(label_list)))
    print('    num of saved bulk Proj              = {}\n\n'.format(saveP_bulk))
    print('  <<<<<<<<                Do_projector time    =  {}           >>>>>>\n'.format(tp1 - tp0))
    
    # ----------------------------------------------------- #
    #       record some information about the simulation    #
    # ----------------------------------------------------- #

    T_rec = {}          
    T_rec['Nk'] = Nk
    T_rec['m']  = m
    T_rec['proj'] = tp1 - tp0

    return label_list, T_rec, num_cpus, saveP_bulk, Partition_Pj

