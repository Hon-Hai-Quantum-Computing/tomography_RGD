
import time
import numpy as np

import sys
sys.path.append('./quQST')
import projectors

from importlib import reload
reload(projectors)

import Get_param
reload(Get_param)

def compare_proj_dict(label_list, Proj_dict1, Proj_dict2):

    diff = []
    for label in label_list:
        xx = Proj_dict1[label].matrix - Proj_dict2[label].matrix

        diff += [len(xx.data)]

    if np.sum(diff) > 0:
        print('   ERROR in producing and combining Projectors')
    else:
        print('   CHECK OK for combining proejectors')


def Check_combine_projectors(proj_path, label_list, Pj_combine):
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

def Get_label_list_by_Projector(Nk, m, New_Pj_shot, proj_path, Pj_method, T_rec):
    # ------------------------------------------------------------------ #
    #           projectors  ==> generate / loading  label_list           #
    # ------------------------------------------------------------------ #

    tp0 = time.time()

    saveP_bulk   = 0              #  the default saved bulk of Pj_list is 0
    Partition_Pj = 0              #  default:  no partition in Projectors

    if New_Pj_shot[0] == 1:     #  generate new projectors
        tp0 = time.time()

        print(' -------------     Generate new labels    ------------')
        label_list = projectors.generate_random_label_list(m, Nk)
        #print(np.array(label_list))

        # --------------------------------- #
        #   generate & store projectors     #
        # --------------------------------- #
        num_cpus, saveP_bulk, Partition_Pj = Get_param.Do_projector(label_list, proj_path, Pj_method)

    elif New_Pj_shot[0] == 2:     #  loading bulk projectors and combine  -->  save in One File

        saveP_bulk, label_list, num_cpus = combine_proj_SaveInOne(proj_path)

    elif New_Pj_shot[0] == 21:     #  loading bulk projectors and combine & compare

        saveP_bulk, label_list, Pj_combine, num_cpus = combine_proj_SaveInOne(proj_path, 1)
        
        proj_combine_dict = projectors.ProjectorStore.load_PoolMap(proj_path, label_list)

        num_cpus, saveP_bulk = Check_combine_projectors(proj_path, label_list, proj_combine_dict)
        num_cpus, saveP_bulk = Check_combine_projectors(proj_path, label_list, Pj_combine)
        compare_proj_dict(label_list, proj_combine_dict, Pj_combine)


    elif New_Pj_shot[0] == 0:     #  loading projector from proj_path
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
    
    T_rec['proj'] = tp1 - tp0

    return label_list, T_rec, num_cpus, saveP_bulk, Partition_Pj




