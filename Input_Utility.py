#
#       each element in controlling the process 
#



from Input_projectors import Get_label_list_by_Projector


def Gen_Projectors(params_setup):

    New_Pj_shot = params_setup['New_Pj_shot']
        
    # ---------------------------------------------- #
    #       create & save with projectors            #
    # ---------------------------------------------- #
    label_list, T_rec, Ncpu_Pj, saveP_bulk, Partition_Pj = \
        Get_label_list_by_Projector(params_setup, New_Pj_shot)
        
    print('  ######       saveP_bulk = {}  for Ncpu_Pj = {}    ########\n'.format(saveP_bulk, Ncpu_Pj))

    if Partition_Pj == 1 and New_Pj_shot[0] !=2:
        print('  ######   saveP_bulk = {} > 1  -->  combining bulks of Pj_list   ############\n'.format(saveP_bulk))
        New_Pj_shot[0] = 2
        label_list2, T_rec, Ncpu_Pj, saveP_bulk, Partition_Pj = \
            Get_label_list_by_Projector(params_setup, New_Pj_shot)

        if not sorted(label_list) == label_list2:
            print('  *******   ERROR  in combining  projectors   *******\n')
            return

    return label_list, T_rec, Ncpu_Pj, saveP_bulk


def write_f1_Measure(f1, T_rec, mea_method, measure_method, meas_path, Ncpu_meas):

    if mea_method == 0:
        f1.write('         {} -->  meas time = {}  (each label file) \n'.format(meas_path, T_rec['measurement']))

    elif mea_method == 1:           #   bunch of Pauli operators stored together

        if measure_method == 1:
            f1.write('         {} -->  meas time = {}  (direct the whole label_list  | num_cpus = {}) \n'.format(meas_path, \
                                                                        T_rec['measurement'], Ncpu_meas))
        elif measure_method == 3:
            f1.write('         {} -->  meas time = {}  (parallel for each label part | num_cpus = {}) \n'.format(meas_path, \
                                                                        T_rec['measurement'], Ncpu_meas))

