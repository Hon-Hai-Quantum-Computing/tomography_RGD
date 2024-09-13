#
#       each element in controlling the process 
#



from Input_projectors import Get_label_list_by_Projector


def Gen_Projectors(params_setup):
    """ To produce the Pauli Projectors

    Args:
        params_setup (dict): dictionary of parameters specifying the system

    Returns:
        list: (label_list) the list of labels representing the sample Pauli operators
        dict: (T_rec) recording running time 
        int: (Ncpu_Pj) number of processors for parallel producing Pauli projectors
        int: (saveP_bulk) number of final saved projector files for later processing
    """

    #New_Pj_shot = params_setup['New_Pj_shot']
    New_Pj = params_setup['Gen New Proj sampling']

    # ---------------------------------------------- #
    #       create & save with projectors            #
    # ---------------------------------------------- #
    label_list, T_rec, Ncpu_Pj, saveP_bulk, Partition_Pj = \
        Get_label_list_by_Projector(params_setup, New_Pj)
        
    print('  ######       saveP_bulk = {}  for Ncpu_Pj = {}    ########\n'.format(saveP_bulk, Ncpu_Pj))

    if Partition_Pj == 1 and New_Pj !=2:
    #if Partition_Pj == 1 and New_Pj_shot[0] !=2:
        print('  ######   saveP_bulk = {} > 1  -->  combining bulks of Pj_list   ############\n'.format(saveP_bulk))
        New_Pj = 2
        label_list2, T_rec, Ncpu_Pj, saveP_bulk, Partition_Pj = \
            Get_label_list_by_Projector(params_setup, New_Pj)

        if not sorted(label_list) == label_list2:
            print('  *******   ERROR  in combining  projectors   *******\n')
            return

    return label_list, T_rec, Ncpu_Pj, saveP_bulk


def write_f1_Measure(f1, T_rec, mea_method, measure_method, meas_path, Ncpu_meas):
    """ to write some information (about running time) of the measurement to file

    Args:
        f1 (file): file to write meassage
        T_rec (dict): to record the running time 
        mea_method (int):  the method to save | load  measurement_dict (count_dict) 
        measure_method (int): = 1: direct label_list,  = 3: parallel cpu
        meas_path (str): the path to the measurement results
        Ncpu_meas (int): number of processors for parallelly running measurements
    """

    if mea_method == 0:
        f1.write('         {} -->  meas time = {}  (each label file) \n'.format(meas_path, T_rec['measurement']))

    elif mea_method == 1:           #   bunch of Pauli operators stored together

        if measure_method == 1:
            f1.write('         {} -->  meas time = {}  (direct the whole label_list  | num_cpus = {}) \n'.format(meas_path, \
                                                                        T_rec['measurement'], Ncpu_meas))
        elif measure_method == 3:
            f1.write('         {} -->  meas time = {}  (parallel for each label part | num_cpus = {}) \n'.format(meas_path, \
                                                                        T_rec['measurement'], Ncpu_meas))

