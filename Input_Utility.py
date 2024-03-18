#
#       each element in controlling the process 
#



from Input_projectors import Get_label_list_by_Projector


#def Gen_Projectors(params_setup, params_control):
def Gen_Projectors(params_setup):

    #Nk         = params_setup['n']
    #m          = params_setup['num_labels']
    #proj_path  = params_setup['projector_store_path']
    #Pj_method  = params_setup['Pj_method']
    #StateName  = params_setup['StateName']

    #New_Pj_shot = params_control['New_Pj_shot']
    New_Pj_shot = params_setup['New_Pj_shot']
        
    # ---------------------------------------------- #
    #       create & save with projectors            #
    # ---------------------------------------------- #
    label_list, T_rec, Ncpu_Pj, saveP_bulk, Partition_Pj = \
        Get_label_list_by_Projector(params_setup, New_Pj_shot)
        #Get_label_list_by_Projector(Nk, m, New_Pj_shot, proj_path, Pj_method)
        
    print('  ######       saveP_bulk = {}  for Ncpu_Pj = {}    ########\n'.format(saveP_bulk, Ncpu_Pj))

    if Partition_Pj == 1 and New_Pj_shot[0] !=2:
        print('  ######   saveP_bulk = {} > 1  -->  combining bulks of Pj_list   ############\n'.format(saveP_bulk))
        New_Pj_shot[0] = 2
        label_list2, T_rec, Ncpu_Pj, saveP_bulk, Partition_Pj = \
            Get_label_list_by_Projector(params_setup, New_Pj_shot)
            #Get_label_list_by_Projector(Nk, m, New_Pj_shot, proj_path, Pj_method)

        #label_sorted, Pj_combine, bulk_Pj, Ncpu_Pj_Bk = projectors.ProjectorStore.combine_proj_bulk(proj_path)


        #print(' label_list     = {}'.format(label_list))
        #print(' label_list2    = {}'.format(label_list2))

        if not sorted(label_list) == label_list2:
            print('  *******   ERROR  in combining  projectors   *******\n')
            return

    #print(' params_control = {}'.format(params_control))
    #print(' Ncpu_Pj        = {}'.format(Ncpu_Pj))
    
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

