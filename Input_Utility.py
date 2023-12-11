#
#       each element in controlling the process 
#



from Input_projectors import Get_label_list_by_Projector
from Input_measurements import Get_measurement_by_labels

from Utility import State_Naming


def Gen_Projectors(params_setup, params_control, T_rec):

    Nk         = params_setup['n']
    m          = params_setup['num_labels']
    proj_path  = params_setup['projector_store_path']
    Pj_method  = params_setup['Pj_method']

    New_Pj_shot = params_control['New_Pj_shot']

    # ---------------------------------------------- #
    #       create & save with projectors            #
    # ---------------------------------------------- #
    label_list, T_rec, Ncpu_Pj, saveP_bulk, Partition_Pj = \
        Get_label_list_by_Projector(Nk, m, New_Pj_shot, proj_path, Pj_method, T_rec)
    
    print('  ######       saveP_bulk = {}  for Ncpu_Pj = {}    ########\n'.format(saveP_bulk, Ncpu_Pj))

    if Partition_Pj == 1 and New_Pj_shot[0] !=2:
        print('  ######   saveP_bulk = {} > 1  -->  combining bulks of Pj_list   ############\n'.format(saveP_bulk))
        New_Pj_shot[0] = 2
        label_list2, T_rec, Ncpu_Pj, saveP_bulk, Partition_Pj = \
            Get_label_list_by_Projector(Nk, m, New_Pj_shot, proj_path, Pj_method, T_rec)

        #label_sorted, Pj_combine, bulk_Pj, Ncpu_Pj_Bk = projectors.ProjectorStore.combine_proj_bulk(proj_path)


        #print(' label_list     = {}'.format(label_list))
        #print(' label_list2    = {}'.format(label_list2))

        if not sorted(label_list) == label_list2:
            print('  *******   ERROR  in combining  projectors   *******\n')
            return

    #print(' params_control = {}'.format(params_control))
    #print(' Ncpu_Pj        = {}'.format(Ncpu_Pj))
    
    return label_list, T_rec, Ncpu_Pj, saveP_bulk




def Data_Gen_proj(Frec, params_control, params_setup):
    
    Pj_method   = params_control['Pj_method']
    mea_method  = params_control['mea_method']    
    New_Pj_shot = params_control['New_Pj_shot'] 

    Nk             = params_setup['n']
    m              = params_setup['num_labels']
    mea            = params_setup['num_shots']

    measure_method = params_setup['measure_method']

    Nr          = params_setup['Nr']
    StateName   = params_setup['StateName']
    Noise       = params_setup['Noise']
    StVer       = params_setup['StVer']

    version     = params_setup['version']
    proj_path   = params_setup['projector_store_path']

    f1 = open(Frec, 'a')
    # ----------------------------- #
    #  first generate  projectors   #
    # ----------------------------- #

    T_proj = {}         #  projector time
    version_new  = []     #  the real generated version
    StVer_new    = []     #  the real StVer


    # ----------------------------------------------------- #
    #       record some information about the simulation    #
    # ----------------------------------------------------- #
    T_rec = {}
    T_rec['Nk'] = Nk
    T_rec['m']  = m
    
    # --------------------------------------------------------------------- #
    #   running only single version, i.e.  only one choice of projectors    #
    # --------------------------------------------------------------------- #

    f1.write('\n ----------------   (projection)   ------------------ \n')       
            
    #version = [ps, 0, Noise]
    #exec(open('Input_param_Loop_v5.py').read())
    #Dir, Name, params_setup, version, proj_path, meas_path, Dir0, StVer = \
    #        State_Naming(Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method)

    print('params_setup = {}'.format(params_setup))


    label_list, T_rec, Ncpu_Pj, saveP_bulk  = Gen_Projectors(params_setup, params_control, T_rec)   #  creating projectors


    print(' ######    calculated version = {}  for New_Pj_shot = {}  #######'.format(version, New_Pj_shot))
    print(' ######       -->  T_rec = {}'.format(T_rec))

    T_proj[proj_path] = T_rec['proj']         # = tp1 - tp0

    version_new.append(version)
    StVer_new.append(StVer)

    if StateName == 'rand':
        print('   rand -->  StVer= {}'.format(StVer))

    if Pj_method == 0:
        f1.write('{}  -->  proj time = {}  (original each Pauli calc | Ncpu = {})\n'.format(proj_path, T_rec['proj'], Ncpu_Pj))
    elif Pj_method == 1:
        f1.write('{}  -->  proj time = {}  (parallel cpu: mpPool_map | Ncpu = {})\n'.format(proj_path, T_rec['proj'], Ncpu_Pj))

    f1.close()
    f1 = open(Frec, 'a')

    f1.close()

    return T_rec, T_proj, params_setup, version_new, StVer_new, label_list



def Data_Gen_measure(Frec, params_control, params_setup, label_list, version_new, StVer_new, T_rec):

    Pj_method   = params_control['Pj_method']
    mea_method  = params_control['mea_method']    

    New_Pj_shot = params_control['New_Pj_shot'] 
    New_Pj_shot[0] = 0                                  #  no genertating projectors

    Nk             = params_setup['n']
    m              = params_setup['num_labels']
    mea            = params_setup['num_shots']

    measure_method = params_setup['measure_method']

    Nr             = params_setup['Nr']

    StateName      = params_setup['StateName']
    StVer          = params_setup['StVer']


    T_meas = {}         #  measurement time
    version_new2 = []
    StVer_new2   = []

    #for ii, version in enumerate(version_new):
    for StVer, version in zip(StVer_new, version_new):
        #print('   StVer   = {}'.format(StVer))
        #print('   version = {}'.format(version))

        f1 = open(Frec, 'a')
        f1.write('\n          -------------   (measurement)  -------------------- \n')

        #if StateName == 'rand':
        #    StVer = StVer_new[ii]
        #    StVer = [StVer[0], 0]              #  Not Generating, but only loading    


        # ------------------------- #
        #   from the second shot    #
        # ------------------------- #

        for shot in range(1):            #  shot better start from 1  (just naming)         

            #exec(open('Input_param_Loop_v5.py').read())


            Dir, Name, params_setup, version, proj_path, meas_path, Dir0, StVer = \
                State_Naming(Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method)


            target_density_matrix, input_S, T_rec, Ncpu_meas = Get_measurement_by_labels(params_setup, label_list, New_Pj_shot, \
                                StVer, T_rec, measure_method)


            print('  New_Pj_shot = {}, version = {}, shot = {}'.format(New_Pj_shot, version, shot))

            T_meas[meas_path] = T_rec['measurement']  # = tm1 - tm0


            #version_new.append(version)
            #StVer_new.append(StVer)
            version_new2.append(version)
            StVer_new2.append(StVer)

            write_f1_Measure(f1, T_rec, mea_method, measure_method, meas_path, Ncpu_meas)


        f1.close()        

    return version_new2, StVer_new2



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

