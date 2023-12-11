
import time
import pickle                   # only for loading qiskit

import Get_param


from importlib import reload
#from Utility import State_Naming, Gen_Rho, Gen_randM, state_measure, \
#                    state_measure_save, state_measure_wID, state_S_coef, split_list, mp_state_measure

from Utility import State_Naming
                    
from Input_Utility import Gen_Projectors 

# ----------------------------------- #
import sys
sys.path.append('./quQST')

import projectors
#import measurements
#from methodsRGD import Amea
reload(projectors)


import numpy as np
import os
#import multiprocessing
#from itertools import repeat
#from functools import reduce

import Input_measurements
reload(Input_measurements)

from Input_measurements import Get_measurement_by_labels

import Input_setting
reload(Input_setting)

from Input_setting import basic_sys_setting




if __name__ == "__main__":

    # ----------------------------------------------- #
    #   (A-0)  parameters for states (density matrix) #  
    #                                                 #
    #   this part can be manually specified           #
    #   just convenience by calling function          #
    #       Nk  =  number of qubits                   #
    #       m   =  number of sampled Pauli matrices   #
    #       mea =  number of shot measurements        #
    # ----------------------------------------------- #

    Choose = 2
    Nk, m, mea = basic_sys_setting(Choose)

    # --------------------------------------------------------------------------------- #
    #   (A-1)  choosing the input state                                                 #
    #                                                                                   #
    #          Nr = number of rank                                                      #
    #                                                                                   #
    #       Noise =  0           for pure states                                        #
    #               [0, 0, 0.1]  for mixed states                                       #
    #                  i.e. [Gaussian model, mean, variance] --> not implemented yet    #
    #                       now only exact value is used (by qutip)                     #
    #                                                                                   #
    #       StVer    = version control for random mixed states                          #
    #       StVer[0] = version of generated random density matrices of rank Nr          #
    #       StVer[1] =  0    loading random density matrices (generated earlier)        #
    #                           --> will check if the directory exists or not           #
    #                   1    generate new random density matrices of rank Nr            # 
    #                           if file already exists for specified version StVer[0]   #
    #                           then the version StVer[0] will +1 automatically         #
    # --------------------------------------------------------------------------------- #

    StateName = 'GHZ'        #  'GHZ' | 'Had' |  'rand'

    if StateName == 'rand':
        Nr    = 3
        Noise = [0, 0, 0.1]  #  actually implementing exact result only
                             #     --> can extend to Gaussian noise or other noise model
        StVer = [1, 1]       #  [version R?, New or Generate]   for (random) State
    else:
        Nr    = 1
        Noise = 0
        StVer = 0            #  by default = 0


    # ----------------------------------------------------------------------------- #
    #   (B) loading or generating sampled Pauli matrices & measurements             #
    #          projectors = labels = Pauli matrices (obs = observables)             #
    #          measurement: using qiskit shot measurement for pure states           #
    #                         or  qutip direct exact calculation for mixed states   # 
    #          Note the measurements for pure states & mixed states are different   #
    #                       for the treatment now                                   #                                 
    #                                                                               #
    #   (B-1)  program control parameters for generating Proj/measurement           #
    #                                                                               #
    #       Pj_method  = 0  (each Pauli matrix is stored separately)                #
    #                    1  (all Pauli matrices saved in chunks)                    #
    #                        if: not too many --> saved in one Pj_List.pickle file  #
    #                        else: will save in several chunks of Pj_list_xx.pickle #
    #                                                                               #
    #       mea_method = 0  (measurement for each Pauli obs is saved separately)    #
    #                    1  (measurement result is stored in chunks)                #
    #       measure_method  = 1   directly calculate the whole label_list           #
    #                         3   parallel measurement for each label part          #
    #                                                                               #
    # ----------------------------------------------------------------------------- #

    Pj_method = 1               #   the method to save | load  projectors
    mea_method = 1              #   the method to save | load  measurement_dict (count_dict) 

    measure_method = 3          #  = 1: direct label_list,  = 3: parallel cpu

    # ------------------------------------------------------------------------------ #
    #   (B-2) version control:  loading or creating data                             #
    #                              of  sampled Pauli matrices & measurement          #
    #                                                                                #
    #       Data_In =  []               defult                                       #
    #                  [Dir, PjFile]    if specified (only for testing)              #
    #                                                                                #
    #    version[0] =  version of sampled Projectors                                 #
    #           [1] =  version of measurement                                        # 
    #           [2] =  specification of Noise                                        #
    #       --> each version will have separate directory                            #
    #             if loading Proj/measurement  --> will check if directory exists    #
    #             if generating new Pj/measure --> will create new directory         #   
    #                                                                                #
    #    New_Pj_shot[0] = 0  loading sampled projectors                              #
    #                     1  creating new sampled projectors                         #
    #                          if the version already exists                         #
    #                          then the version New_Pj_shot[0] will +1 automatically #
    #                     2  combining each Pj_list_xx.pickle into single file first #
    #                                                                                #
    #    New_Pj_shot[1] = 0  loading measurement result from file                    #
    #                     1  doing measurement (qiskit shot or qutip calculation)    #
    #                          if the version already exists                         #
    #                          then the version New_Pj_shot[1] will +1 automatically #
    #                     2  converting shot measured data_dict into the coef needed #
    #                       i.e. get measurement_list = the coef of each Pauli obs   #
    #    (default usage)                                                             #
    #       New_Pj_shot = [1, 1]  = creating new Projector & measurement             #
    # ------------------------------------------------------------------------------ #
    Obtain_Pj = 0
    if Obtain_Pj == 1:

        Dir    = 'data/rand-6-r6/R1/Lk_rand-6-r3_R1_m1638_s1'
        PjFile = 'rand-6-r3_m1638_s1_Proj'

        Data_In = [Dir, PjFile]

        New_Pj_shot = [0, 1]  #  [New Proj?, New shot?]  | [0,0] loading | 

    else:
        Data_In = []

        New_Pj_shot = [1, 1]  #  [New Proj?, New shot?]  | [0,0] loading | 
        #New_Pj_shot = [0, 1]  #  [New Proj?, New shot?]  | [0,0] loading | 
        #New_Pj_shot = [2, 1]  #  [New Proj?, New shot?]  | [0,0] loading | 
        #New_Pj_shot = [2, 0]  #  [New Proj?, New shot?]  | [0,0] loading | 
        #New_Pj_shot = [0, 2]  #  [New Proj?, New shot?]  | [0,0] loading | 
        #New_Pj_shot = [0, 0]  #  [New Proj?, New shot?]  | [0,0] loading | 

    version = [1, 1, Noise]


    # ----------------------------------------------------------------- #
    #   (C) to call each optimization method to calculate               #
    #       if chosen to run,                                           #
    #          each optimization method will be executed sequentially   #
    #               with recording the result & run time                #
    #       MiFGD & RGD can be called separately                        #
    # ----------------------------------------------------------------- #

    Call_MiFGD = 1      #  = 1: call the MiFGD optimization to calculate
    Call_RGD = 1        #  = 1: call the RGD   optimization to calculate
    
    To_save_X0 = 0      #  to save the initial X0 or not? 

    # -------------------------------------------------------------------------------- #
    #   (C-1) parameters for calling MiFGD optimization algorithm                      #
    #                                                                                  #
    #     InitX_MiFGD = 0   generate random matrix according to MiFGD paper & code     #
    #                                                                                  #
    #     [Usage]   worker = methodsMiFGD.BasicWorker(params_dict, input_S)            #
    #               worker.compute(InitX_MiFGD)                                        #
    #     demonstrated by executing  MiFGD_optRun.py                                   #
    # -------------------------------------------------------------------------------- #

    InitX_MiFGD = 0;     #Ld_MiFGD = None

    muList = [3/4, 0.5, 1/3, 0.2, 0.25, 0.125]
    Num_mu = 1      #  number of mu for running MiFGD
                    #  (eg) = 2 --> run cases of muList = [3/4, 0.5]

    # ---------------------------------------------------------------------------------------------- #
    #   (C-2) parameters for calling RGD optimization algorithm                                      #
    #                                                                                                #
    #       InitX_RGD = 0   generate random X0  (same implementation as MiFGD)                       #
    #                   1   generate X0 according to paper = Hr(A^\dagger(y))                        #
    #                                                                                                #
    #       some other control parameters are implemented for testing performance                    #
    #       Now the default & confirmed setting is                                                   #
    #           [Md_tr, Md_alp, Md_sig, Ch_svd] = [0, 0, 0, -1]                                      #
    #       where Ch_svd: the choice of SVD for getting the initial X0                               #
    #                   =  0  the full SVD, i.e. LA.svd                                              #
    #                      1  scipy.sparse.linalg.svds                                               #
    #                      2  calling power_Largest_EigV = power method to have largest sig vec only #
    #                     -1  using the randomized SVD mentioned in the supplementary                #
    #         Note the SVD for each iteration still adopts the normal SVD                            # 
    #                                                                                                #             
    #      [Usage]  after 
    #                    worker = methodsRGD.BasicWorkerRGD(params_dict, input_S)                    #
    #                    worker.computeRGD(InitX_RGD, Ch_svd, Md_tr, Md_alp, Md_sig)                 #
    #               demonstrated by executing  RGD_optRun.py                                         #
    # ---------------------------------------------------------------------------------------------- #

    Md_tr =  0                  #   Method if including trace = 1
    Md_alp = 0                  #   method for scaling alpha
    Md_sig = 0                  #   method for scaling singular value
    Ch_svd = -1                 #   choice for initial SVD  (0: LA.svd, 1: svds;  2: power_Largest_EigV, -1: rSVD)

    InitX_RGD = 1;       #Ld_RGD = None                #   method of choosing initial X0

    Rpm = [InitX_RGD, Md_tr, Md_alp, Md_sig, Ch_svd]


    # ------------------------------------------------------------- #
    #   the above (A) (B) (C) defines ALL the parameters needed     #
    #                                                               #
    #   the below (D) start running the program                     # 
    # ------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    #   define some variables for                                      #
    #       recording running some information about the simulation    #
    #                                                                  #
    #   T_rec: record the projector & measurement & running time       #
    #                                                                  #
    #   Ttot : record  running time for calling                        # 
    #             MiFGD_Fun_Use_v5.py             (MiFGD)              #
    #             RGD_EnRui_SL_read_y_v5.py       (RGD)                #
    #                                                                  #
    #                                                                  #
    #   Wk_dict : record the parameters & optization result            #
    # ---------------------------------------------------------------- #
    T_rec = {}          
    T_rec['Nk'] = Nk
    T_rec['m']  = m

    Ttot = []
    Wk_dict = {}                #  record the worker list

    # ------------------------------------- #
    #    define controlling parameters      #     
    # ------------------------------------- #
    params_control = {'New_Pj_shot':  New_Pj_shot, 
                      'Pj_method':  Pj_method,
                      'mea_method': mea_method
                      }

    # ---------------------------------------------------------- #
    #   (D-0) generate the Naming of the directory to store      #
    #       according to the state & version defined in (A) (B)  # 
    # ---------------------------------------------------------- #

    Dir, Name, params_setup, version, proj_path, meas_path, Dir0, StVer = \
            State_Naming(Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, \
                        Pj_method, mea_method, measure_method, Data_In)

    # ------------------------------------------------- #
    #   (D-1) generate the projectors = Pauli matrices  #
    # ------------------------------------------------- #
    label_list, T_rec, Ncpu_Pj, saveP_bulk  = Gen_Projectors(params_setup, params_control, T_rec)   #  creating projectors

    # ------------------------------------------------------------------------------------------ #
    #   (D-2)  creating the input density matrix                                                 #
    #        & do the measurement                                                                #
    #                                                                                            #
    #     the input density matrix                                                               #
    #        = target_density_matrix (in the matrix form)                                        #
    #        = rho      (None for pure state | qu.rand_dm_ginibre from qutip for random states)  #
    #        = input_S  (pure state from qiskit circuit |  None for random mixed states)         #
    # ------------------------------------------------------------------------------------------ #

    target_density_matrix, input_S, T_rec, Ncpu_meas = Get_measurement_by_labels(params_setup, label_list, New_Pj_shot, \
                                StVer, T_rec, measure_method)

    # -------------------------------------------------------------- #
    #  get params_dict as the input for each optimization method     #
    #           i.e.  used in both MiFGD &  RGD                      #
    # -------------------------------------------------------------- #

    params_dict = Get_param.Initial_params(Dir, params_setup, target_density_matrix)

    if input_S:
        print(' ********************************************************* \n')
        print('   there exists pure state   ')
        params_dict['target_state'] = input_S.get_state_vector()            



    # ------------------------------------------------------ #
    #   (D-3)   start executing the optimization method      #
    #               MiFGD  &/or  RGD   method                #
    # ------------------------------------------------------ #

    if Call_RGD == 1:

        print('\n')
        print(' ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ')
        print(' ++           start calculating  RGD                 ++ ')
        print(' ++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n')

        exec(open('RGD_optRun.py').read())
        Ttot.append(('RGD', tw2b - tw2a))
        T_rec['RGD']  = tw2b - tw2a

        print("The total time for ALL is {}".format(Ttot))

        #wk_list.append(worker)
        #Wk_dict['RGD'] = [Rpm, worker]              #       Rpm = [InitX_RGD, Md_tr, Md_alp, Md_sig]
        Wk_dict['RGD'] = [Rpm, wc]              #       Rpm = [InitX_RGD, Md_tr, Md_alp, Md_sig]
        del worker

    if Call_MiFGD == 1:

        print(' ---------   to call MiFGD_call.py   -----------')

        for ii in range(0, Num_mu):

            mu = muList[ii]
            print('\n ----------------  {}-th mu = {}   start ------------------ \n'.format(ii, mu))

            params_dict['eta'] = 0.01       #  only for MiFGD
            params_dict['mu'] = mu          #  only for MiFGD

            exec(open('MiFGD_optRun.py').read())
            Ttot.append((mu, tw2b - tw2a))
            T_rec[mu] = tw2b - tw2a

            #wk_list.append(worker)
            #Wk_dict[mu] = worker
            Wk_dict[mu] = wc
            del worker

            #print(' ----------------  {}-th mu = {}   done  ----------------\n\n'.format(ii, mu))

        print("The total time for each mu is {}".format(Ttot))
    

    # ----------------------------- #
    #   show some results           #
    # ----------------------------- #
    print(' --------------------- \n   Target_Err:\n')
    #Target_Err = {}
    #for wk, mt in zip(wk_list, Ttot):
    #    print(' ({:.3},  {}), '.format(mt[0], wk.Target_Err_Xk[-1]))
    for result in Ttot:
        method = result[0]
        wk     = Wk_dict[method]
        if method == 'RGD':         #  or  isinstance(wk, list)
            wk = wk[1]              #  Wk_dict['RRD'] = [Rpm, worker]
        print(' ({:.3},  {}), '.format(method, wk.Target_Err_Xk[-1]))
    print('\n')


    if To_save_X0 == 1:             #  to save the initial X0 or not
        w_X0 = {}
        #for wk, mt in zip(wk_list, Ttot):
        for result in Ttot:
            method = result[0]
            wk     = Wk_dict[method]    #  Wk_dict[mu]    = worker
            if method == 'RGD':         #  or  isinstance(wk, list)
                wk = wk[1]              #  Wk_dict['RRD'] = [Rpm, worker]

            if wk.InitX >= 0:
                if wk.method == 'MiFGD':
                    w_X0[method] = (wk.U0, wk.X0)
                elif wk.method == 'RGD':
                    w_X0[method] = (wk.u0, wk.v0, wk.s0, wk.X0)

        F_X0 = '{}/X0.pickle'.format(meas_path)
        with open(F_X0, 'wb') as f:
            pickle.dump(w_X0, f)

    #exec(open('Show_Results.py').read())

    # ------------------------  #
    #   record some results     #
    # ------------------------  #
    Fname = '{}_Rec_info.txt'.format(Dir)
    print('Fname = {}'.format(Fname))
    f = open(Fname, 'a')

    f.write('\n ------------------------------------------------- \n')
    f.write(' Fname = {}\n'.format(Fname))
    f.write('     Nk       = {}\n'.format(Nk))
    f.write('     m        = {}\n'.format(m))
    f.write('     mea      = {}\n'.format(mea))
    f.write('     version  = {}\n'.format(version))
    f.write('     meas_path    = {}\n'.format(meas_path))

    f.write('     Pj_method    = {}\n'.format(Pj_method))
    f.write('     mea_method   = {}\n'.format(mea_method))
    f.write('     measure_method   = {}\n'.format(measure_method))
    f.write(' Time (projection)  = {}  [num_cpus] = {}\n'.format(T_rec['proj'], Ncpu_Pj))
    f.write(' Time (measurement) = {}  [num_cpus] = {}\n\n'.format(T_rec['measurement'], Ncpu_meas))
    
    for ii in range(len(Ttot)):
        md, dt = Ttot[ii]

        if isinstance(md, str):        # RGD  method
            Rpm, wk = Wk_dict[md]
            f.write('   {}, Rpm = [InitX, Md_tr, Md_alp, Md_sig, Ch_svd] = {}  -->  time = {}, '.format(md, Rpm, dt))

        else:                          # MiFGD method
            wk = Wk_dict[md]
            f.write('   mu  =  {:.3}                                                          -->  time = {}, '.format(md, dt))

        Err_rel_Xk    = wk.Err_relative_Xk[-1]
        Err_Target_Xk = wk.Target_Err_Xk[-1]
        Nite          = wk.iteration
        converged     = wk.converged

        f.write(' {} iter converge? {}: relative Err = {}, Target Err = {}\n'.format(Nite, converged, Err_rel_Xk, Err_Target_Xk))

    f.close()
