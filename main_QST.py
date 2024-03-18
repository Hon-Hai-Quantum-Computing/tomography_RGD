
#import time
#import pickle                   # only for loading qiskit

import Get_param

from importlib import reload
#from Utility import State_Naming, Gen_Rho, Gen_randM, state_measure, \
#                    state_measure_save, state_measure_wID, state_S_coef, split_list, mp_state_measure

from Utility import State_Naming
from Utility import Params_define
#from Utility import Ver_Rho_Naming, Ver_Proj_Naming, Ver_meas_Naming                 

from Input_Utility import Gen_Projectors 

from Input_measurements import Get_measurement_by_labels

from Input_setting import basic_sys_setting

from Utility import Decompose_Rho, Tune_Kappa
from Utility import Params_Kappa_change
from Utility import Params_RND_samples
from Utility import Gen_randM
from Utility import Gen_sNew_Power

from RGD_optRun import Run_RGD
from MiFGD_optRun import Run_MiFGD

import pickle

# ----------------------------------- #
#import sys
#sys.path.append('./quQST')

#import projectors
#import measurements
#from methodsRGD import Amea
#reload(projectors)


#import numpy as np
#import os
#import multiprocessing
#from itertools import repeat
#from functools import reduce

#import Input_measurements
#reload(Input_measurements)

#import Input_setting
#reload(Input_setting)



def Print_Target_Err(Ttot, Wk_dict):
    # ----------------------------- #
    #   show some results           #
    # ----------------------------- #
    print(' --------------------- \n   Target_Err:\n')
    #Target_Err = {}
    #for wk, mt in zip(wk_list, Ttot):
    #    print(' ({:.3},  {}), '.format(mt[0], wk.Target_Err_Xk[-1]))
    for result in Ttot:
        #method = result[0]
        method, InitX, Option, dt = result


        if method == 'RGD':         #  or  isinstance(wk, list)
            OptMd = '{}{}'.format(method, InitX)    #  Optimization method = RGD0  or RGD1
            #print('  Rmd  =  {}'.format(OptMd))

            wk  = Wk_dict[OptMd]
            wk = wk[1]              #  Wk_dict['RRD'] = [Rpm, worker]

        else:                       #  MiFGD
            if InitX == 0:
                OptMd = method
            elif InitX == 1:
                OptMd = -method
            else:
                OptMd = '{}+{}'.format(method, InitX)

            OptMd = '{}+{}'.format(OptMd, Option)

            #OptMd  = method
            #wk     = Wk_dict[method]
            wk     = Wk_dict[OptMd]

        #print(' (method, InitX, Err) = ({:.7}, {}, {}), dt = {}'.format(OptMd, \
        #                        InitX, wk.Target_Err_Xk[-1], dt))
        print(' (method, InitX, Option, Err) = ({:.7}, {}, {}, {}), dt = {}'.format(OptMd, \
                                InitX, Option, wk.Target_Err_Xk[-1], dt))
        if Option == 1:
            print(' ---------------------- ')  
        elif Option == 2:
            print(' ++++++++++++++++++++++ ')  
        elif OptMd == 'RGD0':
            print(' ********************** ')


    print('\n')

def FileRec_Tproj_Tmea(File_rec, T_rec, params_setup):

    Nk        = params_setup['n']
    m         = params_setup['num_labels']
    mea       = params_setup['num_shots']
    StVer     = params_setup['StVer']
    version   = params_setup['version']

    Pj_method  = params_setup['Pj_method']
    mea_method = params_setup['mea_method']
    measure_method = params_setup['measure_method']

    #StateName = params_setup['StateName']
    Data_In  = params_setup['Data_In']

    DirRho   = params_setup['DirRho']
    Dir_proj = params_setup['projector_store_path']
    Dir_meas = params_setup['measurement_store_path']

    print('File_rec = {}'.format(File_rec))
    f = open(File_rec, 'a')

    f.write('\n ------------------------------------------------- \n')
    f.write(' File_rec = {}\n'.format(File_rec))
    f.write('     Nk       = {}\n'.format(Nk))
    f.write('     m        = {}\n'.format(m))
    f.write('     mea      = {}\n'.format(mea))
    if len(Data_In) == 0:
        f.write('     StVer    = {}\n'.format(StVer))
        f.write('     version  = {}\n'.format(version))
        #f.write('     meas_path    = {}\n'.format(meas_path))
        f.write('     meas_path    = {}\n'.format(Dir_meas))
    else:
        f.write('     ****  read from Data_In -->  \n  ')
        f.write('     ****     DirRho   = {}\n  '.format(DirRho))
        f.write('     ****     Dir_proj = {}\n  '.format(Dir_proj))
        f.write('     ****     Dir_meas = {}\n\n'.format(Dir_meas))

    f.write('     Pj_method    = {}\n'.format(Pj_method))
    f.write('     mea_method   = {}\n'.format(mea_method))
    f.write('     measure_method   = {}\n'.format(measure_method))
    f.write(' Time (projection)  = {}  [num_cpus] = {}\n'.format(T_rec['proj'], Ncpu_Pj))

    Ncpu_meas = T_rec['Ncpu_meas']
    f.write(' Time (measurement) = {}  [num_cpus] = {}\n\n'.format(T_rec['measurement'], Ncpu_meas))

    f.close()


def FileRec_Tmethod(File_rec, Ttot, Wk_dict):
    # ------------------------------------  #
    #   record some results in File_rec     #
    # ------------------------------------  #
 
    print('  File_rec = {}'.format(File_rec))

    f = open(File_rec, 'a')        #  File_rec = '{}_Rec_info.txt'.format(Dir)   defined earlier
    for ii in range(len(Ttot)):
        md, InitX, Option, dt = Ttot[ii]
        print(' md = {},  InitX = {}, Option = {}, dt = {}'.format(md, InitX, Option, dt))

        if isinstance(md, str):                  #   RGD  method
            method = '{}{}'.format(md, InitX)    # = RGD0  or RGD1

            Rpm, wk = Wk_dict[method]
            f.write('   {}, Rpm = [InitX, Md_tr, Md_alp, Md_sig, Ch_svd] = {}  -->  time = {}, '.format(method, Rpm, dt))

        else:                           # MiFGD method
            if InitX == 0:
                method = md
            elif InitX == 1:
                method = -md   
            else:
                method = '{}+{}'.format(md, InitX) 
               
            method = '{}+{}'.format(method, Option)    

            #wk = Wk_dict[md]
            wk = Wk_dict[method]

            eta   = wk.eta
            #etaTh = wk.etaTh            # only for Option = 1
            #f.write('   mu  =  {:.3}                                                          -->  time = {}, '.format(md, dt))
            #f.write('   mu  =  {:.3},  InitX = {} |=> method = {}                            -->  time = {}, '.format(md, InitX, method, dt))
            #f.write('   mu  =  {:.3},  InitX = {}, eta = {} |=> method = {}                   -->  time = {}, '.format(md, InitX, eta, method, dt))
            #f.write('   mu  =  {:.3},  InitX = {}, eta = {}, etaTh = {} |=> method = {}       -->  time = {}, '.format(md, InitX, eta, etaTh, method, dt))
            f.write('   mu  =  {:.3},  InitX = {}, eta = {:.3}, Op = {} |=> method = {}          -->  time = {}, '.format(md, \
                                                                InitX, eta, Option, method, dt))

        Err_rel_Xk    = wk.Err_relative_Xk[-1]
        Err_Target_Xk = wk.Target_Err_Xk[-1]
        Nite          = wk.iteration
        converged     = wk.converged

        f.write(' {} iter converge? {}: relative Err = {}, Target Err = {}\n'.format(Nite, converged, Err_rel_Xk, Err_Target_Xk))

    f.close()

def Tomography_over_measurement(params_setup, target_DM, input_S, \
                                pm_RGD, pm_MiFGD, T_rec):

    # -------------------------------------------------------------- #
    #  get params_dict as the input for each optimization method     #
    #           i.e.  used in both MiFGD &  RGD                      #
    # -------------------------------------------------------------- #

    #params_dict = Get_param.Initial_params(Dir, params_setup, target_density_matrix)
    params_dict = Get_param.Initial_params(params_setup, target_DM, input_S)

    # --------------------------------------------------------------- #
    #   extract parameters for (each optimization method)/recording   #
    # --------------------------------------------------------------- #

    Call_MiFGD, InitX_MiFGD, muList, Num_mu = pm_MiFGD
    Call_RGD, Rpm   = pm_RGD

    InitX_RGD, Md_tr, Md_alp, Md_sig, Ch_svd = Rpm

    # ---------------------------------------------------------------- #
    #   define some variables for                                      #
    #       recording running some information about the simulation    #
    #                                                                  #
    #   T_rec: record the projector & measurement & running time       #
    #                                                                  #
    #   Ttot : record  running time for calling                        # 
    #             MiFGD_optRun.py             (MiFGD)                  #
    #             RGD_optRun.py       (RGD)                            #
    #                                                                  #
    #   Wk_dict : record the parameters & optization result            #
    # ---------------------------------------------------------------- #

    Ttot = []
    Wk_dict = {}                #  record the worker list

    # ----------------------------------------------------------- #
    #   start executing the optimization method                   #
    #               MiFGD  &/or  RGD   method                     #
    #           &   record method running time = tw2b - tw2a      #
    # ----------------------------------------------------------- #


    if Call_RGD == 1:

        for InitX_RGD in [1, 0]:

            Rpm = InitX_RGD, Md_tr, Md_alp, Md_sig, Ch_svd

            #exec(open('RGD_optRun.py').read())
            Frec_RGD, wc, RunTime = Run_RGD(params_dict, Rpm)

            #Ttot.append(('RGD', RunTime))
            #T_rec['RGD']  = RunTime
            #Wk_dict['RGD'] = [Rpm, wc]              #  Rpm = [InitX_RGD, Md_tr, Md_alp, Md_sig]

            method = 'RGD{}'.format(InitX_RGD)

            #Ttot.append(('RGD', InitX_RGD, RunTime))
            Ttot.append(('RGD', InitX_RGD, 0, RunTime))

            T_rec[method]  = RunTime
            Wk_dict[method] = [Rpm, wc]              #  Rpm = [InitX_RGD, Md_tr, Md_alp, Md_sig]

            print('   RGD method = {}'.format(method))
            print("The total time for ALL is {}".format(Ttot))


    if Call_MiFGD == 1:

        for ii in range(0, Num_mu):

            mu  = muList[ii]
            eta = 0.01

            print('\n ----------------  {}-th mu = {}   start ------------------ \n'.format(ii, mu))

            #for InitX_MiFGD in [1, 0]:
            for InitX_MiFGD in [0]:
            #for InitX_MiFGD in [1]:
                #for Option in [0, 1, 2]: 
                #for Option in [2, 1, 0]: 
                #for Option in [1, 2]: 
                for Option in [2]: 
                #for Option in [0]: 

                    Rpm_MiFGD = [InitX_MiFGD, mu, eta, Option]

                    Frec_MiFGD, wc, RunTime = Run_MiFGD(params_dict, Rpm_MiFGD)

                    #Ttot.append((mu, InitX_MiFGD, RunTime))
                    Ttot.append((mu, InitX_MiFGD, Option, RunTime))                    

                    if InitX_MiFGD == 0:
                        method = mu
                    elif InitX_MiFGD == 1:
                        method = -mu
                    else:
                        method = '{}+{}'.format(mu, InitX_MiFGD)

                    method = '{}+{}'.format(method, Option)

                    #T_rec[mu] = RunTime
                    #Wk_dict[mu] = wc

                    T_rec[method] = RunTime
                    Wk_dict[method] = wc


            #print(' ----------------  {}-th mu = {}   done  ----------------\n\n'.format(ii, mu))

            print("The total time for each mu is {}".format(Ttot))


    # ----------------------------------------- #
    #   (D-2)'     file record information      #
    # ----------------------------------------- #
    StateName = params_setup['StateName']
    DirRho    = params_setup['DirRho']
    Dir2m     = params_setup['Dir2m']

    if StateName == 'KapRnd':
        File_rec = '{}/Rec_info.txt'.format(DirRho)
    else:
        File_rec = '{}_Rec_info.txt'.format(Dir2m)

    # ----------------------------- #
    #   output & record results     #
    # ----------------------------- #  

    FileRec_Tproj_Tmea(File_rec, T_rec, params_setup)

    FileRec_Tmethod(File_rec, Ttot, Wk_dict)


    Print_Target_Err(Ttot, Wk_dict)

    print('    DirRho    = {}'.format(params_setup['DirRho']))
    print('  measurement_store_path = ')
    print('            {}'.format(params_setup['measurement_store_path']))

    return Ttot, Wk_dict

def Default_OptMethod_params():

    # ----------------------------------------------------------------- #
    #   (C) to call each optimization method to calculate               #
    #       if chosen to run,                                           #
    #          each optimization method will be executed sequentially   #
    #               with recording the result & run time                #
    #       MiFGD & RGD can be called separately                        #
    # ----------------------------------------------------------------- #

    Call_MiFGD = 1      #  = 1: call the MiFGD optimization to calculate
    Call_RGD = 0        #  = 1: call the RGD   optimization to calculate

    # -------------------------------------------------------------------------------- #
    #   (C-1) parameters for calling MiFGD optimization algorithm                      #
    #                                                                                  #
    #     InitX_MiFGD = 0   generate random matrix according to MiFGD paper & code     #
    #                                                                                  #
    #     [Usage]   worker = methodsMiFGD.BasicWorker(params_dict, input_S)            #
    #               worker.compute(InitX_MiFGD)                                        #
    #     demonstrated by executing  MiFGD_optRun.py                                   #
    # -------------------------------------------------------------------------------- #

    #InitX_MiFGD = 0;     #Ld_MiFGD = None
    InitX_MiFGD = 1

    #muList = [3/4, 0.5, 1/3, 0.2, 0.25, 0.125]
    #Num_mu = 1     #  number of mu for running MiFGD
                    #  (eg) = 2 --> run cases of muList = [3/4, 0.5]
    #muList = [3/4, 0.25, 4.5e-5]
    #Num_mu = 3      #  number of mu for running MiFGD

    #muList = [0.75]
    muList = [4.5e-5]
    Num_mu = 1      #  number of mu for running MiFGD

    #muList = [3/4, 4.5e-5]
    #muList = [4.5e-5, 3/4]
    #Num_mu = 2      #  number of mu for running MiFGD

    #muList = [0.25]
    #Num_mu = 1      #  number of mu for running MiFGD

    pm_MiFGD = [Call_MiFGD, InitX_MiFGD, muList, Num_mu]

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

    InitX_RGD = 1;              #Ld_RGD = None       #   method of choosing initial X0

    Rpm    = [InitX_RGD, Md_tr, Md_alp, Md_sig, Ch_svd]
    pm_RGD = [Call_RGD, Rpm]

    return pm_MiFGD, pm_RGD

def Default_Setting_Run_Tomography(params_setup, target_density_matrix,\
                                input_S, T_rec):

    pm_MiFGD, pm_RGD = Default_OptMethod_params()

    T_Rho  = {}      #  record Ttot for each Rho
    Wk_Rho = {}      #  record Wk_dict for each Rho

    StateName = params_setup['StateName']
    DirRho    = params_setup['DirRho']

    # --------------------------------------------- #
    #   start doing the tomography optimization     #
    #       using   RGD  &  MiFGD                   #
    # --------------------------------------------- #

    Ttot, Wk_dict = Tomography_over_measurement(params_setup, target_density_matrix, \
                             input_S, pm_RGD, pm_MiFGD, T_rec)

    T_Rho[DirRho]  = Ttot
    Wk_Rho[DirRho] = Wk_dict

    if StateName == 'KapRnd':
        #List_Kappa = [10, 20, 30, 40]
        #List_alpha = [0.25, 0.5, 1, 2, 3]
        #List_Kappa = [20, 40, 60, 80]
        #List_alpha = [0.5, 2, 4]
        #List_Kappa = [80]
        
        #List_alpha = [4, 8, 10, 15, 20]
        #List_alpha = [8, 10, 15, 20]
        List_alpha = [4, 6, 8, 9, 10, 15, 20]   # paper use

        T_Rho2, Wk_Rho2 = Tune_Kappa_Tomography(params_setup, target_density_matrix, \
                    input_S, T_rec, List_alpha)

        T_Rho.update(T_Rho2)
        Wk_Rho.update(Wk_Rho2)

    return T_Rho, Wk_Rho

def Tune_Kappa_Tomography(params_setup, target_density_matrix, input_S, T_rec, \
                List_alpha):
    """ (eg)        List_Kappa = [10, 20, 30]
                    List_alpha = [0.25, 0.5, 1, 2, 3]
    """
    pm_MiFGD, pm_RGD = Default_OptMethod_params()

    T_Rho  = {}      #  record Ttot for each Rho
    Wk_Rho = {}      #  record Wk_dict for each Rho

    Nr        = params_setup['Nr']
    StateName = params_setup['StateName']
    StVer     = params_setup['StVer']
    
    if StateName == 'KapRnd':
        if StVer[1] == 1:            #  need to create New Rho
            u, s, vh = Decompose_Rho(target_density_matrix, Nr, params_setup)
            del target_density_matrix

        for alpha in List_alpha:

            if StVer[1] == 1:       #  create New Rho
                NewRho, s_new, Kappa = Tune_Kappa(u, s, vh, Nr, alpha)
            elif StVer[1] == 0:     #  Rho already existed
                NewRho = []
                s_new, Kappa = Gen_sNew_Power(alpha, Nr)

            params_setup, Dt_meas = Params_Kappa_change(params_setup, Kappa, alpha, NewRho, s_new)
            T_rec['measurement'] = Dt_meas

            if StVer[1] == 0:       #  to read existing Rho
                DirRho = params_setup['DirRho']

                with open('{}/RhoKappa.dat'.format(DirRho), 'rb') as f:
                    Kappa2, NewRho = pickle.load(f)

            Ttot, Wk_dict = Tomography_over_measurement(params_setup, NewRho, \
                                        input_S, pm_RGD, pm_MiFGD, T_rec)

            DirRho               = params_setup['DirRho']

            Wk_Rho[DirRho]       = Wk_dict
            T_Rho[DirRho]        = Ttot

    return T_Rho, Wk_Rho

def Add_More_Kappa(params_setup, T_rec, List_alpha):
    """ (eg)        List_Kappa = [40, 50]
                    List_alpha = [0.3, 1, 2]
    """

    StateName = params_setup['StateName']
    Nk        = params_setup['n']
    Nr        = params_setup['Nr']

    DirRho    = params_setup['DirRho']
    StVer     = params_setup['StVer']
    
    StVer[1] = 0
    params_setup['StVer'] = StVer

    T_rec['Ncpu_meas']    =  1   # same as  Get_measurement_by_labels

    target_density_matrix, rho = Gen_randM(Nk, StateName, Nr, StVer, DirRho)
    input_S = None

    T_Rho, Wk_Rho = Tune_Kappa_Tomography(params_setup, target_density_matrix, \
                    input_S, T_rec, List_alpha)

    return T_Rho, Wk_Rho



def Sample_Rnd_Tomography(params_setup, Samples, Num_Rho, T_rec):
    """  (eg)      Samples = 'sample4'
                  Num_Rho  = 3

    """

    New_Pj_shot = params_setup['New_Pj_shot']
    StateName   = params_setup['StateName']

    #if New_Pj_shot[1] != 1:
    #    print(' No new measurement needed')
    #    return params_setup
    
    if StateName != 'rand':
        print('StateName = {}'.format(StateName))
        print('  -->   NOT considered here \n\n')
        return params_setup
 
    for numID in range(1, Num_Rho+1):
    #for numID in [4, 5]:

        params_setup = Params_RND_samples(params_setup, Samples, numID)

        target_density_matrix, input_S, T_rec, Ncpu_meas = \
            Get_measurement_by_labels(params_setup, label_list, T_rec)

        T_Rho, Wk_Rho = Default_Setting_Run_Tomography(params_setup, \
                    target_density_matrix, input_S, T_rec)

    print(' ---------- Run samples of RND completed --------- \n\n')
    return params_setup


def Default_MeasureState_parm(Obtain_Pj):
    """  (eg)       Obtain_Pj = 2

    """
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
    #       Data_In =  []                 defult                                     #
    #                  [DirRho, Dir_proj] if specified (only for testing)            #
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
    #       -->  Caveat:   Now only this works for call cases                        #
    #                   other choice of New_Pj_shot may not work                     #    
    #                       (eg) New_Pj_shot = [1, 2] for StateName = 'KapRnd'       #
    # ------------------------------------------------------------------------------ #
    
    if Obtain_Pj == 1:      #  direct read in data for tomography

        #StateName = 'GHZ';  Nk = 3; Nr = 1; m = 10; Noise = -1
        #DirRho   = '../Tomography/calc/GHZ-3/GHZ-3_m10_s7'
        #Dir_proj = '../Tomography/calc/GHZ-3/GHZ-3_m10_s7/GHZ-3_m10_s7_Proj'

        #StateName = 'rand';  Nk = 3; Nr = 3; m = 10; Noise = 0
        #DirRho    = '../Tomography/calc/rand-3-r3/R14/rand-3-r3_m10_s1'
        #Dir_proj  = '../Tomography/calc/rand-3-r3/R14/rand-3-r3_m10_s1/rand-3-r3_m10_s1_Proj'
        #Dir_proj  = '../Tomography/calc/rand-3-r3/rand-3-r3_m10_s1_Proj'
        #Dir_proj  = '../Tomography/calc/rand-3-r3/Pj_s2'

        StateName = 'KapRnd';  Nk = 3; Nr = 3; m = 10; Noise = 0
        #DirRho    = '../Tomography/calc/KapRnd-3-r3/R28/Kap0'
        #Dir_proj  = '../Tomography/calc/KapRnd-3-r3/Proj_m10_s27/'
        #DirRho   = '../Tomography/calc/KapRnd-3-r3/R53/Kap0'
        #Dir_proj = '../Tomography/calc/KapRnd-3-r3/Proj_m10_s53/'
        #DirRho   = '../Tomography/data/KapRnd-8-r3/R1/Kap0'
        #Dir_proj = '../Tomography/data/KapRnd-8-r3/Proj_m13107_s1'

        #DirRho   = '../Tomography/data/KapRnd-8-r3/R2/Kap0'
        #Dir_proj = '../Tomography/data/KapRnd-8-r3/Proj_m19660_s1'

        DirRho   = '../Tomography/data/KapRnd-3-r3/R6/Kap0'
        Dir_proj = '../Tomography/data/KapRnd-3-r3/Proj_m10_s6'

        Data_In     = [DirRho, Dir_proj]
        New_Pj_shot = [0, 1]  #  [New Proj?, New shot?]  | [0,0] loading | 

    elif Obtain_Pj == 2:        #  only load Dir_proj

        #Dir_proj  = '../Tomography/calc/rand-3-r3/Pj_s1'
        #Dir_proj  = '../Tomography/calc/rand-8-r3/rand-8-r3_m13107_s1_Proj'
        #Dir_proj  = '../Tomography/calc/rand-10-r3/rand-10-r3_m209715_s1_Proj'
        #Dir_proj  = '../Tomography/calc/rand-12-r3/rand-12-r3_m838860_s1_Proj'

        #Dir_proj  = '../Tomography/data/rand-6-r3/rand-6-r3_m1228_s1_Proj'
        Dir_proj  = '../Tomography/data/rand-8-r3/rand-8-r3_m13107_s1_Proj'
        #Dir_proj  = '../Tomography/data/rand-8-r3/rand-8-r3_m19660_s1_Proj'
        #Dir_proj  = '../Tomography/data/rand-10-r3/rand-10-r3_m209715_s1_Proj'
        #Dir_proj  = '../Tomography/data/rand-12-r3/rand-12-r3_m838860_s1_Proj'

        Data_In     = ['', Dir_proj]
        #New_Pj_shot = [0, 1]  #  [New Proj?, New shot?]  | [0,0] loading | 
        New_Pj_shot = [0, 0]  #  [New Proj?, New shot?]  | [0,0] loading | 

    elif Obtain_Pj == -1:     #  the label_list from single|two sites, instead of sampling
        Data_In = []

        New_Pj_shot = [1, 1]  #  still need to produce new label_list  

    elif Obtain_Pj == -2:
        Data_In = []
        New_Pj_shot = [0, 0]  #  [New Proj?, New shot?]  | [0,0] loading | 

    elif Obtain_Pj == 0:
        Data_In = []

        New_Pj_shot = [1, 1]  #  [New Proj?, New shot?]  | [0,0] loading | 
        #New_Pj_shot = [0, 1]  #  [New Proj?, New shot?]  | [0,0] loading | 
        #New_Pj_shot = [2, 1]  #  [New Proj?, New shot?]  | [0,0] loading | 
        #New_Pj_shot = [2, 0]  #  [New Proj?, New shot?]  | [0,0] loading | 
        #New_Pj_shot = [0, 2]  #  [New Proj?, New shot?]  | [0,0] loading | 

    # ------------------------------------- #
    #    define controlling parameters      #     
    # ------------------------------------- #
    params_control = {'New_Pj_shot':  New_Pj_shot, 
                      'Pj_method':  Pj_method,
                      'mea_method': mea_method
                      }
    return Pj_method, mea_method, measure_method, New_Pj_shot, Data_In


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
    # --------------------------------------------------------------------------------- #
    #   (A-1)  choosing the input state                                                 #
    #                                                                                   #
    #          Nr = number of rank                                                      #
    #                                                                                   #
    #       Noise =  -1          shot measurements (for pure states)                    #
    #                 0          get exact coef   (perfect measurement)                 #
    #                       now only exact value is used (by qutip)                     #
    #               [may implment some noise model in the future]                       #
    #                  i.e. [Gaussian model, mean, variance] --> not implemented yet    #
    #                                                                                   #
    #       StVer    = version control for random mixed states                          #
    #       StVer[0] = version of generated random density matrices of rank Nr          #
    #       StVer[1] =  0    loading random density matrices (generated earlier)        #
    #                           --> will check if the directory exists or not           #
    #                   1    generate new random density matrices of rank Nr            # 
    #                           if file already exists for specified version StVer[0]   #
    #                           then the version StVer[0] will +1 automatically         #
    #  i.e.  StVer  =  [version R?, New or Generate]   for (random) State               #
    # --------------------------------------------------------------------------------- #

    Choose = 7
    Nk, m, mea = basic_sys_setting(Choose)

    # ---------------------------------------------------------- #
    #   StateName choce = 'GHZ' | 'Had' |  'rand' |  'KapRnd'    #
    #                                                            #
    #   Default_MeasureState_parm() --> default params needed    #
    #             for  projectors  &  measurement                # 
    # ---------------------------------------------------------- #

    #StateName = 'KapRnd'; Nr = 3;  Noise = 0;    StVer = [1, 1]
    #StateName = 'KapRnd';  Nr = 3;  Noise = 0;   StVer = [1, 0]    # load existing Rho

    #StateName = 'rand';   Nr = 3;  Noise = 0;    StVer = [1, 1]   # create new Rho
    #StateName = 'rand';   Nr = 3;  Noise = 0;    StVer = [1, 0]   # load existing Rho

    StateName = 'Had';    Nr = 1;  Noise = -1;   StVer = 0      
    #StateName = 'quGH';    Nr = 1;  Noise = 0;   StVer = 0      
    #StateName = 'quWst';    Nr = 1;  Noise = 0;   StVer = 0      

    #Obtain_Pj = 0      #  0: whole new; 1: read Dir_In;  2: read Dir_Proj
    Obtain_Pj = -2      #  0: whole new; 1: read Dir_In;  2: read Dir_Proj
    #Obtain_Pj = 2      #  0: whole new; 1: read Dir_In;  2: read Dir_Proj

    Pj_method, mea_method, measure_method, New_Pj_shot, Data_In = \
                Default_MeasureState_parm(Obtain_Pj)

    version = [1, 1, Noise]   # [Proj version, measurement version, Noise]


    # ---------------------------------------------------------- #
    #   (D) start running the optimization program               #
    #                   for tomography                           # 
    #   (D-0) generate the Naming of the directory to store      #
    #       according to the state & version defined in (A) (B)  # 
    # ---------------------------------------------------------- #

    params_setup = Params_define(Nk, StateName, m, mea, Nr, \
                Pj_method, mea_method, measure_method, Obtain_Pj) 

    params_setup = State_Naming(StVer, version, \
                        params_setup, New_Pj_shot, Data_In)
    #del Nk, m, mea, Nr, Pj_method, mea_method, measure_method
    del StVer, version, New_Pj_shot, Data_In

    # ------------------------------------------------- #
    #   (D-1) generate the projectors = Pauli matrices  #
    # ------------------------------------------------- #

    label_list, T_rec, Ncpu_Pj, saveP_bulk  = Gen_Projectors(params_setup)   #  creating projectors

    # ------------------------------------------------------------------------------------------ #
    #   (D-2)  creating the input density matrix                                                 #
    #        & do the measurement                                                                #
    #                                                                                            #
    #     the input density matrix                                                               #
    #        = target_density_matrix (in the matrix form)                                        #
    #        = rho      (None for pure state | qu.rand_dm_ginibre from qutip for random states)  #
    #        = input_S  (pure state from qiskit circuit |  None for random mixed states)         #
    # ------------------------------------------------------------------------------------------ #
    # --------------------------------------------------------------------- #
    #   (D-3) Using the optimization method (MiFGD  or RGD)                 #
    #       to do the tomography with                                       #
    #             some default parameters needed by MiFGD/RGD respectively  #
    # --------------------------------------------------------------------- #

    if Obtain_Pj == 0 or Obtain_Pj == -1 or Obtain_Pj == -2: # usually run this

        #target_density_matrix, input_S, T_rec, Ncpu_meas, s_label, yProj = \
        target_density_matrix, input_S, T_rec, Ncpu_meas = \
            Get_measurement_by_labels(params_setup, label_list, T_rec)

        T_Rho, Wk_Rho = Default_Setting_Run_Tomography(params_setup, \
                        target_density_matrix, input_S, T_rec)

    elif Obtain_Pj == 1:    #  continue from existing Rho

        if StateName == 'KapRnd':   # continue generate more Kappa
            #List_Kappa = [100]
            #List_alpha = [0.25, 0.5, 1, 2, 3]
            #List_Kappa = [80]

            List_alpha = [8, 14]

            T_Rho, Wk_Rho = Add_More_Kappa(params_setup, T_rec, List_alpha)


    elif Obtain_Pj == 2:    # having more samples (only for 'rand')

        Samples = 'sample1'
        Num_Rho  = 5

        params_setup = \
            Sample_Rnd_Tomography(params_setup, Samples, Num_Rho, T_rec)


    print(' ******   Happy Ending   *****')