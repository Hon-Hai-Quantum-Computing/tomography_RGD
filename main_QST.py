
#import time
#import pickle                   # only for loading qiskit

import Get_param

from importlib import reload

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

def FileRec_Tproj_Tmea(File_rec, T_rec, params_setup, Ncpu_Pj):

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
                                pm_RGD, pm_MiFGD, T_rec, Ncpu_Pj):

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

            for InitX_MiFGD in [1, 0]:
            #for InitX_MiFGD in [0]:
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

    FileRec_Tproj_Tmea(File_rec, T_rec, params_setup, Ncpu_Pj)

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
    Call_RGD = 1        #  = 1: call the RGD   optimization to calculate

    # -------------------------------------------------------------------------------- #
    #   (C-1) parameters for calling MiFGD optimization algorithm                      #
    #                                                                                  #
    #     InitX_MiFGD = 0   generate random matrix according to MiFGD paper & code     #
    #                                                                                  #
    #     [Usage]   worker = methodsMiFGD.BasicWorker(params_dict, input_S)            #
    #               worker.compute(InitX_MiFGD)                                        #
    #     demonstrated by executing  MiFGD_optRun.py                                   #
    # -------------------------------------------------------------------------------- #

    #InitX_MiFGD = 0;     
    InitX_MiFGD = 1       # 0: random start,  1: MiFGD specified init

    #muList = [3/4, 0.5, 1/3, 0.2, 0.25, 0.125]
    #Num_mu = 1     #  number of mu for running MiFGD
                    #  (eg) = 2 --> run cases of muList = [3/4, 0.5]

    #muList = [3/4, 0.25, 4.5e-5]
    #muList = [4.5e-5, 3/4, 0.25]
    #Num_mu = 3      #  number of mu for running MiFGD

    #muList = [0.75]
    muList = [4.5e-5]
    Num_mu = 1      #  number of mu for running MiFGD

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
                                input_S, T_rec, Ncpu_Pj, pm_MiFGD, pm_RGD, List_alpha):

    T_Rho  = {}      #  record Ttot for each Rho
    Wk_Rho = {}      #  record Wk_dict for each Rho

    DirRho    = params_setup['DirRho']

    # --------------------------------------------- #
    #   start doing the tomography optimization     #
    #       using   RGD  &  MiFGD                   #
    # --------------------------------------------- #

    Ttot, Wk_dict = Tomography_over_measurement(params_setup, target_density_matrix, \
                             input_S, pm_RGD, pm_MiFGD, T_rec, Ncpu_Pj)

    T_Rho[DirRho]  = Ttot
    Wk_Rho[DirRho] = Wk_dict

    T_Rho2, Wk_Rho2 = Tune_Kappa_Tomography(params_setup, target_density_matrix, \
            input_S, T_rec, Ncpu_Pj, List_alpha, pm_MiFGD, pm_RGD)
    T_Rho.update(T_Rho2)
    Wk_Rho.update(Wk_Rho2)

    return T_Rho, Wk_Rho

def Tune_Kappa_Tomography(params_setup, target_density_matrix, input_S, T_rec, Ncpu_Pj, \
                List_alpha, pm_MiFGD, pm_RGD):
    """ (eg)        List_Kappa = [10, 20, 30]
                    List_alpha = [0.25, 0.5, 1, 2, 3]

    """

    T_Rho  = {}      #  record Ttot for each Rho
    Wk_Rho = {}      #  record Wk_dict for each Rho                    

    Nr        = params_setup['Nr']
    StateName = params_setup['StateName']
    StVer     = params_setup['StVer']

    if StateName != 'KapRnd':
        print('  NOT "KapRnd" --> not tuning Kappa')
        return T_Rho, Wk_Rho
    
    elif StateName == 'KapRnd':

        Data_In = params_setup['Data_In']

        if len(Data_In) >0 and Data_In[0] != '':        # DirRho specified
            raise TypeError(' The KapRnd not suitable for specifying DirRho in Data_In')
        elif len(Data_In) > 2:                          # DirMease specified
            raise TypeError(' The KapRnd not suitable for specifying DirMeas')
            

        if StVer[1] == 1:            #  need to create New Rho
            u, s, vh = Decompose_Rho(target_density_matrix, Nr, params_setup)
            del target_density_matrix

        for alpha in List_alpha:

            if StVer[1] == 1:       #  create New Rho
                NewRho, s_new, Kappa = Tune_Kappa(u, s, vh, Nr, alpha)
            elif StVer[1] == 0:     #  Rho already existed

                #raise ValueError(' need to Generate Rho for New Kappa, i.e. need StVer[1] = 1')
                
                NewRho = []
                s_new, Kappa = Gen_sNew_Power(alpha, Nr)

            params_setup, Dt_meas = Params_Kappa_change(params_setup, Kappa, alpha, NewRho, s_new)
            T_rec['measurement'] = Dt_meas

            if StVer[1] == 0:       #  to read existing Rho
                DirRho = params_setup['DirRho']

                with open('{}/RhoKappa.dat'.format(DirRho), 'rb') as f:
                    Kappa2, NewRho = pickle.load(f)

            Ttot, Wk_dict = Tomography_over_measurement(params_setup, NewRho, \
                                        input_S, pm_RGD, pm_MiFGD, T_rec, Ncpu_Pj)

            DirRho               = params_setup['DirRho']

            Wk_Rho[DirRho]       = Wk_dict
            T_Rho[DirRho]        = Ttot

    return T_Rho, Wk_Rho

def Add_More_Kappa(params_setup, T_rec, Ncpu_Pj, List_alpha, pm_MiFGD, pm_RGD):
    """ (eg)        List_Kappa = [40, 50]
                    List_alpha = [0.3, 1, 2]
    """

    StateName = params_setup['StateName']
    Nk        = params_setup['n']
    Nr        = params_setup['Nr']

    #DirRho    = params_setup['DirRho']
    Dir2m     = params_setup['Dir2m']
    DirRho    = '{}/Kap0'.format(Dir2m)

    StVer     = params_setup['StVer']

    #StVer[1] = 0
    #params_setup['StVer'] = StVer
    #StVerR0 = [StVer[0], 0]

    if StVer[1] == 0:
        raise ValueError(' Need to generate New Kappa, i.e. must StVer[1] = 1')

    T_rec['Ncpu_meas']    =  1   # same as  Get_measurement_by_labels

    target_density_matrix, rho = Gen_randM(Nk, StateName, Nr, StVer, DirRho)
    input_S = None

    T_Rho, Wk_Rho = Tune_Kappa_Tomography(params_setup, target_density_matrix, \
                    input_S, T_rec, Ncpu_Pj, List_alpha, pm_MiFGD, pm_RGD)

    return T_Rho, Wk_Rho



def Sample_Rnd_Tomography(params_setup, Samples, Num_Rho, T_rec, Ncpu_Pj, pm_MiFGD, pm_RGD, List_alpha):
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
                    target_density_matrix, input_S, T_rec, Ncpu_Pj, pm_MiFGD, pm_RGD, List_alpha)

    print(' ---------- Run samples of RND completed --------- \n\n')
    return params_setup


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

    #exec(open('parm_yaml.py').read())

    Choose = 1
    Nk, m, mea = basic_sys_setting(Choose)

    # ---------------------------------------------------------- #
    #   StateName choce = 'GHZ' | 'Had' |  'rand' |  'KapRnd'    #
    #                                                            #
    #   Default_MeasureState_parm() --> default params needed    #
    #             for  projectors  &  measurement                # 
    # ---------------------------------------------------------- #

    StateName = 'KapRnd'; Nr = 3;  Noise = 0;    StVer = [1, 1]
    #StateName = 'KapRnd';  Nr = 3;  Noise = 0;   StVer = [1, 0]    # load existing Rho

    #StateName = 'rand';   Nr = 3;  Noise = 0;    StVer = [1, 1]   # create new Rho
    #StateName = 'rand';   Nr = 3;  Noise = 0;    StVer = [1, 0]   # load existing Rho

    #StateName = 'GHZ';    Nr = 1;  Noise = -1;   StVer = 0      
    #StateName = 'quWst';    Nr = 1;  Noise = 0;   StVer = 0      

    # --------------------------------- #
    #   some default parameters         #
    # --------------------------------- #

    Pj_method = 1               #   the method to save | load  projectors
    mea_method = 1              #   the method to save | load  measurement_dict (count_dict) 

    measure_method = 3          #  = 1: direct label_list,  = 3: parallel cpu

    # ----------------------------------------------------- #
    #   version parameter for  Projectors & measurement     #
    #                          DirStore to store data       #
    # ----------------------------------------------------- #

    #DirStore = '../Tomography/calc'
    DirStore = '../Tomography/DataTest'


    version_choice = 1
    if version_choice == 1:     #  creating the whole new data set
        Data_In   = []

        New_Pj_shot = [-1, 1]  #  [New Proj?, New shot?]  | [0,0] loading | 

    elif version_choice == 2:   # Specify Dir of [Rho, Proj, Measure] 

        DirRho   = ''
        #Dir_proj = '{}/rand-3-r3/R1/rand-3-r3_m10_s1/rand-3-r3_m10_s1_Proj'.format(DirStore)
        #Dir_meas = '{}/rand-3-r3/R1/rand-3-r3_m10_s1/rand-3-r3_m10_s1_zN0_v1_measurements'.format(DirStore)

        Dir_proj = '{}/KapRnd-3-r3/Proj_m10_s1'.format(DirStore)
        Dir_meas = '{}/KapRnd-3-r3/R1/Kap0/Pm10_s1_zExact_v1'.format(DirStore)

        #Data_In   = ['', Dir_proj, Dir_meas]
        Data_In   = ['', Dir_proj]

        New_Pj_shot = [0, 0]  #  [New Proj?, New shot?]  | [0,0] loading | 

    elif version_choice == 3:   # Specify Dir of [Rho, Proj, Measure] 

        DirRho   = '{}/rand-3-r3/R1'.format(DirStore)
        Dir_proj = '{}/rand-3-r3/R1/rand-3-r3_m10_s1/rand-3-r3_m10_s1_Proj'.format(DirStore)
        Dir_meas = '{}/rand-3-r3/R1/rand-3-r3_m10_s1/rand-3-r3_m10_s1_zN0_v1_measurements'.format(DirStore)

        #Data_In   = [DirRho, Dir_proj, Dir_meas]
        Data_In   = [DirRho, Dir_proj]

        New_Pj_shot = [0, 1]  #  [New Proj?, New shot?]  | [0,0] loading | 


    version = [1, 2, Noise]   # [Proj version, measurement version, Noise]


    # ---------------------------------------------------------- #
    #   (D) start running the optimization program               #
    #                   for tomography                           # 
    #   (D-0) generate the Naming of the directory to store      #
    #       according to the state & version defined in (A) (B)  # 
    # ---------------------------------------------------------- #    

    params_setup = Params_define(Nk, StateName, m, mea, Nr, \
                Pj_method, mea_method, measure_method, DirStore) 

    params_setup = State_Naming(StVer, version, \
                        params_setup, New_Pj_shot, Data_In)
    del Nk, m, mea, Nr, Pj_method, mea_method, measure_method
    #del StVer, version, New_Pj_shot, Data_In

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

    pm_MiFGD, pm_RGD = Default_OptMethod_params()
    List_alpha = [5, 7, 10]     #  only for 'KapRnd' --> to tune Kappa

    Run_meas_Tomography = 1
    if Run_meas_Tomography == 1:      # do measurement & run tomography

        target_density_matrix, input_S, T_rec, Ncpu_meas = \
            Get_measurement_by_labels(params_setup, label_list, T_rec) 

        T_Rho, Wk_Rho = Default_Setting_Run_Tomography(params_setup, \
                        target_density_matrix, input_S, T_rec, Ncpu_Pj, \
                        pm_MiFGD, pm_RGD, List_alpha)



    # ----------------------------------------- #
    #   special cases:  run more cases          #
    #   [default]    Special_Usage = 0          #
    # ----------------------------------------- #

    Special_Usage = 2

    if Special_Usage == 1:    #  adding more Kappa cases  (only for KapRnd)

        if StateName == 'KapRnd':   # continue generate more Kappa

            List_alpha = [11, 12]

            T_Rho, Wk_Rho = Add_More_Kappa(params_setup, T_rec, Ncpu_Pj, List_alpha, pm_MiFGD, pm_RGD)


    elif Special_Usage == 2:    # having more samples (only for 'rand')

        Samples = 'sample1'
        Num_Rho  = 5

        params_setup = Sample_Rnd_Tomography(params_setup, Samples, \
                        Num_Rho, T_rec, Ncpu_Pj, pm_MiFGD, pm_RGD, List_alpha)


    print(' ******   Happy Ending   *****')