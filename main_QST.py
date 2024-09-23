

import Get_param

from importlib import reload

from Utility import State_Naming
from Utility import Params_define

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




def Print_Target_Err(Ttot, Wk_dict):
    """ To print the final target error message for each case

    Args:
        Ttot (dict): dictionary of recording running time for each case
        Wk_dict (dict): dctionary of recording the optimization result of each case 
    """
    # ----------------------------- #
    #   show some results           #
    # ----------------------------- #
    print(' --------------------- \n   Target_Err:\n')

    for result in Ttot:
        method, InitX, Option, dt = result

        if method == 'RGD':         #  or  isinstance(wk, list)
            OptMd = '{}{}'.format(method, InitX)    #  Optimization method = RGD0  or RGD1

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

            wk     = Wk_dict[OptMd]

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
    """ To print some information about the settings, including the measurement method

    Args:
        File_rec (str): File name of recording
        T_rec (dict): dictionary of recording the running time
        params_setup (dict): dictionary of parameters 
        Ncpu_Pj (int): number of processors for parallel generation of Pauli projectors 
    """

    Nk        = params_setup['n']
    m         = params_setup['num_labels']
    mea       = params_setup['num_shots']

    Pj_method  = params_setup['Pj_method']
    mea_method = params_setup['mea_method']
    measure_method = params_setup['measure_method']

    StateName = params_setup['StateName']
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

        if StateName == 'rand' or StateName == 'KapRnd':
            f.write('     (Generate New rand, rand matrix version) = ({}, {})\n'.format(
                params_setup['Generate New rand'], params_setup['rand matrix version']))

        f.write('     (Projector   version, measurement version) = ({}, {})\n'.format(
            params_setup['Proj version'], params_setup['measure version']))

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
    """ To record the information of different methods for comparison, including their
        parameter settings, running time, and target errors

    Args:
        File_rec (str): File name of recording
        Ttot (dict): dictionary of recording the running time 
        Wk_dict (dict): dictionary of recording the optimization results
    """
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

            wk = Wk_dict[method]

            eta   = wk.eta
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
    """ Do the tomography for given the measurement results

    Args:
        params_setup (dict): dictionary of parameters
        target_DM (nd.array): density matrix of the target state 
        input_S (class instance): constructed circuit for the target state
        pm_RGD (dict): dictionary of parameters for running the RGD method
        pm_MiFGD (dict): dictionary of parameters for running the MiFGD method 
        T_rec (dict): recording the running time 
        Ncpu_Pj (int): number of processors for parallel generation of Pauli operators

    Returns:
        list: (Ttot) list of recording running times
        dict: (Wk_dict) dictionary of recording running optimization results
        dict: (params_setup) updated dictionary of parameters for running each case
        dict: (T_rec) dictionary of recording running times
    """

    # -------------------------------------------------------------- #
    #  get params_dict as the input for each optimization method     #
    #           i.e.  used in both MiFGD &  RGD                      #
    # -------------------------------------------------------------- #

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

            method = 'RGD{}'.format(InitX_RGD)

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
                for Option in [2]: 

                    Rpm_MiFGD = [InitX_MiFGD, mu, eta, Option]

                    Frec_MiFGD, wc, RunTime = Run_MiFGD(params_dict, Rpm_MiFGD)

                    Ttot.append((mu, InitX_MiFGD, Option, RunTime))                    

                    if InitX_MiFGD == 0:
                        method = mu
                    elif InitX_MiFGD == 1:
                        method = -mu
                    else:
                        method = '{}+{}'.format(mu, InitX_MiFGD)

                    method = '{}+{}'.format(method, Option)

                    T_rec[method] = RunTime
                    Wk_dict[method] = wc


            #print(' ----------------  {}-th mu = {}   done  ----------------\n\n'.format(ii, mu))

            print("The total time for each mu is {}".format(Ttot))

    return Ttot, Wk_dict, params_setup, T_rec

def Rec_Info(Ttot, Wk_dict, params_setup, T_rec):
    """ To record information about each measurement case and all the methods applied to it

    Args:
        Ttot (dict): dictionary of recording run times
        Wk_dict (dict): dictionary of recording the optimization results 
        params_setup (dict): dictionary of parameters 
        T_rec (dict): dictionary of recording the run times 
    """
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


def Default_OptMethod_params():
    """ to obtain some parameters for running the optimizer (MiFGD and RGD)

    Returns:
        list: (pm_MiFGD) = [Call_MiFGD, InitX_MiFGD, muList, Num_mu]

            Call_MiFGD  = 1/0: to call the MiFGD optimization to calculate or not

            InitX_MiFGD = 1    # 0: random start,  1: MiFGD specified init
            muList      = list of mu parameters, (eg) [3/4, 0.25, 4.5e-5]
            Num_mu      = number of mu for running MiFGD = len(muList)


        list: (pm_RGD) = [Call_RGD, Rpm], 
                        where Rpm = [InitX_RGD, Md_tr, Md_alp, Md_sig, Ch_svd]

            Call_RGD    = 1/0: to call the RGD   optimization to calculate or not

            InitX_RGD = 1;    # 0: random start,  1: RGD specified init
            Md_tr  = Method if including unit trace, only use default 0
            Md_alp = method for scaling alpha, only use default 0
            Md_sig = method for scaling singular value, only use default 0
            Ch_svd = (default -1) choice for initial SVD (0: LA.svd, 1: svds;  2: power_Largest_EigV, -1: rSVD)

    """

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

    InitX_MiFGD = 1       # 0: random start,  1: MiFGD specified init

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
                                input_S, T_rec, Ncpu_Pj, pm_MiFGD, pm_RGD):
    """ Given the measurement results (read from directory specified by params_setup)
        and do the tomography according to the optimization method parameters
            where pm_MiFGD is the parameter for the optimization MiFGD method
              and pm_RGD is the parameter for the optimization RGD method,    
            
        and return the time taken and the output of the tomography.

    Args:
        params_setup (dict): dictionary of parameters
        target_density_matrix (nd.array): density matrix of the target state
        input_S (class instance): circuit implementation of the target state
        T_rec (dict): dictionary to update the running time
        Ncpu_Pj (int): number of cpu to produce the Pauli projectors parallelly 
        pm_MiFGD (dict): parameters specific to run the MiFGD optimization method 
        pm_RGD (dict): parameters specific to run the RGD optimization method 
        List_alpha (list): list of power bases to construct the singular values if provided

    Returns:
        dict: (T_Rho) to record the running time for different running cases
        dict: (Wk_Rho) to record the output for different running cases

    """

    T_Rho  = {}      #  record Ttot for each Rho
    Wk_Rho = {}      #  record Wk_dict for each Rho


    # --------------------------------------------- #
    #   start doing the tomography optimization     #
    #       using   RGD  &  MiFGD                   #
    # --------------------------------------------- #

    Ttot, Wk_dict, params_setup, T_rec = \
        Tomography_over_measurement(params_setup, target_density_matrix, \
                             input_S, pm_RGD, pm_MiFGD, T_rec)
    Rec_Info(Ttot, Wk_dict, params_setup, T_rec)

    DirRho         = params_setup['DirRho']
    T_Rho[DirRho]  = Ttot
    Wk_Rho[DirRho] = Wk_dict

    if params_setup['StateName'] == 'KapRnd':
        T_Rho2, Wk_Rho2 = Tune_Kappa_Tomography(params_setup, target_density_matrix, \
                input_S, T_rec, Ncpu_Pj, pm_MiFGD, pm_RGD)
        T_Rho.update(T_Rho2)
        Wk_Rho.update(Wk_Rho2)

    return T_Rho, Wk_Rho

def Tune_Kappa_Tomography(params_setup, target_density_matrix, input_S, T_rec, Ncpu_Pj, \
                pm_MiFGD, pm_RGD):
    """ change the singular values and the do tomography over measurements.
        params_setup must have the key of List_Kappa that specifies the ratio between 
            singular values. 
                (eg) List_Kappa = [10, 20, 30]
                     List_alpha = [0.25, 0.5, 1, 2, 3]
 
    Args:
        params_setup (dict): dictionary of parameters that must include key of List_Kappa 
        target_density_matrix (nd.array): density matrix of the target state
        input_S (class instance): circuit implementation of the target state 
        T_rec (dict): dictionary of recording running time
        Ncpu_Pj (int): number of parallel processors  
        pm_MiFGD (dict): dictionary of parameters running the MiFGD optimization method
        pm_RGD (dict): dictionary of parameters running the RGD optimization method 

    Returns:
        dict: (T_Rho) dictionary of recording running time for each case with different optimzation methods
        dict: (Wk_Rho) dictionary of recording the optimzation results for each case with different optimzation methods
    """

    T_Rho  = {}      #  record Ttot for each Rho
    Wk_Rho = {}      #  record Wk_dict for each Rho                    

    Nr        = params_setup['Nr']
    StateName = params_setup['StateName']
    GenNewRand = params_setup['Generate New rand']  

    if StateName != 'KapRnd':
        print('  StateName NOT "KapRnd" --> no need to tune Kappa')
        return T_Rho, Wk_Rho
    
    elif StateName == 'KapRnd':
        List_alpha = params_setup['List_alpha']


        Data_In = params_setup['Data_In']

        if len(Data_In) >0 and Data_In[0] != '':        # DirRho specified
            raise TypeError(' The KapRnd not suitable for specifying DirRho in Data_In')
        elif len(Data_In) > 2:                          # DirMease specified
            raise TypeError(' The KapRnd not suitable for specifying DirMeas')
            

        if GenNewRand == 1:            #  need to create New Rho
            u, s, vh = Decompose_Rho(target_density_matrix, Nr, params_setup)
            del target_density_matrix

        for alpha in List_alpha:

            if GenNewRand == 1:       #  create New Rho
                NewRho, s_new, Kappa = Tune_Kappa(u, s, vh, Nr, alpha)
            elif GenNewRand == 0:     #  Rho already existed

                #raise ValueError(' need to Generate Rho for New Kappa, i.e. need StVer[1] = 1')
                
                NewRho = []
                s_new, Kappa = Gen_sNew_Power(alpha, Nr)

            params_setup, Dt_meas = Params_Kappa_change(params_setup, Kappa, alpha, NewRho, s_new)
            T_rec['measurement'] = Dt_meas

            if GenNewRand == 0:       #  to read existing Rho
                DirRho = params_setup['DirRho']

                with open('{}/RhoKappa.dat'.format(DirRho), 'rb') as f:
                    Kappa2, NewRho = pickle.load(f)

            Ttot, Wk_dict, params_setup, T_rec = \
                Tomography_over_measurement(params_setup, NewRho, \
                                        input_S, pm_RGD, pm_MiFGD, T_rec)
            Rec_Info(Ttot, Wk_dict, params_setup, T_rec)

            DirRho               = params_setup['DirRho']

            Wk_Rho[DirRho]       = Wk_dict
            T_Rho[DirRho]        = Ttot

    return T_Rho, Wk_Rho

def Add_More_Kappa(params_setup, T_rec, Ncpu_Pj, List_alpha, pm_MiFGD, pm_RGD):
    """ To continue generating more kappa caes, i.e. more singular value sets that 
        have different Kappa values

    Args:
        params_setup (dict): dictionary of parameters
        T_rec (dict): dictionary of recording running  
        Ncpu_Pj (int): number of processors for parallel executions
        List_alpha (list): list of alpha parameters representing the new ratio between singular values 
                (eg)        List_Kappa = [40, 50]
                            List_alpha = [0.3, 1, 2]
        pm_MiFGD (dict): dictionary of parameters for running the MiFGD tomography optimization
        pm_RGD (dict): dictionary of parameters for running the RGD tomography optimization

    Raises:
        ValueError: params_setup['Generate New rand'] must = 1 for generating new Kappa cases

    Returns:
        dict: (T_Rho) dictionary of recording running time for each case with different optimzation methods
        dict: (Wk_Rho) dictionary of recording the optimzation results for each case with different optimzation methods
    """

    Dir2m     = params_setup['Dir2m']
    DirRho    = '{}/Kap0'.format(Dir2m)

    GenNewRand = params_setup['Generate New rand']  

    if GenNewRand == 0:
        raise ValueError(' to generate New Kappa --> need to Generate new random matrices')

    T_rec['Ncpu_meas']    =  1   # same as  Get_measurement_by_labels

    target_density_matrix, rho = Gen_randM(params_setup, DirRho)
    input_S = None

    params_setup['List_alpha'] = List_alpha
    T_Rho, Wk_Rho = Tune_Kappa_Tomography(params_setup, target_density_matrix, \
                    input_S, T_rec, Ncpu_Pj, pm_MiFGD, pm_RGD)

    return T_Rho, Wk_Rho



def Sample_Rnd_Tomography(params_setup, Samples, Num_Rho, T_rec, Ncpu_Pj, pm_MiFGD, pm_RGD):
    """ to run several samples of random density matrices cases

    Args:
        params_setup (dict): dictionary of parameters
        Samples (str): sub-directory name for storing different samples
        Num_Rho (int): number of samples to generate random density matrices 
        T_rec (dict): recording the running time 
        Ncpu_Pj (int): number of processors to have parallel executions
        pm_MiFGD (dict): dictionary of parameters for running the MiFGD method 
        pm_RGD (dict): dictionary of parameters for running the RGD method

    Returns:
        dict: updated params_setup dictionary
    """
    StateName   = params_setup['StateName']
    
    if StateName != 'rand':
        print('StateName = {}'.format(StateName))
        print('  -->   NO need to generate random matrices \n\n')
        return params_setup
 
    for numID in range(1, Num_Rho+1):

        params_setup = Params_RND_samples(params_setup, Samples, numID)

        target_density_matrix, input_S, T_rec, Ncpu_meas = \
            Get_measurement_by_labels(params_setup, label_list, T_rec)

        T_Rho, Wk_Rho = Default_Setting_Run_Tomography(params_setup, \
                    target_density_matrix, input_S, T_rec, Ncpu_Pj, pm_MiFGD, pm_RGD)
        
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
    #                                                                                   #
    #       The following two keys only needed for StateName = 'rand' or 'KapRnd'       #
    #                                                                                   #
    #       State_version['rand matrix version'']                                       #
    #                = version of generated random density matrices of rank Nr          #
    #       State_version['Generate New rand']                                          #
    #                =  0    loading random density matrices (generated earlier)        #
    #                           --> will check if the directory exists or not           #
    #                   1    generate new random density matrices of rank Nr            # 
    #                           if file already exists for specified version            #
    #                           then the version will +1 automatically                  #
    # --------------------------------------------------------------------------------- #

    #exec(open('parm_yaml.py').read())

    Choose = 1
    Nk, m, mea = basic_sys_setting(Choose)

    # ---------------------------------------------------------- #
    #   StateName choce = 'GHZ' | 'Had' |  'rand' |  'KapRnd'    #
    # ---------------------------------------------------------- #

    State_version = { 'StateName': 'GHZ',        #  = 'GHZ', 'Had', 'rand', 'KapRnd'
                      'Nr': 1,                    #  rank: 'GHZ', 'Had' --> only allow Nr =1 
                      'Generate New rand': 1,     #  (only for 'rand', 'KapRnd') 1: generate new data |  0: load existing data
                      'rand matrix version': 1,   #  (only for 'rand', 'KapRnd') version for random matrices  
                      'List_alpha': [5, 7, 10]    #  only for 'KapRnd' to tune Kappa (not needed for other states)
    }


    Version_Def = { 
        'Gen New Proj sampling': 1,     #  -1: specify fixed list, 1: generate sampling, 0: load existing sampling
        'Proj version': 1,              #  counting projector version number
        'Gen New Measure': 1,           #  1: generate new shot/calculated measure, 0: load existing 
        'measure version': 2,           #  counting measure version number
    }



    # --------------------------------- #
    #   some default parameters         #
    # --------------------------------- #

    Pj_method = 1               #   the method to save | load  projectors
    mea_method = 1              #   the method to save | load  measurement_dict (count_dict) 

    measure_method = 3          #  = 1: direct label_list,  = 3: parallel cpu
    #measure_method = 1          #  = 1: direct label_list,  = 3: parallel cpu

    # ----------------------------------------------------- #
    #   version parameter for  Projectors & measurement     #
    #                          DirStore to store data       #
    # ----------------------------------------------------- #

    DirStore = './DataTest'
    #Data_In  = []


    # ---------------------------------------------------------- #
    #   (D) start running the optimization program               #
    #                   for tomography                           # 
    #   (D-0) generate the Naming of the directory to store      #
    #       according to the state & version defined in (A) (B)  # 
    # ---------------------------------------------------------- #    

    params_setup = Params_define(State_version, Nk, m, mea,\
                Pj_method, mea_method, measure_method, DirStore) 

    params_setup = State_Naming(params_setup, \
#                        State_version, Version_Def, Data_In)
                        State_version, Version_Def)

    del Nk, m, mea, Pj_method, mea_method, measure_method

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

    Run_meas_Tomography = 1
    if Run_meas_Tomography == 1:      # do measurement & run tomography

        target_density_matrix, input_S, T_rec, Ncpu_meas = \
            Get_measurement_by_labels(params_setup, label_list, T_rec) 

#def fjeifj():

        T_Rho, Wk_Rho = Default_Setting_Run_Tomography(params_setup, \
                        target_density_matrix, input_S, T_rec, Ncpu_Pj, \
                        pm_MiFGD, pm_RGD)

#def fjei():
    # ----------------------------------------- #
    #   special cases:  run more cases          #
    #   [default]    Special_Usage = 0          #
    # ----------------------------------------- #

    Special_Usage = 2

    if Special_Usage == 1:    #  adding more Kappa cases  (only for KapRnd)

        if params_setup['StateName'] == 'KapRnd':   # continue generate more Kappa

            List_alpha = [11, 12]

            T_Rho, Wk_Rho = Add_More_Kappa(params_setup, T_rec, Ncpu_Pj, List_alpha, pm_MiFGD, pm_RGD)


    elif Special_Usage == 2:    # having more samples (only for 'rand')

        Samples = 'sample1'
        Num_Rho  = 5

        params_setup = Sample_Rnd_Tomography(params_setup, Samples, \
                        Num_Rho, T_rec, Ncpu_Pj, pm_MiFGD, pm_RGD)

    print(' ******   Happy Ending   *****')