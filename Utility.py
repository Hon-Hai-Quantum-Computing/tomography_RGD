
import numpy as np
from numpy import linalg as LA

import os

import pickle

import sys
import warnings

warnings.filterwarnings("ignore")

sys.path.append('./quQST')

import projectors
import states
import measurements
import methodsRGD

import qutip as qu
from qutip.states import basis

import Get_param

# -------------------------- #
import time

import multiprocessing
from itertools import repeat
from functools import reduce

from Gen_localPauli import counting_elements

def Params_define(Nk, StateName, num_labels, num_shots, Nr, \
                Pj_method, mea_method, measure_method, DirStore):

    Name = '{}-{}'.format(StateName, Nk)
    if StateName == 'rand':
        Name = '{}-r{}'.format(Name, Nr)
    elif StateName == 'KapRnd':
        Name = '{}-r{}'.format(Name, Nr)

    Dir0 = '{}/{}'.format(DirStore, Name)

    params_setup = {'n': Nk,
                   'Name': Name,
                   'StateName': StateName, 
                   'num_labels': num_labels,
                   'num_shots': num_shots,
                   'Nr': Nr,
                   #'version': version,
                   #'StVer': StVer,
                   #'Noise': Noise,
                   'Dir0': Dir0,
                   #'Dir': Dir,
                   #'measurement_store_path': Dir_meas,
                   #'projector_store_path': Dir_proj,
                   'Pj_method': Pj_method,
                   'mea_method': mea_method,
                   'measure_method': measure_method
                   #'Obtain_Pj': Obtain_Pj
                   }

    return params_setup

def Params_RND_samples(params_setup, Samples, numID):
    """  (eg)      Samples = 'samples3'
    """

    Dir0     = params_setup['Dir0']
    Name     = params_setup['Name']
    Dir_proj = params_setup['projector_store_path']
    PjFile   = Dir_proj.split("/")[-1]

    Dir0     = '{}/{}'.format(Dir0, Samples)

    DirRho   = '{}/RND{}'.format(Dir0, numID)
    Dir2m    = '{}/{}'.format(DirRho, PjFile)
    Dir_meas = '{}_zExact'.format(Dir2m)
    
    print('DirRho   = {}'.format(DirRho))
    print('Dir_meas = {}'.format(Dir_meas))

    if not os.path.exists(DirRho):
        os.makedirs(DirRho)

    params_setup['DirRho']                 = DirRho
    params_setup['measurement_store_path'] = Dir_meas
    params_setup['Dir2m']                  = Dir2m

    return params_setup

def Params_Kappa_change(params_setup, Kappa, alpha, NewRho, s_new):

    tm0 = time.time()

    Dir2m       = params_setup['Dir2m']
    mea_sub_pth = params_setup['measure_sub_path']
    
    DirRho  = '{}/Kap{}_alp{}'.format(Dir2m, Kappa, alpha)
    #DirRho  = '{}/alp{}_Kap{}'.format(Dir2m, alpha, Kappa)

    Dir_mea = '{}/{}'.format(DirRho, mea_sub_pth)

    params_setup['Kappa']  = Kappa
    params_setup['DirRho'] = DirRho
    params_setup['measurement_store_path'] = Dir_mea

    print(' Dir2m     = {}'.format(Dir2m))
    print(' DirRho    = {}'.format(params_setup['DirRho']))
    print(' meas_path = {}'.format(params_setup['measurement_store_path']))

    if not os.path.exists(DirRho):
        os.makedirs(DirRho)

    # ------------------------------------ #
    #   record NewRho / s_new info         #
    # ------------------------------------ #

    Nk       = params_setup['n']
    m        = params_setup['num_labels']
    Dir_proj = params_setup['projector_store_path']

    StVer    = params_setup['StVer']

    if StVer[1] == 0:       #  reading existing Rho
        print(' New Rho with with New Kappa alrady existed!')

        with open('{}/RhoKappa.dat'.format(DirRho), 'rb') as f:
            Kappa2, NewRho = pickle.load(f)
    
        if not np.isclose(Kappa, Kappa2):
            raise ValueError(' Kappa Not consistent in the DirRho')

    elif StVer[1] == 1:     #  recording new Rho & s_new

        with open('{}/RhoKappa.dat'.format(DirRho), 'wb') as f:
            pickle.dump([Kappa, NewRho], f)        

        F_rec = '{}/s_new.txt'.format(DirRho)
        f = open(F_rec, 'a')

        f.write(' ----------------------------------\n')
        f.write(' DirRho   = {}\n'.format(DirRho))
        f.write(' Dir_proj = {}\n\n'.format(Dir_proj))

        f.write('     Nk       = {}\n'.format(Nk))
        f.write('     m        = {}\n\n'.format(m))
        
        f.write(' sum(s_new) = {}\n'.format(sum(s_new)))

        f.write(' s_new  = {}\n'.format(s_new))
        
        RatioMinMax = max(s_new) / min(s_new)
        RatioHdTL   = s_new[0]   / s_new[-1]

        f.write('  max(s_new) = {}\n'.format(max(s_new)))
        f.write('  min(s_new) = {}\n\n'.format(min(s_new)))

        f.write(' RatioMinMax = max(s_new) / min(s_new) = {}\n'.format(RatioMinMax))
        f.write(' RatioHdTL   = s_new[0]   / s_new[-1]  = {}\n'.format(RatioHdTL))
        f.close()

    # ----------------------------- #
    #   obtain exact coef           #
    # ----------------------------- #
    
    New_Pj_shot = params_setup['New_Pj_shot']

    if New_Pj_shot[1] == 1:     #  do new measurement
        s_label, yProj, zN, yProj_Exact = \
            state_S_coef(params_setup, NewRho)      #  directly create


    tm1 = time.time()
    dt_meas = tm1 - tm0

    return params_setup, dt_meas


def Check_Dir_version(Dir0, StVer, Appendix=''):
        
    print(' StVer = {}'.format(StVer))
    
    ver_Data = StVer[0]
    Dir_Data = '{}/R{}{}'.format(Dir0, ver_Data, Appendix)

    Fname = '{}/Rho.dat'.format(Dir_Data)
    print('  is Fname = {}  exists?  {}'.format(Fname, os.path.isfile(Fname)))

    if StVer[1] == 0:             #  loading (state | DM) data

        if not os.path.isfile(Fname):
            print('  File NOT exists')
            return
        
        Dir0 = Dir_Data

    elif StVer[1] == 1:           #  generate (state | DM) data

        while os.path.isfile(Fname):
            ver_Data += 1
            Dir_Data = '{}/R{}{}'.format(Dir0, ver_Data, Appendix)
            Fname = '{}/Rho.dat'.format(Dir_Data)
            print('  Dir_Data = {}'.format(Dir_Data))            
        StVer = [ver_Data, StVer[1]]

        if not os.path.exists(Dir_Data):
            os.makedirs(Dir_Data)
        Dir0 = Dir_Data

    print(' Dir0 = {}'.format(Dir0))


    return Dir0, StVer

def Ver_Rho_Naming(params_setup, StVer, GetDir0=0):

    StateName = params_setup['StateName']

    StateName_List = ['GHZ', 'Had', 'rand', 'KapRnd', 'quWst']

    if StateName not in StateName_List:
        raise ValueError(' Specified StateName not implemented')

    # ----------------------------- #
    #   creating Dir0               #
    # ----------------------------- #

    Dir0 = params_setup['Dir0']

    if not os.path.exists(Dir0):
        os.makedirs(Dir0)

    #if GetDir0 == 1:
    #    return Dir0

    # ----------------------------- #
    #   check if existing DataRho   #
    # ----------------------------- #

    Data_In   = params_setup['Data_In']

    if len(Data_In) !=0 and Data_In[0] != '':  # existing DataRho specified in Data_In[0]

        DirRho   = Data_In[0]

        params_setup['DirRho']        = DirRho  
        params_setup['StVer']         = StVer

        if not os.path.exists(DirRho):
            os.makedirs(DirRho)

        return

    # ----------------------------------------- #
    #   extra treatment for random states       #
    #                                           #
    #   to generated or loading random state    #
    # ----------------------------------------- #

    if StateName == 'KapRnd':           #  to change the kappa

        #Dir2p  = Dir0                    #  to convert to proj later
        #DirRho = Check_Dir_version(Dir0, StVer)     #  no change Dir0
        #DirRho = Check_Dir_version(Dir0, StVer, '_qutip')  #  no change Dir0
        DirRho, StVer = Check_Dir_version(Dir0, StVer, '/Kap0')  #  no change Dir0

    elif StateName == 'rand':             #  random matrices
        
        DirRho, StVer = Check_Dir_version(Dir0, StVer)
        #Dir2p         = DirRho

    else:                               #  pure states (GHZ, Had, quWst)
        DirRho = Dir0           #  will be changed in later Ver_Proj_Naming
        #Dir2p  = Dir0
        

    params_setup['Dir0']   = Dir0
    params_setup['DirRho'] = DirRho     #  may be changed in later Ver_Proj_Naming
    params_setup['StVer']  = StVer


def Ver_Proj_Naming(params_setup, version):
    """ version[0]     = the version of the projector/Pauli matrices

        New_Pj_shot[0] = 1 -> to create new projectors/Pauli matrices 
    """
    StateName  = params_setup['StateName']
    Name       = params_setup['Name']
    num_labels = params_setup['num_labels'] 
    Pj_method  = params_setup['Pj_method']

    Data_In    = params_setup['Data_In']
    DirRho     = params_setup['DirRho']
    Dir0       = params_setup['Dir0']

    StVer      = params_setup['StVer']

    if StateName == 'KapRnd':           #  to change the kappa
        Dir2p  = Dir0                   #  to convert to proj later
    elif StateName == 'rand':           #  random matrices
        Dir2p  = DirRho
    else:                               #  pure states (GHZ, Had)
        Dir2p  = Dir0
    params_setup['Dir2p']    =  Dir2p

    # --------------------------------- #
    #   check if existing read in data  #
    # --------------------------------- #
    if len(Data_In) >0 and Data_In[1] != '':          #  specified DirProj

        Dir_proj = Data_In[1]

        xx = Dir_proj.split('/')[-1].split('_')
        if xx[-1] == 'Proj':
            Dp = xx[-2][1:]
        elif xx[-1][0] == 's':
            Dp = xx[-1][1:]
        else:
            #print('xx = {}'.format(xx))

            Dp = 0
            print(' Not defined! --> Let version[0] = Ver_Proj = {}'.format(Dp))
        version[0] = int(Dp)

        Noise    = version[2]

        params_setup['projector_store_path'] = Dir_proj
        params_setup['version']              = version
        params_setup['Noise']                = Noise 

        return

    # ----------------------------- #
    #   deal with projectors        #
    # ----------------------------- #

    New_Pj_shot = params_setup['New_Pj_shot']


    if Pj_method == 0:      #  original MiFGD method
        Pj_Name = 'projectors'
    elif Pj_method == 1:    #  new method
        Pj_Name = 'Proj'

    ver_Prj = version[0]

    if StateName == 'KapRnd':
        Dir_proj = '{}/{}_m{}_s{}/'.format(Dir2p, Pj_Name, num_labels, ver_Prj)
        #Dir2m = DirRho
        Dir2m = '{}/R{}'.format(Dir2p, StVer[0])

    elif StateName == 'rand':  
        Dir2m = '{}/{}_m{}_s{}/{}_m{}_s{}'.format(DirRho, Name, num_labels, ver_Prj, 
                                                          Name, num_labels, ver_Prj)
        Dir_proj = '{}_{}'.format(Dir2m, Pj_Name)

    else:                      #  pure state (GHZ, Had) -> same as 'rand'
        Dir2m = '{}/{}_m{}_s{}/{}_m{}_s{}'.format(Dir2p, Name, num_labels, ver_Prj, 
                                                         Name, num_labels, ver_Prj)
        Dir_proj = '{}_{}'.format(Dir2m, Pj_Name)
        DirRhoNew = '{}/{}_m{}_s{}'.format(DirRho, Name, num_labels, ver_Prj)


    if New_Pj_shot[0] == 1:                 # i.e. Create New Projectors
        while os.path.exists(Dir_proj):
            ver_Prj = ver_Prj + 1

            if StateName == 'KapRnd':
                Dir_proj = '{}/{}_m{}_s{}/'.format(Dir2p, Pj_Name, num_labels, ver_Prj)
                #Dir2m = DirRho

            elif StateName == 'rand':
                Dir2m = '{}/{}_m{}_s{}/{}_m{}_s{}'.format(DirRho, Name, num_labels, ver_Prj, 
                                                                  Name, num_labels, ver_Prj)
                Dir_proj = '{}_{}'.format(Dir2m, Pj_Name)

            else:           #  pure state (GHZ, Had) -> same as 'rand'
                Dir2m = '{}/{}_m{}_s{}/{}_m{}_s{}'.format(Dir2p, Name, num_labels, ver_Prj,
                                                                 Name, num_labels, ver_Prj)
                Dir_proj = '{}_{}'.format(Dir2m, Pj_Name)
            
                DirRhoNew = '{}/{}_m{}_s{}'.format(DirRho, Name, num_labels, ver_Prj)


            print(' version (of projectors) = {}'.format(ver_Prj))
            version[0] = ver_Prj

    Noise = version[2]

    params_setup['Noise']                = Noise 
    params_setup['version']              = version
    #params_setup['Dir']                 = Dir
    params_setup['projector_store_path'] = Dir_proj

    if StateName == 'GHZ' or StateName == 'Had' or StateName == 'quWst':
        params_setup['DirRho'] = DirRhoNew     # important


def Ver_meas_Naming(params_setup):
    """ version[1]     = the version of the measurements

        New_Pj_shot[0] = 1 -> to create new measurements
    """

    StateName  = params_setup['StateName']
    Name       = params_setup['Name']

    num_labels = params_setup['num_labels']
    num_shots  = params_setup['num_shots']
    mea_method = params_setup['mea_method']

    New_Pj_shot = params_setup['New_Pj_shot']

    version     = params_setup['version']
    ver_Prj     = version[0]      # = version of projectors/Pauli matrices

    StVer       = params_setup['StVer']

    # ---------------------- #
    #   specify the Dir2m    #
    # ---------------------- #
    DirRho  = params_setup['DirRho']
    Dir2p   = params_setup['Dir2p']  

    if StateName == 'KapRnd':
        Dir2m = '{}/R{}'.format(Dir2p, StVer[0])

    elif StateName == 'rand':  
        Dir2m = '{}/{}_m{}_s{}/{}_m{}_s{}'.format(DirRho, Name, num_labels, ver_Prj, 
                                                          Name, num_labels, ver_Prj)
    else:                      #  pure state (GHZ, Had) 

        Dir2m    = '{}/{}_m{}_s{}'.format(DirRho, Name, num_labels, ver_Prj)

    params_setup['Dir2m']              = Dir2m

    # ------------------------------ #
    #   check if existing DirMeas    #
    # ------------------------------ #

    Data_In = params_setup['Data_In']

    if len(Data_In) > 2:
        Dir_meas = Data_In[2]
        Dir2m    = '{}/{}_m{}_s{}'.format(Dir_meas, Name, num_labels, ver_Prj)

        params_setup['measurement_store_path'] = Dir_meas
        params_setup['Dir2m']                  = Dir2m
        return
    
    # ----------------------------- #
    #   deal with shots or noise    #
    # ----------------------------- #

    zNoise   = version[2]      #  specify the noise model (0: exact)

    if StateName == 'KapRnd':
        print('    StateName = KapRnd')

        ver_meas = version[1]      #  version for this model

        if zNoise == 0:         #  i.e.  exact case
            zText = 'zExact'    #  noise model text
        #Dir_meas = 'zN{}_v{}_Exact'.format(zNoise, ver_Mode1)
            
        Dir_meas = 'Pm{}_s{}_{}_v{}'.format(num_labels, ver_Prj, zText, ver_meas)

        if New_Pj_shot[1] == 1:        #  create new measurement
            while os.path.exists('{}/R0/{}'.format(Dir2m, Dir_meas)):
                
                ver_meas = ver_meas + 1
                Dir_meas = 'Pm{}_s{}_{}_v{}'.format(num_labels, ver_Prj, zText, ver_meas)

                version = [ver_Prj, ver_meas, version[2]]
                print("     version now = {}".format(version))

            version = [ver_Prj, ver_meas, version[2]]

    elif zNoise == 0:                 #  the exact coef for 'rand'
        print(' -----  exact value for the noise model  --------- ')

        ver_meas = version[1]         #  version for this measurement

        Dir_meas = '{}_zN{}_v{}_measurements'.format(Dir2m, zNoise, ver_meas)

        print('  Dir_meas  = {}'.format(Dir_meas))

        if New_Pj_shot[1] == 1:        #  create new measurement
            while os.path.exists(Dir_meas):
                
                ver_meas  = ver_meas + 1
                Dir_meas  = '{}_zN{}_v{}_measurements'.format(Dir2m, zNoise, ver_meas)

                #version = [ver_Prj, ver_zMode1, version[2]]
                version = [ver_Prj, ver_meas, zNoise]

                print("     version now = {}".format(version))


            version = [ver_Prj, ver_meas, zNoise]


    elif zNoise == -1:         # (GHZ or Had) using shot measurements from qiskit

        if mea_method == 0:        #  original method
            MeasureName = 'measurements'
        elif mea_method == 1:      #  new method
            MeasureName = 'Measure'

        ver_Shot = version[1]      # = version of measurements 
        Dir_meas = '{}_shot{}_v{}_{}'.format(Dir2m, num_shots, ver_Shot, MeasureName)

        if New_Pj_shot[1] == 1:

            #while os.path.exists(DirQiskit):
            #    version = version + 1
            #    Dir = '{}/{}_m{}_shot{}_v{}'.format(Dir0, Name, num_labels, num_shots, version)
            #    DirQiskit = '{}_qiskit.dat'.format(Dir)
            #    print("version now = {}".format(version))
        
            while os.path.exists(Dir_meas):
                ver_Shot = ver_Shot + 1
                Dir_meas = '{}_shot{}_v{}_{}'.format(Dir2m, num_shots, ver_Shot, MeasureName)
                version = [ver_Prj, ver_Shot, zNoise]
                print("version now = {}".format(version))
    
            version = [ver_Prj, ver_Shot, zNoise]


    print('         version = {}'.format(version))
    params_setup['version']            =  version

    if StateName == 'KapRnd':
        params_setup['measurement_store_path'] = '{}/{}'.format(DirRho, Dir_meas)
        params_setup['measure_sub_path']       = Dir_meas
    else:
        params_setup['measurement_store_path'] = Dir_meas


    #return version, Dir_meas 


def State_Naming(StVer, version, \
        params_setup, New_Pj_shot, Data_In=[]):

    # --------------------------------------- #
    #   check if special New_Pj_shot          #
    # --------------------------------------- #
    Nk         = params_setup['n']          # number of qubits
    num_labels = params_setup['num_labels']

    if New_Pj_shot[0] == -1:   # specify Pauli list, instead of sampling

        sym_Pauli        = ['X', 'Y', 'Z']
        NumTerm_Pauli    = counting_elements(Nk, sym_Pauli)
        num_labels_count = NumTerm_Pauli[1] + NumTerm_Pauli[2]

        if num_labels != num_labels_count:
            #raise ValueError('   num_labels not specified correctly')
        
            print(' ****************************************************** ')
            print(' ***   New_Pj_shot[0] == -1  --> specify Pauli list  ** ')
            print(' ***   --> num_labels is calculated and changed      ** ')
            print(' ****************************************************** ')

            params_setup['num_labels'] = num_labels_count

    # ----------------------------- #
    #   if specified Dir & PjFile   #
    # ----------------------------- #

    params_setup['New_Pj_shot'] = New_Pj_shot
    params_setup['Data_In']     = Data_In

    Ver_Rho_Naming(params_setup, StVer)         # DirRho

    Ver_Proj_Naming(params_setup, version)      # DirProj

    Ver_meas_Naming(params_setup)               # DirMeas

    #return Name, params_setup, version, Dir_proj, Dir_meas, DirRho, StVer
    return params_setup



def Gen_Rho(Nk, StateName):
    """  (input)
            version to record which version we are recording
    """

    if StateName=='GHZ':             #  StateName = 'GHZ'    
        target = states.GHZState(Nk)        # Nk = num_qubits

    elif StateName=='Had':           #  StateName = 'Had'          
        target = states.HadamardState(Nk)   #  states from states.py

    target.create_circuit()
    target_density_matrix = target.get_state_matrix()

    rho = None
    return target, target_density_matrix, rho


def Gen_Wst(Nk, Dir0):
    '''     Generate W state
    '''
    pick = [0]*Nk    
    pick[0] = 1 
    Wst = basis([2]*Nk, pick)
    print(' pick = {}  ------'.format(pick))    
    #print('  Wst = {}'.format(Wst))

    for ii in range(1, Nk):
        pick = [0] * Nk
        pick[ii] = 1

        st = basis([2]*Nk, pick)
        Wst = Wst + st

        print(' pick = {}  ------'.format(pick))    
        #print(' st  = {}'.format(st))    
        #print(' Wst = {}'.format(Wst))

    Wst  = Wst /  np.sqrt(Nk)
    norm = LA.norm(Wst)

    #print(' Wst = {}'.format(Wst))
    print(' norm(Wst) = {}'.format(LA.norm(Wst)))

    if np.isclose(norm, 1.0):
        print('  the normalization is correct')
    else:
        print('  sth wrong in the normalization')
        return
    
    Wst = np.array(Wst)

    RhoWst = Wst @ Wst.T.conj()

    # ------------------------- #
    #       Dir Naming          #
    # ------------------------- #
    Fname = '{}/Rho.dat'.format(Dir0)

    print(' Fname = {}'.format(Fname))
    print(' is file exists?  {}'.format(os.path.isfile(Fname)))

    with open(Fname, 'wb') as f:
        pickle.dump([Wst, RhoWst], f)

    return Wst, RhoWst


def Gen_randM(Nk, StateName, Nr, StVer, Dir0):
    """
        Nr = 1                  #  pure state -->  rank = 1    
    """

    Fname = '{}/Rho.dat'.format(Dir0)

    if StateName=='rand' or StateName == 'KapRnd':  #  randomized matrix

        if StVer[1] == 1:           #  to generate data
            rho = qu.rand_dm_ginibre(2**Nk, dims=[[2]*Nk , [2]*Nk], rank=Nr)
            target_density_matrix = rho.full()

            with open(Fname, 'wb') as f:
                pickle.dump([rho, target_density_matrix], f)

        elif StVer[1] == 0:         #  loading data

            print(' Fname = {}'.format(Fname))
            print(' is file exists?  {}'.format(os.path.isfile(Fname)))
            
            with open(Fname, 'rb') as f:
                rho, target_density_matrix = pickle.load(f)

    return target_density_matrix, rho

def Decompose_Rho(Rho, Nr, params_setup):
    u, s, vh = LA.svd(Rho, hermitian=True)

    u  = u[:, :Nr]
    s  = s[:Nr]
    vh = vh[:Nr, :]

    NewRho = u @ np.diag(s) @ vh

    if not np.isclose(np.sum(s), 1.0):
        print('\n\n **** the sum of singular values is not ONE ***** \n\n')
        return
    if not np.allclose(NewRho, Rho):
        print('\n\n ***** SVD not correct for Rho  ***** \n\n')
        return

    print('\n\n     np.sum(s)  =  1    -->  correct   ')
    print('      u @ np.diag(s) @ vh  = Rho    -->  correct \n\n')

    DirRho = params_setup['DirRho']
    print('  DirRho = {}'.format(DirRho))

    f = open('{}/sR0.txt'.format(DirRho), 'a')
    f.write(' DirRho = {}\n'.format(DirRho))
    f.write(' sR0 = {}\n'.format(s))
    f.write(' Kappa = {}\n'.format(s[0]/s[-1]))
    f.close()

    return u, s, vh


def Gen_sNew_Power(alpha, Nr):
    # --------------------------------------------------------- #
    #   implement: s_new =  [1, alpha, alpha^2, ... alpha^Nr]   #
    #       normalized to  sum(s_new) = 1                       #
    # --------------------------------------------------------- #

    power = np.arange(Nr)

    s_new = np.ones(Nr) * alpha
    s_new = s_new ** power          # = [1, alpha, ..., alpha^Nr]

    Norm  = sum(s_new)
    Norm2 = (1.0 - alpha**Nr) / (1.0 - alpha)
    
    if np.isclose(Norm, Norm2):
        print('  norm equal to analytical formula')
        print('  norm  =  {}'.format(Norm))
    else:
        print('\n  ***  WRONG in assigning s_new  ***')
        print('       Norm = {}  !=  Norm2 = {}\n'.format(Norm, Norm2))
        return 
    print('  s_new = {}'.format(s_new))

    s_new = s_new / Norm
    print('  sum(s_new) = {}'.format(sum(s_new)))
    
    if np.isclose(sum(s_new), 1.0):
        if alpha > 1.0:
            s_new = s_new[::-1]

        Kappa = s_new[0] /  s_new[-1]

        print('  -->  sum   = 1.0  as a prob distribution ')
        print('       s_new = {}\n'.format(s_new))

        Kappa2 = alpha**(Nr-1)
        if alpha < 1.0:
            Kappa2 = 1/Kappa2
        print('  -->  Kappa = {}  =?  {}'.format(Kappa, Kappa2))

        return s_new, Kappa

    else:
        print('\n ***    sum  !=  1.0   -->  need check    ***\n')
        return 

def Tune_Kappa(u, s, vh, Nr, alpha):

    s_new, Kappa = Gen_sNew_Power(alpha, Nr)

    print('   s_new = {}'.format(s_new))
    print('   Kappa = {}\n'.format(Kappa))

    # ----------------- #
    #   new Rho         #
    # ----------------- #
    
    NewRho = u @ np.diag(s_new) @ vh
    return NewRho, s_new, Kappa


# ------------------------- #
#   for parallel CPU usage  #
# ------------------------- #

def split_list(x, num_parts):
	n = len(x)
	size = n // num_parts
	parts = [x[i * size: (i+1) * size] for i in range(num_parts - 1 )]
	parts.append(x[(num_parts - 1) * size:])
	return parts

# --------------------------------------------- #
#   for paralle calc qiskit state measurement   #
# --------------------------------------------- #
def mp_state_measure(params_setup, label_list):
    params_mp = {}
    params_mp['n']         = params_setup['n']
    params_mp['StateName'] = params_setup['StateName']
    params_mp['num_shots'] = params_setup['num_shots']

    num_CPU = multiprocessing.cpu_count()
    if num_CPU > len(label_list):
        num_CPU = 3

    ID_list         = np.arange(num_CPU)
    label_part      = split_list(label_list, num_CPU)
    print('label_part = {}'.format(label_part))

    tt3 = time.time()

    print('  *********   start  doing parallel qiskit measurement  ******* \n')
    pool = multiprocessing.Pool(num_CPU)
    ID_data_list = pool.starmap(state_measure_wID, zip(ID_list, repeat(params_mp), label_part))

    print(ID_data_list)
    pool.close()
    #pool.join()
    print('  *********   END of qiskit parallel measurement  ******** ')

    tt4 = time.time()
    print('\n  *****  starmap parallel time = {}  ***** \n'.format(tt4- tt3))

    cnt_list, measure_part, data_part_list, bEnd = zip(*ID_data_list)
    data_dict_list_ALL = reduce(lambda x, y: x+y, data_part_list)            
    backend = bEnd[0]

    measure_ALL = {}
    for ii in ID_list:
        measure_ALL.update(measure_part[ii])

    return measure_ALL, data_dict_list_ALL, backend                 


# ------------------------- #
#   state measurement       #
# ------------------------- #

def state_measure_save(Dir, params_setup, labels, 
                       measurement_dict, data_dict_list, backend,
                       state, target_density_matrix, rho, Name=None):

    num_shots  = params_setup['num_shots']
    version    = params_setup['version']
    ver_Shot   = version[1]

    if Name == None:
        Fname = '{}_shot{}_v{}_qiskit.dat'.format(Dir, num_shots, ver_Shot)
    else:
        Fname = '{}_shot{}_v{}_qiskit_{}.dat'.format(Dir, num_shots, ver_Shot, Name)

    #while not os.path.exists(Fname):
    #    with open(Fname, 'wb') as f:
    #        pickle.dump([measurement_dict, data_dict_list, \
    #                     labels, params_setup, backend, target_density_matrix, state, rho], f)

    if os.path.exists(Fname):
        os.remove(Fname)

    with open(Fname, 'wb') as f:
        pickle.dump([measurement_dict, data_dict_list, \
                    labels, params_setup, backend, target_density_matrix, state, rho], f)
        
def state_measure_wID_wrapper(args):
    ID, params_mp, labels = args
    return state_measure_wID(ID, params_mp, labels)



def state_measure_wID(ID, params_mp, labels):
    Nk         = params_mp['n']
    StateName  = params_mp['StateName']
    num_shots  = params_mp['num_shots']

    input_S, target_density_matrix, rho = Gen_Rho(Nk, StateName)

    measurement_dict, data_dict_list, backend = state_measure(num_shots, labels, input_S)

    return [ID, measurement_dict, data_dict_list, backend]


def state_measure(num_shots, labels, state):
    """  (input)
            version to record which version we are recording
         (eg)
                num_qubits = 3
                num_labels = 50
                num_shots = 100
    """

    print('      *****  start state_measure for #labels={}\n'.format(len(labels)))

    backend    = 'qasm_simulator'                  ## Circuits ##

    ## Measurements ##
    ##
    state.create_circuit()
    data_dict_list = state.execute_measurement_circuits(labels,
                                                        backend=backend,
                                                        num_shots=num_shots)
    print('      ******  Already qiskit DONE creating data_dict_list  **** \n')

    measurement_dict = {item['label']: item['count_dict'] for item in data_dict_list}
    print('      ------  DONE converting data_dict_list INTO measurement_dict -----\n')

    return measurement_dict, data_dict_list, backend


def state_S_coef(params_setup, target_density_matrix, rho=[]):
    """  (input)
            version to record which version we are recording
         (eg)
                num_qubits = 3
                num_labels = 50
                num_shots = 100

         rho is only for comparison in method = 1  (no longer used now)       
    """

    #label_list = projectors.generate_random_label_list(m, Nk)
    #label_list = ['ZZZ', 'XYY', 'YXZ', 'ZXX', 'III']
    #labels = ['YYY', 'YXI', 'XII', 'XXY', 'ZXY']

    meas_path = params_setup['measurement_store_path']
    proj_path = params_setup['projector_store_path']
    Pj_method = params_setup['Pj_method']

    print('Pj_method = {}'.format(Pj_method))

    if Pj_method == 0:          #   original method, saving each Pj separately
        labels = [fname.split('.')[0] for fname in os.listdir(proj_path)]    
    #elif 'Pj_method' == 1:        #   new method
    #    label_file = '{}/labels.pickle'.format(proj_path)
    #    with open(label_file, 'rb') as f:
    #        labels = pickle.load(f)
    #    print('label_file = {}'.format(label_file))
    #print('labels = {}'.format(labels))


    ## Measurements from direct trace of Pauli S ##
    ##

    Nk    = params_setup['n']              #  =  num_qubits
    m     = params_setup['num_labels']     #  =  num_labels
    coefI = np.sqrt(m / 2**Nk)

    method = 2
    if method == 2:
        # ------------------------------------------------------ #
        #        loading projectors --> get coef from trace      #
        # ------------------------------------------------------ #
        if Pj_method == 0:          #  original method
            projector_dict = projectors.ProjectorStore.load(proj_path, labels[0:m])
        elif Pj_method == 1:        #  new method
            labels = projectors.ProjectorStore.load_labels_from_file(proj_path)

            #projector_dict = projectors.ProjectorStore.load_PoolMap(proj_path)
            projector_dict = projectors.ProjectorStore.load_PoolMap(proj_path, labels)

        projector_list = [projector_dict[label] for label in labels]

        #yAmea = methodsRGD.Amea(projector_list, rho, m, 2**Nk) # if rho --> ValueError: cannot set WRITEABLE flag to True of this array
        #yProj_Exact = yAmea * coefI     

        yProj_Exact = methodsRGD.Amea(projector_list, target_density_matrix, m, 1)        #   argv[-1] = coef


    # --------------------------- #
    #   generate nosie  N         #
    # --------------------------- #
    Noise = params_setup['Noise']

    #Model = Noise[0]
    zModel = Noise           #  noise model
    if Noise == 0:           #  exact coef case --> no noise
        zN = np.zeros(m)


    # ------------------------------------ #
    #   get result from exact + noise N    #  
    # ------------------------------------ #
    yProj = yProj_Exact + zN
    print(' Noise       = {}'.format(Noise))
    print(' yProj[:10]  = {}'.format(yProj[:10]))        

    # --------------------------------- #
    #   save the results                #
    # --------------------------------- #

    #mea      = params_setup['num_shots']
    ver_meas = params_setup['version'][1]

    #Fname1 = '{}_zN{}_v{}_Noise'.format(Dir, Model, ver_meas)
    #Fname2 = '{}_zN{}_v{}_measurements'.format(Dir, Model, ver_meas)

    Fname1 = '{}/zN{}_v{}_Noise'.format(meas_path, zModel, ver_meas)
    Fname2 = '{}/zN{}_v{}_measurements'.format(meas_path, zModel, ver_meas)


    if not os.path.exists(meas_path):
        os.makedirs(meas_path)


    print(' Fname1 = {}'.format(Fname1))
    print(' Fname2 = {}'.format(Fname2))

    if os.path.exists(Fname1):
        os.remove(Fname1)

    while not os.path.exists(Fname1):
        with open(Fname1, 'wb') as f:
            pickle.dump([yProj_Exact, zN, yProj, Noise, \
                    labels, params_setup, target_density_matrix, rho], f)

        with open(Fname2, 'wb') as f:
            pickle.dump([labels, yProj, zN], f)


    #return np.asarray(labels), coefS_list, yProj
    return np.asarray(labels), yProj, zN, yProj_Exact

    
# --------------------------------- #
#   check the reconstructed error   #
# --------------------------------- #
def check_reconstruct_error(X):
    """ (example input)
    #X = target_density_matrix
    X = rhoM    
    """

    u0, s0, v0h = np.linalg.svd(X, full_matrices=True) 

    u, s, v = methodsRGD.hard(X, r=1)

    M2 = u @ s @ v.T.conj()

    print(' X type = {}'.format(type(X)))
    if isinstance(X, qu.qobj.Qobj):
        print('   reconstructed ele max error = {}'.format(np.max(M2 - X.full())))
    elif isinstance(X, np.ndarray):
        print('   reconstructed ele max error = {}'.format(np.max(M2 - X)))
    print(' ************************************** \n\n')



# ----------------------------------------- #
#           below the old version           #
# ----------------------------------------- #

# ------------------------------------------------- #
#              to compare results                   #
# ------------------------------------------------- #

def Load_qiskit(Dir, params_setup):
    # ------------------------- #
    #   to load qiskit          #
    # ------------------------- #

    mea      = params_setup['num_shots']
    ver_Shot = params_setup['version'][1] 

    Fname = '{}_shot{}_v{}_qiskit.dat'.format(Dir, mea, ver_Shot)

    print(' to load qiskit measurement data from {}'.format(Fname))
    with open(Fname, 'rb') as f:
        measurement_dict, data_dict_list, s_label, params_setup, backend, \
        target_density_matrix, input_S, rho = pickle.load(f)

    return s_label, rho, measurement_dict

def Get_measurement_list(measurement_store_path):

    #measurement_store_path = '{}_measurement'.format(Dir)
    #measurement_store_path = params_setup['measurement_store_path']

    label_list  = measurements.MeasurementStore.load_labels(measurement_store_path)

    start = 0
    end = m
    measurement_dict = measurements.MeasurementStore.load(measurement_store_path, label_list[start:end])

    # ---------------------------------------- #
    #   from  measurement_dict                 #
    #    to   convert into measurement_list    #
    # ---------------------------------------- #

    print(label_list)
    count_dict_list = [measurement_dict[label] for label in label_list]
    measurement_object_list = [measurements.Measurement(label, count_dict) for (label, count_dict) in 
							   zip(*[label_list, count_dict_list])]

    parity_flavor='effective'
    beta = None
    measurement_list = [measurement_object.get_pauli_correlation_measurement(beta, parity_flavor)[label] for
						(label, measurement_object) in zip(*[label_list,
															measurement_object_list])]

    measurement_list = np.asarray(measurement_list)

    return label_list, measurement_list


# -------------------------------------------- #
#      to deal with the Identity label         #
# -------------------------------------------- #

def Delete_Iden(label_list, Nk, m):
    print(' ******   to delete Identity from label_list   *****')

    Iden = ''.join(['I' for i in range(Nk)])

    # ------------------------------------- #
    #   drop the identity  Iden             #
    # ------------------------------------- #

    Go_Del = 1
    while Go_Del == 1:
        try:
            label_list.remove(Iden);    L1 = len(label_list)
        except:
            L1     = len(label_list)
            Go_Del = 0
    
    print('L1 = {}, label_list = {}'.format(L1, label_list))


    while(len(label_list)< m):
        label_add = projectors.generate_random_label_list(3, Nk)

        try:
            label_add.remove(Iden)
            Del_Iden = 1
        except:
            Del_Iden = 0

        # ------------------------------------- #
        #   pad up label_list to num_labels     #
        # ------------------------------------- #
        label_list = label_list + label_add

        print('Del_Iden = {} -> len ={}'.format(Del_Iden, len(label_list)))


    label_list = label_list[0:m]

    print('  New label_list  =  {}'.format(label_list))
    return label_list

def To_Find_Iden_Idx(label_list, Nk):
    # ------------------------------------------------- #
    #   check if there exists identity I in label_list  #
    # ------------------------------------------------- #
    print(' *******    to find out the Identity index    ********')

    #FirstI = [label[0] == 'I' for label in label_list]
    #print('FirstI = {}'.format(FirstI))

    Iden = ''.join(['I' for i in range(Nk)])

    # ----------------------------------------------- #
    #       this .index() can only find one time      #
    # ----------------------------------------------- #
    try:
        ndx = label_list.index(Iden)
    except:
        ndx = -1
    print(ndx)

    Idx = np.where(np.array(label_list) == Iden)[0]         #  where is the Identity

    print('Idx  = {}'.format(Idx))

    return ndx, Idx

# --------------------------------- #
#   to generate all symbols         #
# --------------------------------- #
def Generate_All_labels(Nk, symbols = ['I', 'X', 'Y', 'Z']):

    symList = symbols

    for i in range(1,Nk):
        print('   the {}-th qubit'.format(i))

        sym_Generated = []
        for symNow in symList:
            sym_Generated = sym_Generated + [''.join([symNow, s]) for s in symbols]
            #print(sym_Generated)   
        symList = sym_Generated
    #print(symList)
    print('  totol number of labels {}'.format(len(symList)))

    return symList


def test_Delete_Identity():
    # ------------------------------------------------- #
    #   <       testing        >                        #
    #                                                   #
    #   check if there exists identity I in label_list  #
    # ------------------------------------------------- #

    label_list.append('III')
    #label_list.append('IXY')
    label_list.remove('IZZ')
    label_list.append('III')
    #label_list.append('IXX')
    print(label_list)

    ndx, Idx = To_Find_Iden_Idx(label_list, Nk)

    label_list = Delete_Iden(label_list, Nk, m)


if __name__ == '__main__':


    StateName = 'GHZ'          #  'GHZ' | 'Had' | 'rand'
    Nr        = 1


    New_Pj_shot = [0, 0]

    if StateName == 'rand':
        Noise = [0, 0, 0.1]     # (eg) [method = 1 (Gaussian), parameters = (mu, sigma)]  
        StVer = [1, 0]          #  [version, New or Generate]   for (random) State        
    else:
        Noise = 0            #  = 0 -> the usual shot model 
        StVer = 0            #  by default = 0
        
    version     = [2, 1, Noise]

    New_Pj_shot = [1, 1]  #  [New Proj?, New shot?]  | [0,0] loading | 


    Nk = 3                      # number of qubit
    m  = 5                      # number of labels

    mea = 100                   # mea =  number of shots


    #State_choice = 1            ##  1: GHZ  |  2: Had
    #input_S, target_density_matrix, rho, Dir, Name, m_store, p_store = Init_Rho_Gen(Nk, State_choice)

    Dir, Name, params_setup, version, proj_path, meas_path, Dir0 = \
            State_Naming(Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer)
    

    # ----------------------------------------- #
    #   generate rho | target_density_matrix    #
    # ----------------------------------------- #

    if StateName == 'rand':
        target_density_matrix, rho = Gen_randM(Nk, StateName, Nr, New_Pj_shot, StVer)
    else:
        input_S, target_density_matrix, rho = Gen_Rho(Nk, StateName)

    #X = target_density_matrix
    X = rho

    check_reconstruct_error(X)


    # ------------------------------------------------------------------ #
    #           projectors  ==> generate / loading  label_list           #
    # ------------------------------------------------------------------ #
    if New_Pj_shot[0] == 1:     #  new projectors
        tp0 = time.time()

        print(' -------------     Generate new labels    ------------')
        label_list = projectors.generate_random_label_list(m, Nk)
        print(np.array(label_list))

        # --------------------------------- #
        #   generate & store projectors     #
        # --------------------------------- #
        Get_param.Do_projector(label_list, proj_path)

    elif New_Pj_shot[0] == 0:     #  loading projectors
        # ------------------------- #
        #   loading projectors      #
        # ------------------------- #
        label_list = [fname.split('.')[0] for fname in os.listdir(proj_path)]

    print('  label_list = {} \n'.format(label_list))

    # ---------------------------------------------- #
    #   can only (A) or (B),  cannot be both         #
    #       otherwise -->  labels more than m        #
    #   Generate (A) ->  shot measurement            #
    #            (B) ->  exact + noise               #
    #   Loading  (C) ->  load from noise model (B)   #
    # ---------------------------------------------- #

    if StateName != 'rand':
        # ------------------------------------------------------------------------- #
        #   (A) Loading qiskit (having labels & measurement_dict)                   #
        #      1. -->  Get_param.Do_measure (save measurement_store_path)           #
        #          --> Get_measurement_list (convert from saved measurement_dict)   #
        #      2. (parallel)   Get_param.Do_projector(Dir, labels)                  #
        #          --> to save projectors for labels                                #
        #                                                                           #
        # ------------------------------------------------------------------------- #

        #  (1)  generate new labels
        measurement_dict0, data_dict_list, backend = \
            state_measure(Dir, params_setup, label_list, input_S, target_density_matrix, rho)

        #  (2)  loading labels from earlier qiskit
        s_label, rho, measurement_dict = Load_qiskit(Dir, params_setup)

        Get_param.Do_measure(meas_path, measurement_dict)
        label_list, measurement_list = Get_measurement_list(meas_path)

    elif StateName == 'rand':
        # ------------------------------------------------------------------------- #
        #  (B)  generate new labels from projectors.py                              #
        #          i.e.  generate_random_label(n, symbols=['I', 'X', 'Y', 'Z'])     #
        #                generate_random_label_list(m, Nk)                          #
        #                                                                           #
        #      -->  Get_param.Do_projector(Dir, labels)                             #
        #                    to save projectors                                     #
        #      -->  loading projectors from projector_store_path                    #
        #      -->  get measurement_list for X, i.e.  do exact measurement for X    #      
        #               from direct trace coef for X                                #
        # ------------------------------------------------------------------------- #

        params_setup['Noise'] = Noise

        s_label, yProj, zN, yProj_Exact = \
            state_S_coef(Dir, params_setup, target_density_matrix, rho)

        # --------------------------------- #
        #   loading random matrix           #
        # --------------------------------- #
        Noise    = version[2]
        Model    = Noise[0]
        ver_meas = params_setup['version'][1]

        Fname1 = '{}_zN{}_v{}_Noise'.format(Dir, Model, ver_meas)
        Fname2 = '{}_zN{}_v{}_measurement'.format(Dir, Model, ver_meas)

        print(' Fname = {}'.format(Fname1))
 
        #with open(Fname1, 'rb') as f:
        #    yProj_Exact, zN, yProj, Noise, s_label, params_setup, target_density_matrix, rho = pickle.load(f)

        with open(Fname2, 'rb') as f:
            s_label_Ld, yProj_Ld, zN_Ld = pickle.load(f)


    # --------------------------------- #
    #   generate all label_list         #
    # --------------------------------- #
    labels = Generate_All_labels(Nk)
    
    #test_Delete_Identity()

