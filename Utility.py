
import numpy as np
import os

import pickle

import sys
import warnings

warnings.filterwarnings("ignore")

#sys.path.append('./qutomo')
#import methods


sys.path.append('./quQST')

#from quQST import projectors
#from quQST import states
#from quQST import measurements
#from quQST import methodsRGD

import projectors
import states
import measurements
import methodsRGD


#sys.path.append('./qutomo')
#import methods
#import methodsMiFGD as methods

from importlib import reload        ##  for reloading
reload(projectors)

import qutip as qu
#from Utility_RGD import Sm_dxd_list, A_1
import Get_param

# -------------------------- #
import time

import multiprocessing
from itertools import repeat
from functools import reduce


def Rnd_StVer_Check(StateName, Dir0, StVer):

    # ----------------------------------------- #
    #   to generated or loading random state    #
    # ----------------------------------------- #
    if StateName == 'rand':
        
        print(' StVer = {}'.format(StVer))

        ver_Data = StVer[0]
        Dir_Data = '{}/R{}'.format(Dir0, ver_Data)

        Fname = '{}/Rho.dat'.format(Dir_Data)
        print('  is Fname = {}  exists?  {}'.format(Fname, os.path.isfile(Fname)))

        if StVer[1] == 1:               #  generate (state | DM) data

            while os.path.isfile(Fname):
                ver_Data += 1
                Dir_Data = '{}/R{}'.format(Dir0, ver_Data)
                Fname = '{}/Rho.dat'.format(Dir_Data)
                print('  Dir_Data = {}'.format(Dir_Data))            
            StVer = [ver_Data, StVer[1]]

            if not os.path.exists(Dir_Data):
                os.makedirs(Dir_Data)
            Dir0 = Dir_Data
        
            #print('  Dir0  = {}'.format(Dir0))
            #print('  StVer = {}'.format(StVer))

        elif StVer[1] == 0:             #  loading (state | DM) data

            if not os.path.isfile(Fname):
                print('  File NOT exists')
                return
            
            Dir0 = Dir_Data

    return Dir0, StVer


def State_Naming(Nk, StateName, num_labels, num_shots, Nr, version, New_Pj_shot, StVer, \
                Pj_method, mea_method, measure_method, Data_In=[]):


    Name = '{}-{}'.format(StateName, Nk)
    if StateName == 'rand':
        Name = '{}-r{}'.format(Name, Nr)


    #Dir0 = 'data/{}'.format(Name)
    Dir0 = 'calc/{}'.format(Name)
    #Dir0 = 'paper/{}'.format(Name)

    if not os.path.exists(Dir0):
        os.makedirs(Dir0)

    # ----------------------------------------- #
    #   to generated or loading random state    #
    # ----------------------------------------- #
    
    #print(' Dir0  = {}'.format(Dir0))
    #print(' StVer = {}'.format(StVer))

    Dir0, StVer = Rnd_StVer_Check(StateName, Dir0, StVer)
    #print(' Dir0  = {}'.format(Dir0))
    #print(' StVer = {}'.format(StVer))

    # ----------------------------- #
    #   record some information     #
    # ----------------------------- #

    params_setup = {'n': Nk,
                   #'Name': Name,
                   'StateName': StateName, 
                   'num_labels': num_labels,
                   'num_shots': num_shots,
                   'Nr': Nr,
                   #'version': version,
                   #'StVer': StVer,
                   #'Noise': Noise,
                   #'Dir0': Dir0,
                   #'Dir': Dir,
                   #'measurement_store_path': Dir_meas,
                   #'projector_store_path': Dir_proj,
                   'Pj_method': Pj_method,
                   'mea_method': mea_method,
                   'measure_method': measure_method
                   }

    if len(Data_In) > 0:
        Dir, params_setup, Dir_proj, Dir_meas = \
            Naming_Data_In(params_setup, Data_In, Name, version, StVer, Dir0)

        return Dir, Name, params_setup, version, Dir_proj, Dir_meas, Dir0, StVer

    # ----------------------------- #
    #   deal with projectors        #
    # ----------------------------- #
    ver_Prj = version[0]
    #Dir = '{}/{}_m{}_s{}'.format(Dir0, Name, num_labels, ver_Prj)
    Dir = '{}/{}_m{}_s{}/{}_m{}_s{}'.format(Dir0, Name, num_labels, ver_Prj, 
                                            Name, num_labels, ver_Prj)
    if Pj_method == 0:      #  original method
        Pj_Name = 'projectors'
    elif Pj_method == 1:    #  new method
        Pj_Name = 'Proj'

    Dir_proj = '{}_{}'.format(Dir, Pj_Name)
    #Dir_proj = '{}/{}_m{}_s{}_projectors'.format(Dir, Name, num_labels, ver_Prj)

    if New_Pj_shot[0] == 1:
        while os.path.exists(Dir_proj):
            ver_Prj = ver_Prj + 1
            Dir = '{}/{}_m{}_s{}/{}_m{}_s{}'.format(Dir0, Name, num_labels, ver_Prj,
                                                    Name, num_labels, ver_Prj)
            Dir_proj = '{}_{}'.format(Dir, Pj_Name)
            #Dir_proj = '{}/{}_m{}_s{}_projectors'.format(Dir, Name, num_labels,  ver_Prj)
            print(' version (of projectors) = {}'.format(ver_Prj))
            version[0] = ver_Prj

    # ----------------------------- #
    #   deal with shots or noise    #
    # ----------------------------- #

    if type(version[2]) == int:         #  the normal case (with shots)

        if mea_method == 0:        #  original method
            MeasureName = 'measurements'
        elif mea_method == 1:      #  new method
            MeasureName = 'Measure'

        ver_Shot = version[1]
        Dir_meas = '{}_shot{}_v{}_{}'.format(Dir, num_shots, ver_Shot, MeasureName)

        if New_Pj_shot[1] == 1:

            #while os.path.exists(DirQiskit):
            #    version = version + 1
            #    Dir = '{}/{}_m{}_shot{}_v{}'.format(Dir0, Name, num_labels, num_shots, version)
            #    DirQiskit = '{}_qiskit.dat'.format(Dir)
            #    print("version now = {}".format(version))
        
            while os.path.exists(Dir_meas):
                ver_Shot = ver_Shot + 1
                Dir_meas = '{}_shot{}_v{}_{}'.format(Dir, num_shots, ver_Shot, MeasureName)
                version = [ver_Prj, ver_Shot, version[2]]
                print("version now = {}".format(version))
    
            version = [ver_Prj, ver_Shot, version[2]]

    elif type(version[2]) == list:     #  the noise model case
        print(' -----  noise model  --------- ')

        ver_Mode1 = version[1]         #  version for this model
        zModel    = version[2][0]      #  specify the noise model
        Dir_meas = '{}_zN{}_v{}_measurements'.format(Dir, zModel, ver_Mode1)
        print('  Dir_meas  = {}'.format(Dir_meas))

        if New_Pj_shot[1] == 1:
            while os.path.exists(Dir_meas):
                
                ver_Mode1 = ver_Mode1 + 1
                Dir_meas  = '{}_zN{}_v{}_measurements'.format(Dir, zModel, ver_Mode1)
                version = [ver_Prj, ver_Mode1, version[2]]

                print("     version now = {}".format(version))

            version = [ver_Prj, ver_Mode1, version[2]]

    print('         version = {}'.format(version))

    Noise = version[2]

    params_setup['Name']                   = Name
    params_setup['version']                = version
    params_setup['StVer']                  = StVer
    params_setup['Noise']                  = Noise    
    params_setup['Dir0']                   = Dir0
    params_setup['Dir']                    = Dir
    params_setup['measurement_store_path'] = Dir_meas
    params_setup['projector_store_path']   = Dir_proj


    #params_setup = {'n': Nk,
    #               'Name': Name,
    #               'StateName': StateName, 
    #               'num_labels': num_labels,
    #               'num_shots': num_shots,
    #               'Nr': Nr,
    #               'version': version,
    #               'StVer': StVer,
    #               'Noise': Noise,
    #               'Dir0': Dir0,
    #               'Dir': Dir,
    #               'measurement_store_path': Dir_meas,
    #               'projector_store_path': Dir_proj,
    #               'Pj_method': Pj_method,
    #               'mea_method': mea_method,
    #               'measure_method': measure_method
    #               }


    return Dir, Name, params_setup, version, Dir_proj, Dir_meas, Dir0, StVer




def Naming_Data_In(params_setup, Data_In, Name, version, StVer, Dir0):
    Dir      = Data_In[0]
    pjFile   = Data_In[1]

    Dir_proj = '{}/{}'.format(Dir, pjFile)

    v_meas = version[1]
    Noise = version[2]
    zModel = Noise[0]

    Dir_meas = '{}/zN{}_v{}_measurements'.format(Dir, zModel, v_meas)


    params_setup['Name']                   = Name
    params_setup['version']                = version
    params_setup['StVer']                  = StVer
    params_setup['Noise']                  = Noise    
    params_setup['Dir0']                   = Dir0
    params_setup['Dir']                    = Dir
    params_setup['measurement_store_path'] = Dir_meas
    params_setup['projector_store_path']   = Dir_proj

    return Dir, params_setup, Dir_proj, Dir_meas


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



#def Init_Rho_Gen(Nk, State_Choice, version):
def Gen_Rho_v0(Nk, StateName):
    """  (input)
            version to record which version we are recording
    """

    if StateName=='GHZ':             #  StateName = 'GHZ'    
        target = states.GHZState(Nk)        # Nk = num_qubits

    elif StateName=='Had':           #  StateName = 'Had'          
        target = states.HadamardState(Nk)   #  states from states.py

    # ----------------- #
    #   use qu          #
    # ----------------- #
    if StateName=='GHZ':             #  StateName = 'GHZ'    
        input_ket = qu.ghz_state(Nk)        #  should = input_S   ( = Nk-qubit  GHZ state )

    elif StateName=='Had':           #  StateName = 'Had'                  
        input_ket = (qu.basis(2, 0) + qu.basis(2, 1)).unit()  # .unit()  --> to normalized
        for ii in range(Nk - 1):
            input_ket = qu.tensor(input_ket, (qu.basis(2, 0) + qu.basis(2, 1)).unit())

    rho = qu.ket2dm(input_ket)          #  method 2: using qutip

    #target.create_circuit()
    target_density_matrix = target.get_state_matrix()

    return target, target_density_matrix, rho


def Gen_randM(Nk, StateName, Nr, StVer, Dir0):
    """
        Nr = 1                  #  pure state -->  rank = 1    
    """

    Fname = '{}/Rho.dat'.format(Dir0)

    if StateName=='rand':         #  randomized matrix

        if StVer[1] == 1:           #  to generate data
            rho = qu.rand_dm_ginibre(2**Nk, dims=[[2]*Nk , [2]*Nk], rank=Nr)
            target_density_matrix = rho.full()

            with open(Fname, 'wb') as f:
                pickle.dump([rho, target_density_matrix], f)

        elif StVer[1] == 0:         #  loading data

            print(' is file exists?  {}'.format(os.path.isfile(Fname)))
            
            with open(Fname, 'rb') as f:
                rho, target_density_matrix = pickle.load(f)

    return target_density_matrix, rho


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
    #num_qubits = params_setup['n']
    #num_labels = params_setup['num_labels']

    #print(' -------------     Generate new labels    ------------')
    #labels = projectors.generate_random_label_list(num_labels, num_qubits)
    #print(np.array(labels))

    #num_shots  = params_setup['num_shots']
    #version    = params_setup['version']
    #ver_Shot   = version[1]

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


#def state_measure(Dir, params_setup, state, target_density_matrix, rho):
def state_S_coef(Dir, params_setup, target_density_matrix, rho):
    """  (input)
            version to record which version we are recording
         (eg)
                num_qubits = 3
                num_labels = 50
                num_shots = 100
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
    if method == 1:
        # ----------------------------------------- #
        #   constructing Pauli matrix from qutip    #
        # ----------------------------------------- #
        S = Sm_dxd_list(labels)

        measurement_scale = A_1(S, Nk, m, qu.Qobj(rho, dims=[[2]*Nk, [2]*Nk]))

        coefS_list = measurement_scale * coefI

        while not os.path.exists(Fname):
            with open(Fname, 'wb') as f:
                pickle.dump([coefS_list, \
                        labels, params_setup, target_density_matrix, rho], f)

        #return np.asarray(labels), coefS_list
        yProj_Exact = coefS_list

    elif method == 2:
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
    Model = Noise[0]

    if Model == 0:               # no noise
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

    Fname1 = '{}/zN{}_v{}_Noise'.format(meas_path, Model, ver_meas)
    Fname2 = '{}/zN{}_v{}_measurements'.format(meas_path, Model, ver_meas)


    if not os.path.exists(meas_path):
        os.mkdir(meas_path)


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

def qiskit_measure(Dir, Name, num_qubits, num_labels, num_shots, state, target_density_matrix, rho, version):
    """  (input)
            version to record which version we are recording
         (eg)
                num_qubits = 3
                num_labels = 50
                num_shots = 100

    """
 
    labels = projectors.generate_random_label_list(num_labels, num_qubits)
    #np.array(labels)

    backend = 'qasm_simulator'                  ## Circuits ##


    state.create_circuit()
    data_dict_list = state.execute_measurement_circuits(labels,
                                                        backend=backend,
                                                        num_shots=num_shots)

    ## Measurements ##
    ##
    measurement_dict = {item['label']: item['count_dict'] for item in data_dict_list}

    params_qiskit = {'n': num_qubits,
                   'num_labels': num_labels,
                   'backend': backend,
                   'num_shots': num_shots}

    #Fname = 'data/GHZ_Nk{}_m{}_v{}.dat'.format(Nk, m, version)
    #Fname = 'data/{}_m{}_shot{}_qiskit_v{}.dat'.format(Name, num_labels, num_shots, version)
    Fname = '{}/{}_m{}_shot{}_qiskit_v{}.dat'.format(Dir, Name, num_labels, num_shots, version)

    while os.path.exists(Fname):
        version = version + 1
        Fname = '{}/{}_m{}_shot{}_qiskit_v{}.dat'.format(Dir, Name, num_labels, num_shots, version)
        print("version now = {}".format(version))

    with open(Fname, 'wb') as f:
        #pickle.dump([num_qubits, num_labels, num_shots, state, measurement_dict, data_dict_list], f)
        pickle.dump([num_qubits, num_labels, num_shots, state, measurement_dict, data_dict_list, \
                     labels, params_qiskit, target_density_matrix, rho], f)

    return np.array(labels), measurement_dict, data_dict_list, params_qiskit, version


#def Init_Rho_Gen(Nk, State_Choice, version):
def Init_Rho_Gen(Nk, State_Choice):
    """  (input)
            version to record which version we are recording
    """
    ##  ------  the target: ground truth density matrix  ---------  ##
    ##   ground truth density matrix  from qiskit   ##
    ##

    # -------------- initial preparation of density matrix  ------------- #

    Nr = 1                  #  pure state -->  rank = 1

    # rho = gen_rho(d,Nk, Nr, Prob); rho = qu.Qobj(rho, dims=[[2]*Nk, [2]*Nk])   #  method 1: my own method Gen_DM.py
    # rho = qu.rand_dm_ginibre(d, dims=[[2]*Nk , [2]*Nk], rank=Nr)

    if State_Choice==1:
        StateName = 'GHZ'          #  state Name

        #rho = qu.ket2dm(qu.ghz_state(Nk))      #  method 2: using qutip
        input_ket = qu.ghz_state(Nk)        #  should = input_S   ( = Nk-qubit  GHZ state )

        target = states.GHZState(Nk)        # Nk = num_qubits

    elif State_Choice==2:

        StateName = 'Had'          #  state Name

        # ---  to construct Hadamard state ----- #
        #input_ket = qu.Qobj([[1]]*(2**Nk))  / (np.sqrt(2)**Nk)    #  Hadamdard state  ( 2**Nk dim -> no good)

        input_ket = (qu.basis(2, 0) + qu.basis(2, 1)).unit()  # .unit()  --> to normalized
        for ii in range(Nk - 1):
            input_ket = qu.tensor(input_ket, (qu.basis(2, 0) + qu.basis(2, 1)).unit())

        target = states.HadamardState(Nk)   #  states from states.py

    rho = qu.ket2dm(input_ket)          #  method 2: using qutip

    target.create_circuit()
    target_density_matrix = target.get_state_matrix()

    # rho = target.get_state_matrix()
    #rho = qu.Qobj(target.get_state_matrix(), dims=[[2] * Nk, [2] * Nk])

    # -----------  Naming for recording  -------------------- #

    #m_store = 'data/measurements_GHZ{}_MiFGD-v{}'.format(Nk, version)
    #p_store = 'data/projectors_GHZ{}_MiFGD-v{}'.format(Nk, version)
    #m_store = 'data/measurements_{}{}_MiFGD'.format(StateName, Nk)
    #p_store = 'data/projectors_{}{}_MiFGD'.format(StateName, Nk)

    Name = '{}-{}'.format(StateName, Nk)
    Dir = 'data/{}'.format(Name)

    m_store = '{}/measurements_{}{}'.format(Dir, StateName, Nk)
    p_store = '{}/projectors_{}{}'.format(Dir, StateName, Nk)

    if not os.path.exists(Dir):
        os.makedirs(Dir)

    return target, target_density_matrix, rho, Dir, Name, m_store, p_store


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


def Get_Pauli_coef(Dir, label_list, X, m, Nk):

    # ------------------------------------------------------------- #
    #   get measurement_list from direct trace of Pauli matrix  S   #
    # ------------------------------------------------------------- #

    S = Sm_dxd_list(label_list)

    mS_scale = A_1(S, Nk, m, qu.Qobj(X, dims=[[2]*Nk, [2]*Nk]))
    #mS2 = methodsRGD.A_1(S, Nk, m, qu.Qobj(X, dims=[[2]*Nk, [2]*Nk]))

    # --------------------------------- #
    #   loading the projectors          #
    # --------------------------------- #
    projector_store_path = '{}_projectors'.format(Dir)

    start = 0
    end = m
    projector_dict = projectors.ProjectorStore.load(projector_store_path, label_list[start:end])
    #projector_dict = projectors.ProjectorStore.load(projector_store_path, label_list)

    projector_list = [projector_dict[label] for label in label_list]

    yAmea = methodsRGD.Amea(projector_list, X, m, 2**Nk)


    coefI = np.sqrt(m / 2**Nk)

    mS    = mS_scale * coefI
    yProj = yAmea    * coefI

    crite = 1e-10
    if np.max(np.abs(mS - yProj)) > crite:
        print(' ERROR btw  mS  &  yProj')
        return

    return mS, yProj


def Reorder_Pj_Get_coef(Dir, label_list, proj_path, X, m, Nk):
    # --------------------------------------------------------- #
    #   generate & store projectors  from input label_list      #
    # --------------------------------------------------------- #

    label_list_P = [fname.split('.')[0] for fname in os.listdir(proj_path)]

    #print(label_list)
    #print(label_list_P)
    if set(label_list) != set(label_list_P):
        print(' ERROR: inconsistent labels for loading projectors')
        return 


    mS, yProj = Get_Pauli_coef(Dir, label_list_P, X, m, Nk)

    return mS, yProj, label_list_P

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
        Noise = [0, 0, 0.1]         # (eg) [method = 1 (Gaussian), parameters = (mu, sigma)]  
        StVer = [1, 0]       #  [version, New or Generate]   for (random) State        
    else:
        Noise = 0                   #  = 0 -> the usual shot model 
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

    # ------------------------------------------------------ #
    #   [comparison | verification]                          #
    #   loading Proj --> get trace directly for comparison   #
    # ------------------------------------------------------ #
    mS_E1, yProj_E1, labels_P = Reorder_Pj_Get_coef(Dir, s_label, proj_path, rho, m, Nk)


    # --------------------------------- #
    #   generate all label_list         #
    # --------------------------------- #
    labels = Generate_All_labels(Nk)
    
    #test_Delete_Identity()

