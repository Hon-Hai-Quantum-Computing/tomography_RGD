
import time
import sys
import warnings

#sys.path.append('./qutomo')
sys.path.append('./quQST')
warnings.filterwarnings("ignore")

import methodsRGD


from importlib import reload
reload(methodsRGD)

import pickle

import Utility

from worker_saver import worker_container


if __name__ == '__main__':
    
    #  ---------------------------------------- #
    #   some input data
    #
    #   Nk = num_qubits  |   m = num_labels
    # ----------------------------------------- #

    #exec(open('Input_param.py').read()) ## read (version = 1, In_Nk = 3, In_m = 50, In_shot = 100, In_mu = 0.2)
    #exec(open('Input_param_Loop_v1.py').read())

    # ----------------------------- #
    #   using methodsRGD            #
    # ----------------------------- #

    tw2a = time.time()

    worker = methodsRGD.BasicWorkerRGD(params_dict, input_S)
    worker.computeRGD(InitX_RGD, Ch_svd, Md_tr, Md_alp, Md_sig)

    tw2b = time.time()
    #print(' --------- RGD worker time = {}  --------'.format(tw2b - tw2a))

    # --------------------------------------- ##
    #   record the RGD result                 ##
    # --------------------------------------- ##
    ver_mea = version[1] 

    if worker.Noise == None:            #  the normal shot measurements
        Fname = '{}_shot{}_v{}-RGD_Ix{}_Tr{}_Ap{}_sg{}'.format(Dir, mea, ver_mea, InitX_RGD, Md_tr, Md_alp, Md_sig)
        print('   RGD -->  Fname = {}'.format(Fname))
    else:                               #  the noise model case
        zModel  = version[2][0]
        Fname = '{}_zN{}_v{}-RGD_Ix{}_Tr{}_Ap{}_sg{}'.format(Dir, zModel, ver_mea, InitX_RGD, Md_tr, Md_alp, Md_sig)


    # --------------------------------------- #
    save_wk_type = 1
    if save_wk_type == 0:
        Fname = '{}.dat'.format(Fname)

        with open(Fname, 'wb') as f:
            pickle.dump([params_dict, 'RGD', worker, tw2b-tw2a], f)

    elif save_wk_type == 1:
        Fname = '{}_wrapper.dat'.format(Fname)

        wc = worker_container(worker)

        with open(Fname, 'wb') as f:
            pickle.dump([params_dict, 'RGD', wc, tw2b-tw2a], f)
        
        print('  RGD: worker container saved in Fname done = {}\n'.format(Fname))




