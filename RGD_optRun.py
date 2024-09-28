
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

from worker_saver import worker_container


def Run_RGD(params_dict, Rpm):
    #  ---------------------------------------- #
    #   some variables
    #
    #   Nk = num_qubits  |   m = num_labels
    # ----------------------------------------- #
    print('\n')
    print(' ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ')
    print(' ++           start calculating  RGD                 ++ ')
    print(' ++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n')

    # ----------------------------- #
    #   using methodsRGD            #
    # ----------------------------- #

    InitX_RGD, Md_tr, Md_alp, Md_sig, Ch_svd = Rpm

    tw2a = time.time()

    #worker = methodsRGD.BasicWorkerRGD(params_dict, input_S)
    worker = methodsRGD.BasicWorkerRGD(params_dict)

    worker.computeRGD(InitX_RGD, Ch_svd, Md_tr, Md_alp, Md_sig)

    tw2b = time.time()
    #print(' --------- RGD worker time = {}  --------'.format(tw2b - tw2a))

    RunTime = tw2b - tw2a

    # ------------------------------------------------------ ##
    #   read out needed parameters for specifying FileName   ##
    # ------------------------------------------------------ ##
    version = params_dict['version']
    ver_mea = version[1] 

    if params_dict['Noise'] == 'shot':
        mea = params_dict['num_shots']

    if params_dict['StateName'] == 'KapRnd':
        Dir = '{}/'.format(params_dict['measurement_store_path'])
    else:
        Dir2m = params_dict['Dir2m']
        Dir   = '{}_'.format(Dir2m)

    # ------------------------------------------------ ##
    #   specify the FileName to record the RGD result  ##
    # ------------------------------------------------ ##

    #if worker.Noise == None:            #  the normal shot measurements
    if worker.Noise == 'shot':            #  the normal shot measurements
        Fname = '{}shot{}_v{}-RGD_Ix{}_Tr{}_Ap{}_sg{}'.format(Dir, mea, ver_mea, InitX_RGD, Md_tr, Md_alp, Md_sig)
        print('   RGD -->  Fname = {}'.format(Fname))

    elif worker.Noise == 'Exact':       #  the Exact coef -> no noise
        zNoise  = version[2]            #  should be 0
        #Fname = '{}_zN{}_v{}-RGD_Ix{}_Tr{}_Ap{}_sg{}'.format(Dir, zModel, ver_mea, InitX_RGD, Md_tr, Md_alp, Md_sig)
        #Fname = '{}_zN0_v{}-RGD_Ix{}_Tr{}_Ap{}_sg{}'.format(Dir, ver_mea, InitX_RGD, Md_tr, Md_alp, Md_sig)
        Fname = '{}zExact_v{}-RGD_Ix{}_Tr{}_Ap{}_sg{}'.format(Dir, ver_mea, InitX_RGD, Md_tr, Md_alp, Md_sig)
    else:
        Fname = 'NotExist'
        print('Noise model does NOT exist!')

    # ----------------------------------------------- #
    #   start to save calculated worker/wc into File  #
    # ----------------------------------------------- #
    save_wk_type = 1
    if save_wk_type == 0:
        Fname = '{}.dat'.format(Fname)

        with open(Fname, 'wb') as f:
            pickle.dump([params_dict, 'RGD', worker, tw2b-tw2a], f)

    elif save_wk_type == 1:
        Fname = '{}_wrapper.dat'.format(Fname)

        wc = worker_container(worker)
        #del worker

        with open(Fname, 'wb') as f:
            pickle.dump([params_dict, 'RGD', wc, tw2b-tw2a], f)
        
        print('  RGD: worker container saved in Fname done = {}\n'.format(Fname))


    return Fname, wc, RunTime


if __name__ == '__main__':

    Frec_RGD, wc, RunTime = Run_RGD(params_dict, Rpm)




