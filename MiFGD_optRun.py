
import time
import sys
import warnings

#sys.path.append('./qutomo')
sys.path.append('./quQST')
warnings.filterwarnings("ignore")

import methodsMiFGD

from importlib import reload    ## for reloading
reload(methodsMiFGD)

import pickle

from worker_saver import worker_container


#def Run_MiFGD(params_dict, mu, InitX_MiFGD):
def Run_MiFGD(params_dict, Rpm):

    InitX_MiFGD, mu, eta, Option = Rpm

    print('\n')
    print(' ****************************************************** ')
    print(' **      start calculating MiFGD with mu = {}     **'.format(mu))
    print(' ****************************************************** \n')

    # ----------------------------- #
    #   using methodsRGD            #
    # ----------------------------- #
    params_dict['eta']    = eta       #  only for MiFGD   (eg)  eta = 0.01
    params_dict['mu']     = mu        #  only for MiFGD
    params_dict['Option'] = Option    #  only for MiFGD

    tw2a = time.time()

    #worker = methodsMiFGD.BasicWorker(params_dict, input_S)
    worker = methodsMiFGD.BasicWorker(params_dict)

    worker.compute(InitX_MiFGD)
    
    tw2b = time.time()
    #print(' *****       mu ={}    -->   elapsed time = {}    ***** '.format(mu, tw2b - tw2a))

    RunTime = tw2b - tw2a

    # -------------------------------------------------------- ##
    #   read out needed parameters for specifying File name    ##
    # -------------------------------------------------------- ##
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
    if worker.Noise == 'shot':           #  the normal shot measurements
        Fname = '{}shot{}_v{}-MiFGD-mu{:.3}'.format(Dir, mea, ver_mea, mu)

    elif worker.Noise == 'Exact':       #  the Exact coef -> no noise
        #zModel  = version[2]           #  should be 0
        #Fname = '{}_zN{}_v{}-MiFGD-mu{:.3}'.format(Dir, zModel, ver_mea, mu)
        #Fname = '{}_zN0_v{}-MiFGD-mu{:.3}'.format(Dir2m, ver_mea, mu)
        Fname = '{}zExact_v{}-MiFGD-mu{:.3}'.format(Dir, ver_mea, mu)
    else:
        Fname = 'NotExist'
        print(' +++++   Noise model does NOT exist !  ++++')

    if worker.InitX != 0:
        Fname = '{}-Ix{}'.format(Fname, worker.InitX)

    if worker.Option != 0:
        Fname = '{}-eT{}'.format(Fname, worker.Option)

    # ---------------------------------------------------- #
    #   start to save calculated worker/wc into File       #
    # ---------------------------------------------------- #
    save_wk_type = 1
    if save_wk_type == 0:
        Fname = '{}.dat'.format(Fname)

        with open(Fname, 'wb') as f:
            pickle.dump([params_dict, 'MiFGD', worker, tw2b-tw2a], f)

        print('  MiFGD worker saved in Fname done = {}\n'.format(Fname))
        return Fname, worker, RunTime

    elif save_wk_type == 1:
        Fname = '{}_wrapper.dat'.format(Fname)

        wc  = worker_container(worker)
        #del worker

        with open(Fname, 'wb') as f:
            pickle.dump([params_dict, 'MiFGD', wc, tw2b-tw2a], f)

        print('  MiFGD: worker container saved in Fname done = {}\n'.format(Fname))

    return Fname, wc, RunTime

if __name__ == '__main__':


    Frec_MiFGD, wc, RunTime = Run_MiFGD(params_dict, mu, InitX_MiFGD)



