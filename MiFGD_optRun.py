
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

if __name__ == '__main__':

    tw2a = time.time()

    worker = methodsMiFGD.BasicWorker(params_dict, input_S)
    worker.compute(InitX_MiFGD)

    tw2b = time.time()
    #print(' *****       mu ={}    -->   elapsed time = {}    ***** '.format(mu, tw2b - tw2a))

    # ----------------------------------------- ##
    #   record the MiFGD result                 ##
    # ----------------------------------------- ##
    ver_mea = version[1] 

    if InitX_MiFGD >= 0:

        if worker.Noise == None:            #  the normal shot measurements
            Fname = '{}_shot{}_v{}-MiFGD-mu{:.3}'.format(Dir, mea, ver_mea, mu)

        else:                               #  the noise model case
            zModel  = version[2][0]
            Fname = '{}_zN{}_v{}-MiFGD-mu{:.3}'.format(Dir, zModel, ver_mea, mu)


        save_wk_type = 1
        if save_wk_type == 0:
            Fname = '{}.dat'.format(Fname)

            with open(Fname, 'wb') as f:
                pickle.dump([params_dict, 'MiFGD', worker, tw2b-tw2a], f)

            print('  MiFGD worker saved in Fname done = {}\n'.format(Fname))

        elif save_wk_type == 1:
            Fname = '{}_wrapper.dat'.format(Fname)

            wc  = worker_container(worker)

            with open(Fname, 'wb') as f:
                pickle.dump([params_dict, 'MiFGD', wc, tw2b-tw2a], f)

            print('  MiFGD: worker container saved in Fname done = {}\n'.format(Fname))




