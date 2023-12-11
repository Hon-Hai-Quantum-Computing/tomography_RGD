#%%
import matplotlib.pyplot as plt

#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import numpy as np

from Utility import State_Naming
import pickle
import os

import sys



def Get_wk_RGD_MiFGD(InitX, Md_tr, Md_alp, Md_sig, mu_List, Data_pm):

    wk_list = []
    dt_list = []
    pm_list = []
    label_list = []
    # ------------------------------------- #
    #   common setting for measurement      #
    # ------------------------------------- #

    Setup = Data_Info(Data_pm)

    # -----------------     RGD     --------------------- #
    pm_RGD = 'RGD', InitX, Md_tr, Md_alp, Md_sig
    wk_list, dt_list, pm_list = Get_wk_list(Setup, version, wk_list, dt_list, pm_list, pm_RGD)

    label_list.append('RGD')

    # -----------------     MiFGD   ---------------------- #
    for mu in mu_List:
        pm_MiF = 'MiFGD', mu
        wk_list, dt_list, pm_list = Get_wk_list(Setup, version, wk_list, dt_list, pm_list, pm_MiF)

        label     = wk_list[-1].mu
        label_list.append(label)


    print(' ------------------------------------------------------------- ')
    return wk_list, dt_list, label_list



def Get_wk_diff_m_single_method(m_List, Data_pm, pm_method):

    wk_list = []
    dt_list = []
    pm_list = []
    label_list = []

    Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method = Data_pm

    for m in m_List:
        Data_pm = [Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method]

        Setup = Data_Info(Data_pm)

        # -----------------     RGD / MiFGD    --------------------- #
        wk_list, dt_list, pm_list = Get_wk_list(Setup, version, wk_list, dt_list, pm_list, pm_method)

        label_list.append(' {}, m = {},  shot = {}'.format(pm_method[0], m, mea))

    return wk_list, dt_list, label_list


def Get_wk_diff_m(m_List, InitX, Md_tr, Md_alp, Md_sig, mu_List, Data_pm):

    pm_RGD = 'RGD', InitX, Md_tr, Md_alp, Md_sig
    wk_listR, dt_listR, label_listR = Get_wk_diff_m_single_method(m_List, Data_pm, pm_RGD)
    #plt_TargetErr_RunT(wk_listR, dt_listR, label_listR, tit)

    for mu in mu_List:
        pm_MiFGD = 'MiFGD', mu
        wk_listMi, dt_listMi, label_listMi = Get_wk_diff_m_single_method(m_List, Data_pm, pm_MiFGD)
        #plt_TargetErr_RunT(wk_listMi, dt_listMi, label_listMi, tit)

    if len(mu_List) > 0:
        wk_list    = wk_listR    + wk_listMi
        dt_list    = dt_listR    + dt_listMi
        label_list = label_listR + label_listMi
    else:
        wk_list    = wk_listR   
        dt_list    = dt_listR    
        label_list = label_listR 

    return wk_list, dt_list, label_list


def Get_wk_diff_shots(mea_List, InitX, Md_tr, Md_alp, Md_sig, mu_List, Data_pm):

    Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method = Data_pm

    wk_list = []
    dt_list = []
    pm_list = []
    label_list = []

    for mea in mea_List:
        Data_pm = [Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method]

        Setup = Data_Info(Data_pm)


        # -----------------     RGD     --------------------- #
        label = ' RGD shot = {}'.format(mea)
        print('  label  --> {}'.format(label))
        label_list.append(label)

        pm_RGD = 'RGD', InitX, Md_tr, Md_alp, Md_sig
        wk_list, dt_list, pm_list = Get_wk_list(Setup, version, wk_list, dt_list, pm_list, pm_RGD)


    for mea in mea_List:
        Data_pm = [Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method]

        Setup = Data_Info(Data_pm)

        # -----------------     MiFGD   ---------------------- #
        for mu in mu_List:
            label = '  MiFGD mu = {}, shot = {}'.format(mu, mea)
            print('  label  --> {}'.format(label))
            label_list.append(label)

            pm_MiF = 'MiFGD', mu
            wk_list, dt_list, pm_list = Get_wk_list(Setup, version, wk_list, dt_list, pm_list, pm_MiF)

    return wk_list, dt_list, label_list



def Get_Wk_diff_Nk_m(Data_pm, pm_List, InitX, Md_tr, Md_alp, Md_sig, mu_List):

    Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method = Data_pm

    wk_list = []
    dt_list = []
    pm_list = []
    label_list = []

    for Nk, m in pm_List:
        print(' Nk = {}, m = {}'.format(Nk, m))

        Data_pm = [Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method]
        Setup = Data_Info(Data_pm)

        # -----------------     RGD     --------------------- #
        pm_RGD = 'RGD', InitX, Md_tr, Md_alp, Md_sig
        wk_list, dt_list, pm_list = Get_wk_list(Setup, version, wk_list, dt_list, pm_list, pm_RGD)

        label = ' RGD: {}-{}, m = {}'.format(StateName, Nk, m)
        label_list.append(label)

    for Nk, m in pm_List:
        print(' Nk = {}, m = {}'.format(Nk, m))

        Data_pm = [Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method]
        Setup = Data_Info(Data_pm)

        # -----------------     MiFGD   ---------------------- #
        for mu in mu_List:
            pm_MiF = 'MiFGD', mu
            wk_list, dt_list, pm_list = Get_wk_list(Setup, version, wk_list, dt_list, pm_list, pm_MiF)

            label   = 'MiFGD mu = {}: {}-{}, m = {}'.format(wk_list[-1].mu, StateName, Nk, m)
            label_list.append(label)

    return wk_list, dt_list, label_list


def Get_wk_MiFGD_diff_param_then_mu(parm_List, note, Data_pm, mu_List, Option):

    # ------------------------------ #
    #           for MiFGD            #
    # ------------------------------ #

    print(' *************************************** ')
    print(' ***      Now working for MiFGD      *** ')
    print(' *************************************** ')

    wk_list = []
    dt_list = []
    pm_list = []
    label_list = []

    for elem in parm_List:
        print('\n  -------------  elem = {}   ------------'.format(elem))

        Data_pm, label_0 = Option_update_params(Data_pm, elem, Option)
        Setup = Data_Info(Data_pm)
        print('      Data_pm = {}'.format(Data_pm))

        # -----------------     MiFGD   ---------------------- #
        for mu in mu_List:
            pm_MiF = 'MiFGD', mu
            wk_list, dt_list, pm_list = Get_wk_list(Setup, version, wk_list, dt_list, pm_list, pm_MiF)

            if label_0 == None:
                label = 'MiFGD mu = {}'.format(mu)
            else:
                label  = 'MiFGD mu = {}, {}'.format(mu, label_0)

            label_list.append(label)

    return wk_list, dt_list, pm_list, label_list







# ------------------------------------- #
#   the above old method                #
# ------------------------------------- #

def Compare_RunT_work(wk_list, label_list, Case_list, mu_List):
        print('\n\n')
        print(' ----------------------------------------- ')
        print(' **   check the run time & convergence   **')
        print(' ----------------------------------------- ')

        List_T_RGD = [[], []]               #  Run Time for RGD   [InitX = 0, InitX = 1]
        List_T_MiF = [[] for ii in range(len(mu_List))]     #  each elem [] for each mu
    
        List_runT  = []
        List_Err   = []
        List_Iter  = []
        List_conv  = []      #  converged or not
        List_meth  = []      #  RGD --> InitX = 0 or 1   |  MiFGD --> mu
        List_note  = []      #  elem = [RGD, 0], [RGD, 1], [MiFGD, mu]

        # --------------------------------------------- #
        #       find the runT for each method case      #
        # --------------------------------------------- #
        for worker, label in zip(wk_list, label_list):
            method        = worker.method
            converged     = worker.converged
            iteration     = worker.iteration
            Final_rho_Err = worker.Target_Err_Xk[-1]
            step_Time     = worker.step_Time

            Time_iter     = step_Time.values()
            T_run         = sum(Time_iter)

            if method == 'RGD':
                InitX = worker.InitX
                note  = 'IniX = {}'.format(InitX)
                List_note.append(['RGD', InitX])

                List_meth.append(InitX)
                List_T_RGD[InitX].append(T_run)

            elif method == 'MiFGD':
                mu = worker.mu
                note  = 'mu = {:.2}'.format(mu)
                List_note.append(['MiFGD', mu])

                List_meth.append(mu)
                List_T_MiF[mu_List.index(mu)].append(T_run)

            print(' method = {}, iter = {}, RunT = {}, converged = {}, Err = {}, {}'.format( \
                        method, iteration, T_run, converged, Final_rho_Err, note))

            List_runT.append(T_run)
            List_Err.append(Final_rho_Err)
            List_Iter.append(iteration)
            List_conv.append(converged)

        # ------------------------------------------ #
        #       find index for each method case      #
        # ------------------------------------------ #
        matches_R0  = [ii  for ii, elem in enumerate(Case_list) if elem[0]==['RGD', 0]]
        matches_R1  = [ii  for ii, elem in enumerate(Case_list) if elem[0]==['RGD', 1]]
        matches_MiF = []
        for mu in mu_List:
            matches_mu = [ii  for ii, elem in enumerate(Case_list) \
                            if elem[0]==['MiFGD', 0] and elem[1][0] == mu]
            matches_MiF.append(matches_mu)

        # ------------------------------------------------- #
        #       find the final target err for each case     #
        # ------------------------------------------------- #
        Err_R0  = np.array(List_Err)[matches_R0]
        Err_R1  = np.array(List_Err)[matches_R1]
        Err_MiF = []
        for mu, matches_mu in zip(mu_List, matches_MiF):
            Err_mu  = np.array(List_Err)[matches_mu]
            Err_MiF.append(Err_mu)

        # ------------------------------------------------- #
        #       find the converge | Iter for each case      #
        # ------------------------------------------------- #
        ALL_Iter_RGD = [[],[]]
        ALL_conv_RGD = [[],[]]
        ALL_Iter_MiF = [[] for mu in mu_List]
        ALL_conv_MiF = [[] for mu in mu_List]

        #ALL_Iter_RGD[0] = np.array(List_Iter)[matches_R0]
        ALL_Iter_RGD[0] = [List_Iter[mm] for mm in matches_R0]
        ALL_Iter_RGD[1] = [List_Iter[mm] for mm in matches_R1]

        ALL_conv_RGD[0] = [List_conv[mm] for mm in matches_R0]
        ALL_conv_RGD[1] = [List_conv[mm] for mm in matches_R1]

        for kk, mu in enumerate(mu_List):
            #ALL_Iter_MiF[kk] = np.array(List_Iter)[matches_MiF[kk]]
            #ALL_conv_MiF[kk] = np.array(List_conv)[matches_MiF[kk]]

            ALL_Iter_MiF[kk] = [List_Iter[ii] for ii in matches_MiF[kk]]
            ALL_conv_MiF[kk] = [List_conv[ii] for ii in matches_MiF[kk]]

        ALL_Iter = [ALL_Iter_RGD, ALL_Iter_MiF]
        ALL_conv = [ALL_conv_RGD, ALL_conv_MiF]

        #return List_runT, List_T_RGD, List_T_MiF, Err_R0, Err_R1, Err_MiF, List_conv, List_Iter, List_note
        return List_runT, List_T_RGD, List_T_MiF, Err_R0, Err_R1, Err_MiF, ALL_conv, ALL_Iter, List_note


def Plt_scaling_compare(parm_List, List_T_RGD, List_T_MiF, Err_R0, Err_R1, Err_MiF, \
                        ALL_Iter, ALL_conv, axNow=None):
    # ------------------------------------------ #
    #   plot the scaling of RunT versuse cases   #
    # ------------------------------------------ #

    New_Fig = 0
    if axNow == None:
        #plt.figure()
        #plt.scatter(List_Err, List_runT)
        figNew, axNew = plt.subplots(1,3, figsize=(8,6))

        New_Fig = 1
        axNow   = axNew[0]
    print(' Create New figure?  New_Fig = {}'.format(New_Fig))

    #axNow.plot(parm_List, List_T_RGD[0], 'o-', label='RGD, InitX = 0')
    #axNow.plot(parm_List, List_T_RGD[1], '^-', label='RGD, InitX = 1')

    axNow.plot(parm_List, List_T_RGD[1], '-+', label=r'RGD, $X_0 = \mathcal{H}_r(\mathcal{A}^\dagger(y))$')
    axNow.plot(parm_List, List_T_RGD[0], '-.^', label='RGD, random $X_0$')


    #mkln_list = ['+-.', 'x--', '>:', '-']

    lnmk_list = ['--o','--x', ':>','-<',2, 3]
    #mk_list = ['+','^','o','x', '>','<',2, 3]
    #ln_list = ['-', '-.', '--', '--', ':', '-']
    
    for ii, mu in enumerate(mu_List):
        axNow.plot(parm_List, List_T_MiF[ii], lnmk_list[ii], \
                    label=r'MiFGD, $\mu$ = {}'.format(mu))
        
    axNow.set_xlim([min(parm_List)*0.9, max(parm_List)*1.1])
    axNow.set_xticks(parm_List)

    #axNow.legend()

    # --------------------------------------------- #
    #   plot the final target error for each case   #
    # --------------------------------------------- #

    if New_Fig == 1:
        #plt.figure()   
        axNow = axNew[1]

        axNow.plot(parm_List, Err_R0, marker='o', label='RGD, InitX = 0')
        axNow.plot(parm_List, Err_R1, marker='^', label='RGD, InitX = 1')

        mk_list = ['+','x', '>','<',2, 3]
        ln_list = ['-.', '--', ':', '-']
        for ii, mu in enumerate(mu_List):
            axNow.plot(parm_List, Err_MiF[ii], marker=mk_list[ii], linestyle=ln_list[ii] , \
                        label='MiFGD, mu = {}'.format(mu))

        axNow.set_xlim([min(parm_List)*0.9, max(parm_List)*1.1])
        axNow.set_xticks(parm_List)
        axNow.legend()

        # ------------------------------- #
        axNow = axNew[2]
        axNow.plot(parm_List, ALL_Iter[0][0], marker='o', label='RGD, InitX = 0')
        axNow.plot(parm_List, ALL_Iter[0][1], marker='^', label='RGD, InitX = 1')

        for ii, mu in enumerate(mu_List):
            axNow.plot(parm_List, ALL_Iter[1][ii], marker=mk_list[ii], linestyle=ln_list[ii] , \
                        label='MiFGD, mu = {}'.format(mu))

        plt.xlim([min(parm_List)*0.9, max(parm_List)*1.1])
        plt.xticks(parm_List)
        plt.show()

def calc_RunT(wk_list, dt_list, label_list):
    T_sum = []
    T_rec = []
    Target_Err_wk = []
    RunT_wk       = []

    Final_Err_wk  = []           #   the final Error
    Converged     = []

    for worker, dt, label in zip(wk_list, dt_list, label_list):
        print('   ******** method = {}  -->  {}'.format(worker.method, label))
        Target_Err_Xk = worker.Target_Err_Xk
        step_Time     = worker.step_Time

        print('     Target_Err_Xk[:10] = {}'.format(Target_Err_Xk[:10]))
        #print('     step_Time     = {}'.format(step_Time))
        print('       dt               = {}'.format(dt))

        RunT = []
        Ttot = 0
        for ti in range(-1, len(step_Time)-1):
            Each_Step = step_Time[ti]
            Ttot += Each_Step
            RunT.append(Ttot)

            #Err  = Target_Err_Xk[ti+1]
            #print('            ti = {},  Time = {}, Ttot ={}  -> Err = {}'.format(ti, Each_Step, Ttot, Err))    
        print('      RunT[:10]         = {}'.format(RunT[:10]))
        print(' ----------------------------------------------------------------------')

        T_sum.append(Ttot)
        T_rec.append(dt)

        Final_Err_wk.append(Target_Err_Xk[-1])
        Converged.append(worker.converged)

        RunT_wk.append(RunT)
        Target_Err_wk.append(Target_Err_Xk)


    return T_sum, T_rec, Target_Err_wk, RunT_wk, Final_Err_wk, Converged

def call_inset(axNow, RunT_wk, Target_Err_wk, label_list, xRng, yRng, yTick, xTick=[]):
    mk_list = ['+','^','o','x', '>','<',2, 3]
    ln_list = ['-', '-.', '--', '--', ':', '-']

    method = 0
    if method == 0:
        #ax2 = inset_axes(axNow, width="30%", height="50%", loc=4, borderpad=4)

        ax2 = inset_axes(axNow, width="50%", height="75%",
                #bbox_to_anchor=(.5, .24, .6, .45),             # rand, [r]
                bbox_to_anchor=(.55, .18, .7, .6),             # rand, [Nk, m]
                #bbox_to_anchor=(.53, .3, .72, .5),           # Had, [Nk, m]
                #bbox_to_anchor=(.5, .3, .76, .5),           # GHZ, [Nk, m]
                #bbox_to_anchor=(.58, .3, .72, .5),           # GHZ, diff m
                #bbox_to_anchor=(.5, .33, .72, .5),           # Had, diff m
                #bbox_to_anchor=(.5, .18, .72, .5),           # rand, diff m
                bbox_transform=axNow.transAxes, loc=3)

        mark_inset(axNow, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

    elif method == 1:

        ax2 = plt.axes([0,0,1,1])


        #ip  = InsetPosition(axNow, [0.55, 0.3, 0.4, 0.3])
        ip  = InsetPosition(axNow, [0.55, 0.3, 0.4, 0.3])

        ax2.set_axes_locator(ip)
        mark_inset(axNow, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

    for kk, (RunT, Target_Err_Xk, label) in enumerate(zip(RunT_wk, Target_Err_wk, label_list)):
        #if kk > 1:
        #    break
        ax2.plot(RunT, Target_Err_Xk, marker=mk_list[kk], linestyle=ln_list[kk],
                label='{}'.format(label))

    #ax2.legend(loc=0)
    #ax2.set_yticks(np.arange(0.0, 0.1, 0.01))
    ax2.set_yticks(yTick)

    if len(xTick) > 0: 
        ax2.set_xticks(xTick)

    #xRng  = [3.5, 8.0]
    #yRng  = [0.028, 0.045]

    #ax2.set_xlim([0.3, 0.55])
    #ax2.set_ylim([0.0, 0.1])
    ax2.set_xlim(xRng)
    ax2.set_ylim(yRng) 


def plt_TargetErr_RunT(wk_list, dt_list, label_list, tit=None, combine=[]):

    print('\n ------------------------------------------------------------------- ')
    print(' ##          start to plot Target Err  (versus)  Run Time         ## ')
    print(' ------------------------------------------------------------------- \n')

    if len(combine) == 0:
        fig, axNow = plt.subplots(1, 1, figsize=(8,6))
    else:
        axNow = combine[0]
        lab_x = combine[1]
        lab_y = combine[2]

    T_sum, T_rec, Target_Err_wk, RunT_wk, Final_Err_wk, Converged = calc_RunT(wk_list, dt_list, label_list)

    lnmk_list = ['-+','-.^','--o','--x', ':>','-<',2, 3]

    mk_list = ['+','^','o','x', '>','<',2, 3]
    ln_list = ['-', '-.', '--', '--', ':', '-']

    #for RunT, Target_Err_Xk, label in zip(RunT_wk, Target_Err_wk, label_list):
    for ii, (RunT, Target_Err_Xk, label) in enumerate(zip(RunT_wk, Target_Err_wk, label_list)):

        #plt.plot(RunT, Target_Err_Xk, label='{}'.format(label))
        #plt.title(tit)

        #axNow.plot(RunT, Target_Err_Xk, label='{}'.format(label))
        #axNow.plot(RunT, Target_Err_Xk, lnmk_list[cnt],
        #           label='{}'.format(label))

        #if ii > 0:
        #    break

        axNow.plot(RunT, Target_Err_Xk, marker=mk_list[ii], linestyle=ln_list[ii],
                   label='{}'.format(label))

    axNow.set_title(tit, y=0.8, fontsize=14)

    # ------------------------------ #
    #               inset            #
    # ------------------------------ #
    Plt_inset = 1
    if Plt_inset == 1:

        if lab_x == 0 and lab_y == 1:       #  ax[0, 0] corner upper left

            #xRng  = [33, 65];   yRng  = [0.05, 0.15]      #  GHZ, diff m
            #yTick = np.arange(0.0, 0.18, 0.03)                #  GHZ, diff m

            #xRng  = [18, 65];   yRng  = [0.065, 0.08]      #  Had, diff m
            #yTick = np.arange(0.0, 0.18, 0.01)                #  Had, diff m

            #xRng  = [38, 102];   yRng  = [-0.05, 0.25]      #  rand, diff m
            #yTick = np.arange(-0.1, 0.38, 0.1)                #  rand, diff m
            # ------------------------------------- #

            #xRng  = [0.24, 0.68];   yRng  = [-0.03, 0.27]      #  GHZ, [Nk, m]
            #yTick = np.arange(0.0, 0.3, 0.1)                #  GHZ, [Nk, m]
            #xTick = [0.3, 0.6]
            #call_inset(axNow, RunT_wk, Target_Err_wk, label_list, xRng, yRng, yTick, xTick)

            #xRng  = [0.25, 0.68];   yRng  = [0.017, 0.056]      #  Had, [Nk, m]
            #yTick = np.arange(0.0, 0.1, 0.02)                #  Had, [Nk, m]
            #xTick = [0.3, 0.6]
            #call_inset(axNow, RunT_wk, Target_Err_wk, label_list, xRng, yRng, yTick, xTick)

            xRng  = [-0.05, 2.4];   yRng  = [-0.06, 0.35]      #  rand, [Nk, m]
            yTick = np.arange(0.0, 0.7, 0.2)                #  rand, [Nk, m]

            #xRng  = [275, 380];   yRng  = [-0.01, 0.08]        #  rand, [r]
            #yTick = np.arange(0.0, 0.1, 0.03)                   #  rand, [r]

            call_inset(axNow, RunT_wk, Target_Err_wk, label_list, xRng, yRng, yTick)


    Plt_inset = 1
    if Plt_inset == 1:

        if lab_x == 0 and lab_y == 0:       #  ax[0, 1] corner upper right

            #xRng  = [70, 105];   yRng  = [0.04, 0.08]      #  GHZ, diff m
            #yTick = np.arange(0.0, 0.1, 0.02)                #  GHZ, diff m

            #xRng  = [68, 105];   yRng  = [0.038, 0.115]      #  Had, diff m
            #yTick = np.arange(0.0, 0.1, 0.03)                #  Had, diff m

            #xRng  = [90, 165];   yRng  = [-0.025, 0.155]      #  rand, diff m
            #yTick = np.arange(0.0, 0.2, 0.05)                #  rand, diff m
            # ------------------------------------- #

            #xRng  = [4.7, 7.9];   yRng  = [0.029, 0.041]        #  GHZ, [Nk, m]
            #yTick = np.arange(0.0, 0.1, 0.01)                   #  GHZ, [Nk, m]

            #xRng  = [3.5, 8.5];   yRng  = [0.028, 0.045]        #  Had, [Nk, m]
            #yTick = np.arange(0.0, 0.1, 0.01)                   #  Had, [Nk, m]

            xRng  = [3.9, 13];   yRng  = [-0.04, 0.20]        #  rand, [Nk, m]
            yTick = np.arange(0.0, 0.21, 0.06)                   #  rand, [Nk, m]

            #xRng  = [350, 540];   yRng  = [-0.01, 0.08]        #  rand, [r]
            #yTick = np.arange(0.0, 0.1, 0.03)                   #  rand, [r]
            call_inset(axNow, RunT_wk, Target_Err_wk, label_list, xRng, yRng, yTick)

    Plt_inset = 1
    if Plt_inset == 1:
        if lab_x == 1 and lab_y == 1:       #  ax[1, 0] corner lower left

            #xRng  = [145, 185];   yRng  = [0.025, 0.06]      #  GHZ, diff m
            #yTick = np.arange(0.0, 0.1, 0.03)                #  GHZ, diff m

            #xRng  = [140, 185];   yRng  = [0.03, 0.043]      #  Had, diff m
            #yTick = np.arange(0.0, 0.1, 0.01)                #  Had, diff m

            #xRng  = [180, 270];   yRng  = [-0.02, 0.1]      #  rand, diff m
            #yTick = np.arange(0.0, 0.1, 0.04)                #  rand, diff m
            # ------------------------------------- #

            #xRng  = [140, 185];    yRng  = [0.028, 0.06]    #  GHZ, [Nk, m]
            #yTick = np.arange(0.0, 0.1, 0.02)               #  GHZ, [Nk, m]

            #xRng  = [135, 185];    yRng  = [0.028, 0.05]    #  Had, [Nk, m]
            #yTick = np.arange(0.0, 0.1, 0.01)               #  Had, [Nk, m]
            #call_inset(axNow, RunT_wk, Target_Err_wk, label_list, xRng, yRng, yTick)

            xRng  = [155, 280];    yRng  = [-0.02, 0.10]    #  rand, [Nk, m]
            yTick = np.arange(0.0, 0.13, 0.04)               #  rand, [Nk, m]
            xTick = [180, 260]
            call_inset(axNow, RunT_wk, Target_Err_wk, label_list, xRng, yRng, yTick, xTick)
            
            #xRng  = [415, 610];    yRng  = [-0.009, 0.08]    #  rand, [r]
            #yTick = np.arange(0.0, 0.1, 0.03)               #  rand, [r]
            #xTick = [180, 260]
            #call_inset(axNow, RunT_wk, Target_Err_wk, label_list, xRng, yRng, yTick)
 
    Plt_inset = 1
    if Plt_inset == 1:
        if lab_x == 1 and lab_y == 0:       #  ax[1, 1] corner lower right

            #xRng  = [215, 270];   yRng  = [0.022, 0.045]      #  GHZ, diff m
            #yTick = np.arange(0.0, 0.1, 0.01)                #  GHZ, diff m

            #xRng  = [215, 270];   yRng  = [-0.01, 0.075]      #  Had, diff m
            #yTick = np.arange(0.0, 0.1, 0.03)                #  Had, diff m

            #xRng  = [278, 380];   yRng  = [-0.01, 0.075]      #  rand, diff m
            #yTick = np.arange(0.0, 0.1, 0.03)                #  rand, diff m
            # ------------------------------------- #

            #xRng  = [1600, 2500];    yRng  = [0.06, 0.11]   #  GHZ, [Nk, m]
            #yTick = np.arange(0.0, 0.12, 0.02)               #  GHZ, [Nk, m]

            #xRng  = [1700, 2600];    yRng  = [0.06, 0.08]   #  Had, [Nk, m]
            #yTick = np.arange(0.0, 0.1, 0.02)               #  Had, [Nk, m]
            #call_inset(axNow, RunT_wk, Target_Err_wk, label_list, xRng, yRng, yTick)
            
            xRng  = [2600, 4300];    yRng  = [-0.009, 0.115]   #  rand, [Nk, m]
            yTick = np.arange(0.0, 0.1, 0.04)               #  rand, [Nk, m]
            xTick = [2800, 4200]

            #xRng  = [460, 710];    yRng  = [-0.009, 0.085]   #  rand, [r]
            #yTick = np.arange(0.0, 0.1, 0.03)               #  rand, [r]

            #xTick = [2800, 4200]
            call_inset(axNow, RunT_wk, Target_Err_wk, label_list, xRng, yRng, yTick, xTick)
            #call_inset(axNow, RunT_wk, Target_Err_wk, label_list, xRng, yRng, yTick)


    # ------------------------------------- #
    #               labels                  #  
    # ------------------------------------- #
    #plt.xlabel('  Run time (sec)')
    #plt.ylabel('  || Xk - rho ||_F')
    if len(combine) == 0:

        axNow.set_xlabel('  Run time (sec)')
        #axNow.set_ylabel('  $\|| Xk - rho \||_F$')
        axNow.set_ylabel(r'$\left\Vert X_k -\rho \right\Vert_F$')

    else:
        if lab_x == 1:            
            axNow.set_xlabel('  Run time (sec)')
        if lab_y == 1:
            #axNow.set_ylabel('  || Xk - rho ||_F')
            axNow.set_ylabel(r'$\left\Vert X_k -\rho \right\Vert_F$', fontsize=16)

    if len(combine) == 0:
            axNow.legend()
    else:
        if lab_x == 0 and lab_y == 1:          #  uper left corner ax[0][0]
        #if lab_x == 1 and lab_y == 0:           #  lower right corner ax[1][1]
        #if lab_x == 1 and lab_y == 1:           #  lower left corner ax[1][0]

            #axNow.legend(loc='upper right')
            #axNow.legend(loc='center right', bbox_to_anchor=(0.95, 0.5))
            #axNow.legend(loc='center right', fontsize=13)
            #axNow.legend(loc='lower right', fontsize=13)

            #axNow.legend(loc='center right', fontsize=13)
            #axNow.legend(loc='best', fontsize=13)
            #axNow.legend(bbox_to_anchor=(0.5, 0.15, 0.5, 0.5), fontsize=13)
            #axNow.legend(loc='best', bbox_to_anchor=(0.5, 0.149, 0.5, 0.5), fontsize=13)  # rand, diff m
            #axNow.legend(loc='best', bbox_to_anchor=(0.5, 0.13, 0.5, 0.5), fontsize=13)  # rand, diff r
            #axNow.legend(loc='best',bbox_to_anchor=(0.08, 0.26, 0.5, 0.5), fontsize=13)  # Had, diff m
            #axNow.legend(loc='best',bbox_to_anchor=(0.08, 0.23, 0.5, 0.5), fontsize=13)  # GHZ, diff m

            axNow.legend(loc='upper left', bbox_to_anchor=(-0.25, 1.25), ncol=3, 
                         fontsize=13, borderaxespad=0)  # rand, diff r


            #axNow.set_xlim([-300, 8200])
            #axNow.set_xlim([-300, 6000])

        #if lab_x == 1 and lab_y == 1:           #  lower left corner ax[1][0]
        #    axNow.set_xlim([-50, 1600])
        #    axNow.set_xlim([-50, 2000])         #  rand, diff m
        #    axNow.set_xticks(np.arange(0, 2500, 500))  #  GHZ, [Nk, m]
        #    axNow.set_xlim([-100, 1950])                #  GHZ, [Nk, m]

        #if lab_x == 0 and lab_y == 1:          #  uper left corner ax[0][0]
        #    axNow.set_xlim([-0.08, 2.1])         #  Had, [Nk, m]
        #    axNow.set_xlim([-0.4, 69])          #  Had, diff m

    #plt.legend(loc='upper right')
    plt.plot()


def Get_Time_Err_diff_parm(parm_List, Data_pm, RG_info, mu_List, note):
    """     RG_info = InitX, Md_tr, Md_alp, Md_sig

    """

    print('\n ------------------------------------------------------------------- ')
    print(' ##          start to plot convergence Time (versus) m       ## ')
    print(' ------------------------------------------------------------------- \n')

    Option=3

    T_sum_dict     = {}
    T_rec_dict     = {}
    Final_Err_dict = {}
    note_dict      = {}

    # ------------  RGD  ---------------------------- #
    #pm_method = 'RGD', InitX, Md_tr, Md_alp, Md_sig
    #wk_list, dt_list, label_list = Get_wk_diff_m_single_method(m_List, Data_pm, pm_method)

    wk_list, dt_list, pm_list, label_list = Get_wk_RGD_diff_param(parm_List, note, Data_pm, RG_info, Option)

    T_sum, T_rec, Target_Err_wk, RunT_wk, Final_Err_wk, Converged = calc_RunT(wk_list, dt_list, label_list)

    note_RGD = 'RGD: {}'.format(note)

    T_sum_dict['RGD']     = T_sum
    T_rec_dict['RGD']     = T_rec
    Final_Err_dict['RGD'] = Final_Err_wk
    note_dict['RGD']      = note_RGD

    print(' *************************************** ')
    print(' ***      Now working for MiFGD      *** ')
    print(' *************************************** ')

    # ------------  MiFGD ------------------------- #
    for mu in mu_List:
        #pm_method = 'MiFGD', mu
        #wk_list, dt_list, label_list = Get_wk_diff_m_single_method(m_List, Data_pm, pm_method)

        wk_list, dt_list, pm_list, label_list = Get_wk_MiFGD_mu_diff_param(parm_List, Data_pm, mu, Option)

        T_sum, T_rec, Target_Err_wk, RunT_wk, Final_Err_wk, Converged = calc_RunT(wk_list, dt_list, label_list)
        note_mu = 'MiFGD mu = {}: {}'.format(mu, note)

        T_sum_dict[mu]     = T_sum
        T_rec_dict[mu]     = T_rec
        Final_Err_dict[mu] = Final_Err_wk
        note_dict[mu]      = note_mu

    return T_sum_dict, T_rec_dict, Final_Err_dict, note_dict


def Plt_convergence_Time_Err_diff_parm(parm_List, T_sum_dict, T_rec_dict, Final_Err_dict, note_dict):

    print(' keys = {}'.format(T_sum_dict.keys()))

    fig, ax = plt.subplots(1, 2, figsize=(8,6))

    for key in T_sum_dict.keys():

        T_sum        = T_sum_dict[key]  
        T_rec        = T_rec_dict[key] 
        Final_Err_wk = Final_Err_dict[key]
        note_label   = note_dict[key] 

        plt.subplot(1,2,1)
        plt.plot(parm_List, T_sum, label='{}, Tsum'.format(note_label))
        plt.plot(parm_List, T_rec, label='{}, Trec'.format(note_label))

        plt.subplot(1,2,2)
        plt.plot(parm_List, Final_Err_wk, label='{}, Err'.format(note_label))


    tit = '{}-{}, Nr = {}, shot = {}'.format(StateName, Nk, Nr, mea)

    plt.subplot(1,2,1)

    plt.legend()
    plt.ylabel('Time for convergence')
    plt.xlabel(' m = num_labels ')
    plt.title(tit)

    plt.subplot(1,2,2)
    plt.legend()

    plt.ylabel('Final Error')
    plt.xlabel(' m = num_labels ')


    plt.plot()



def update_DataPm(Data_pm, elem, Name_List, note, label_list):
        """         this function not working for changing locals()  -->  don't know why??
        """

        print('  elem = {}   ------------'.format(elem))
        label = 'RGD: {}'.format(note)

        Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method = Data_pm

        Nk = 8

        for ii, name in enumerate(Name_List):
            print('     ii = {}, var = {},  value = {}'.format(ii, name, elem[ii]))
            label = '{}, {} = {}'.format(label, name, elem[ii])

            print("         globals()['{}']  = {}".format(name, globals()[name]))
            print("         locals()['{}']  = {}".format(name, locals()[name]))

            Nk  = 9
            exec("{} = {}".format(name, elem[ii]))    
            locals()[name]  = elem[ii]
            globals()[name] = elem[ii]

            print("         globals()['{}']  = {}".format(name, globals()[name]))
            print("         locals()['{}']  = {}".format(name, locals()[name]))

            print('  hex(id(  Nk  ))  = {}'.format(hex(id(Nk))))
            print('  hex(id( locals()[name]  ))  = {}'.format(hex(id( locals()[name] ))))
            print('  hex(id( globals()[name]  ))  = {}'.format(hex(id( globals()[name] ))))

            print('     ----- (loop)  Nk = {}, m = {}, version = {}\n'.format(Nk, m, version))      #  just example

        label_list.append(label)
        print('\n          label = {}'.format(label))
        print('    ***********   Nk = {}, m = {}, version = {}\n'.format(Nk, m, version))      #  just example

        print('     m        = {}'.format(m))
        print('     version  = {}'.format(version))

        Data_pm = [Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method]
        print('     Data_pm = {}'.format(Data_pm))

        Setup = Data_Info(Data_pm)

        return Setup, label_list


def Option_update_params(Data_pm, elem, Option):
        
        Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method = Data_pm

        if Option == 1:                       #  do nothing
            
            label = None
            print('  ****    NO update any parameters    ****   ')

        elif Option == 2:                     #  update   [mea]

            mea = elem

            label = 'mea = {}'.format(mea)
            print('      mea = {}\n'.format(mea))      #  just example

        elif Option == 3:                     #  update   [m]

            m = elem

            label = 'm = {}'.format(m)
            #label = ''
            print('      m = {}\n'.format(m))      #  just example

        elif Option == 5:                     #  update   [Nk, m]

            Nk, m = elem
            label = 'Nk = {}, m = {}'.format(Nk, m)
            print('      Nk = {}, m = {}\n'.format(Nk, m))      #  just example

        elif Option == 6:                     #  update   [Nk, m, version]

            Nk, m, version = elem
            label = 'Nk = {}, m = {}, version = {}'.format(Nk, m, version)
            #print('      Nk = {}, m = {}, version = {}\n'.format(Nk, m, version))      #  just example


        elif Option == 10:                     #  update   [StateName, m, version]

            StateName, m, version, StVer = elem
            label = 'StateName = {}, m = {}, version = {}, StVer={}'.format(StateName, m, version, StVer)
            #print('      Nk = {}, m = {}, version = {}\n'.format(Nk, m, version))      #  just example


        elif Option == 7:                     #  update   [Nr]

            Nr = elem
            label = 'Nr = {}'.format(Nr)

        elif Option == 8:                     #  update   [Nk]

            Nk = elem
            label = 'Nk = {}'.format(Nk)

        Data_pm = [Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method]

        print('      label = {}'.format(label))
        #print('      Data_pm = {}'.format(Data_pm))

        return Data_pm, label


def Get_wk_RGD_diff_param(parm_List, note, Data_pm, RG_info, Option):

    InitX, Md_tr, Md_alp, Md_sig = RG_info

    print(' *************************************** ')
    print(' ***      Now working for RGD        *** ')
    print(' *************************************** ')

    wk_list = []
    dt_list = []
    pm_list = []
    label_list = []

    for elem in parm_List:
        print('\n  -------------  elem = {}   ------------'.format(elem))

        Data_pm, label = Option_update_params(Data_pm, elem, Option)
        Setup = Data_Info(Data_pm)
        print('      Data_pm = {}'.format(Data_pm))

        if label == None:
            label = 'RGD: {}'.format(note)
        else:
            #label = 'RGD: {}, {}'.format(note, label)
            #label = 'RGD: {}'.format(label)
            label = 'RGD'

        label_list.append(label)

        # -----------------     RGD     --------------------- #
        pm_RGD = 'RGD', InitX, Md_tr, Md_alp, Md_sig
        wk_list, dt_list, pm_list = Get_wk_list(Setup, version, wk_list, dt_list, pm_list, pm_RGD)

    return wk_list, dt_list, pm_list, label_list

def Get_wk_MiFGD_mu_diff_param(parm_List, Data_pm, mu, Option):

    wk_list = []
    dt_list = []
    pm_list = []
    label_list = []

    pm_MiF = 'MiFGD', mu
    print('\n ***      Now working for MiFGD   mu = {}   *** \n'.format(mu))
    
    # ----------------   different parameters  ---------------- #
    for elem in parm_List:
        print('\n  -------------  elem = {}   ------------'.format(elem))

        Data_pm, label_0 = Option_update_params(Data_pm, elem, Option)
        Setup = Data_Info(Data_pm)
        print('      Data_pm = {}'.format(Data_pm))

        wk_list, dt_list, pm_list = Get_wk_list(Setup, version, wk_list, dt_list, pm_list, pm_MiF)

        if label_0 == None:
            label = 'MiFGD mu = {}'.format(mu)
        else:
            #label  = 'MiFGD mu = {}, {}'.format(mu, label_0)
            label  = 'MiFGD $\mu$ = {}'.format(mu)

        label_list.append(label)

    return wk_list, dt_list, pm_list, label_list


def Get_wk_MiFGD_diff_param(parm_List, note, Data_pm, mu_List, Option):

    # ------------------------------ #
    #           for MiFGD            #
    # ------------------------------ #

    print(' *************************************** ')
    print(' ***      Now working for MiFGD      *** ')
    print(' *************************************** ')

    wk_list = []
    dt_list = []
    pm_list = []
    label_list = []

    # -----------------     MiFGD   ---------------------- #
    for mu in mu_List:
        wk_list_1, dt_list_1, pm_list_1, label_list_1 = Get_wk_MiFGD_mu_diff_param(parm_List, Data_pm, mu, Option)

        wk_list    += wk_list_1
        dt_list    += dt_list_1
        pm_list    += pm_list_1
        label_list += label_list_1

    return wk_list, dt_list, pm_list, label_list


def Get_wk_diff_parameters(parm_List, note, Data_pm, mu_List, RG_info, Option=6, method_List=[]):

    #Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method = Data_pm

    if len(method_List) == 0:
        wk_list_1, dt_list_1, pm_list_1, label_list_1 = \
            Get_wk_RGD_diff_param(parm_List, note, Data_pm, RG_info, Option)

        wk_list_2, dt_list_2, pm_list_2, label_list_2 = \
            Get_wk_MiFGD_diff_param(parm_List, note, Data_pm, mu_List, Option)
            #Get_wk_MiFGD_diff_param_then_mu(parm_List, note, Data_pm, mu_List, Option)

        wk_list    = wk_list_1    + wk_list_2
        dt_list    = dt_list_1    + dt_list_2
        pm_list    = pm_list_1    + pm_list_2
        label_list = label_list_1 + label_list_2

        return wk_list, dt_list, label_list 
    
    else:
        wk_list    = []
        dt_list    = []
        pm_list    = []
        label_list = []
        Case_list  = []

        for ii, method in enumerate(method_List):

            if method[0] == 'RGD':

                InitX, Md_tr, Md_alp, Md_sig = RG_info
                InitX = method[1]
                RG_info = [InitX, Md_tr, Md_alp, Md_sig]

                print('\n  *********  {}-th method: {} --> RG_info = {}\n'.format(ii, method, RG_info))

                wk_list_0, dt_list_0, pm_list_0, label_list_0 = \
                    Get_wk_RGD_diff_param(parm_List, note, Data_pm, RG_info, Option)

                #label_list_0 = ['{}: InitX = {}'.format(label, InitX) for label in label_list_0]

                if InitX == 0:
                    label_list_0 = ['{}: random $X_0$'.format(label) for label in label_list_0]
                elif InitX == 1:
                    #label_list_0 = [r'$X_0 = \mathcal{H}_r(\mathcal{A}^\dagger(y))$' for label in label_list_0]

                    X0tex = r'$X_0 = \mathcal{H}_r(\mathcal{A}^\dagger(y))$'
                    label_list_0 = ['{}: '.format(label)+X0tex for label in label_list_0]


                CaseNow = [[method, pm] for pm in parm_List]

            elif method[0] == 'MiFGD':

                print('\n  *********  {}-th method: {} --> MiFGD InitX = {}\n'.format(ii, method, method[1]))

                wk_list_0, dt_list_0, pm_list_0, label_list_0 = \
                    Get_wk_MiFGD_diff_param(parm_List, note, Data_pm, mu_List, Option)

                mu_pm = [[mu, pm] for mu in mu_List for pm in parm_List]
                CaseNow = [[method, xx] for xx in mu_pm]

            Case_list  += CaseNow

            wk_list    += wk_list_0
            dt_list    += dt_list_0
            pm_list    += pm_list_0
            label_list += label_list_0

        return wk_list, dt_list, label_list, Case_list 
    


def Data_Info(argv):

    Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method = argv
    
    print('\n     ----------------------      Data_Info  (To find Setup = Data Info)     ------------------------------ ')    
    print('        argv = {}'.format(argv))
    print('            Nk = {}, StateName = {}, m = {}, mea = {}, Nr = {}, version = {}'.format(Nk, StateName, m, mea, Nr, version))
    print('            New_Pj_shot = {}, StVer = {}, Pj_method = {}, mea_method = {}, measure_method = {}\n'.format(
                New_Pj_shot, StVer, Pj_method, mea_method, measure_method))

    # ----------------------------- #
    #   to load  each worker        #
    # ----------------------------- #

    Dir, Name, params_setup, version, proj_path, meas_path, Dir0, StVer = \
            State_Naming(Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method)

    ver_Prj = version[0]
    ver_mea = version[1]
    zModel  = version[2]


    Setup = Dir
    #print(Setup)

    if zModel < 0:          #  the normal shot measurements
        #Setup = 'data/{}/{}_m{}_s{}_shot{}_v{}'.format(StName, StName, m, ver_Prj, mea, ver_mea)
        #Setup = 'data/{}/{}_m{}_s{}/{}_m{}_s{}_shot{}_v{}'.format(StName, StName, m, ver_Prj, 
        #                                                          StName, m, ver_Prj, mea, ver_mea)
        Setup = '{}_shot{}_v{}'.format(Setup, mea, ver_mea)


    else:                       #  the noize model
        #Setup = 'data/{}/{}_m{}_s{}_zN{}_v{}'.format(StName, StName, m, ver_Prj, zModel, ver_mea)
        #Setup = 'data/{}/{}_m{}_s{}/{}_m{}_s{}_zN{}_v{}'.format(StName, StName, m, ver_Prj, 
        #                                                        StName, m, ver_Prj, zModel, ver_mea)
        Setup = '{}_zN{}_v{}'.format(Setup, zModel, ver_mea)

    print('     Setup = {}'.format(Setup))
    print('     *******************************************************  \n')
    return Setup


def Get_wk_list(Setup, version, wk_list, dt_list, pm_list, param_argv):
        
        worker, dt, params_dict = method_2_worker(Setup, version, *param_argv)
        #worker, dt, params_dict = method_2_worker(Setup, version, 'RGD', InitX, Md_tr, Md_alp, Md_sig)
        #worker, dt, params_dict = method_2_worker(Setup, version, 'MiFGD', mu)

        wk_list.append(worker)
        dt_list.append(dt)
        pm_list.append(params_dict)

        return wk_list, dt_list, pm_list


def method_2_worker(Setup, version, *method_argv):

    method = method_argv[0]
    zModel = version[2]

    print('\n      arvg = {}'.format(method_argv))
    print('          method = {}'.format(method))

    if method == 'MiFGD':
        mu    = method_argv[1]

        Fname = '{}-MiFGD-mu{:.3}'.format(Setup, mu)
        label = '{}, mu = {:.3}'.format(method, mu)

    elif method == 'RGD':    
        InitX, Md_tr, Md_alp, Md_sig = method_argv[1:]
        print('      RGD:  (InitX, Md_tr, Md_alp, Md_sig) = ({}, {}, {}, {})'.format(InitX, Md_tr, Md_alp, Md_sig))

        Fname = '{}-RGD_Ix{}_Tr{}_Ap{}_sg{}'.format(Setup, InitX, Md_tr, Md_alp, Md_sig)

    # --------------------------------------------- #
    #   directly stored worker or container         #
    # --------------------------------------------- #

    if os.path.exists('{}_wrapper.dat'.format(Fname)):
        Fname = '{}_wrapper.dat'.format(Fname)

    elif os.path.exists('{}.dat'.format(Fname)):
        Fname = '{}.dat'.format(Fname)

    print("          method = {}  --> Fname = {}".format(method, Fname))

    # --------------------------------- #
    #           loading worker          #
    # --------------------------------- #
    with open(Fname, 'rb') as f:
        params_dict, methodRec, worker, dt = pickle.load(f)
        
        if methodRec != method:
            print(' Error:  Not the called method')

        #if m != worker.num_labels:
        #    print(' m (the num of label list) is not correct')

        #if zModel < 0:          #  the normal shot measurements
        #    if mea != worker.num_shots:
        #        print(' mea (the number of shots) is not correct')

    print("                  --> Fname = {} has been loaded (DONE)".format(method, Fname))

    return worker, dt, params_dict


if __name__ == "__main__":

    # ----------------------------- #
    #   basic   setting             #
    # ----------------------------- #

    #Nk, m, mea = 6, 819, 2400
    #Nk, m, mea = 6, 819, 8600
    #Nk, m, mea = 8, 6553, 8600          #  10% = 6553,  20% = 13107
    Nk, m, mea = 10, 314572, 8600           #  for rand ->  mea NOT matter
    #Nk, m, mea = 12, 838860, 8600

    StateName = 'rand'
    if StateName == 'rand':
        version = [1, 1, 0]   # [proj version s, measure version, noise model]
        StVer = [1, 0]
        Nr = 3
    else:
        version = [1, 1, -1]   # [proj version s, measure version, noise model]
        StVer = 0
        Nr = 1

    # ------------------------- #
    #   for producing data      #
    # ------------------------- #
    Pj_method = 1                   #   the method to save | load  projectors
    mea_method = 1                  #   the method to save | load  measurement_dict (count_dict) 

    measure_method = 1              #  = 1: direct label_list,  = 3: parallel cpu

    New_Pj_shot = [0, 0]  #  [New Proj?, New shot?]  | [0,0] loading | 

    Data_pm = [Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method]
    # ------------------------ #
    #       RGD / MiFGD        #
    # ------------------------ #

    #method = 'RGD'
    InitX  = 1                  #   method of choosing initial X0
    Md_tr  = 0                  #   Method if including trace = 1
    Md_alp = 0                  #   method for scaling alpha
    Md_sig = 0                  #   method for scaling singular value

    RG_info = InitX, Md_tr, Md_alp, Md_sig


    Show_Case = 52
    if Show_Case == -1:         #  single case test
        #Setup = 'data/GHZ-6/GHZ-6_m819_s1/GHZ-6_m819_s1_shot3600_v1'
        Setup = 'data/GHZ-12/GHZ-12_m838860_s3/GHZ-12_m838860_s3_shot8600_v1'

        version = [1, 1, -1]

        pm_RGD = 'RGD', InitX, Md_tr, Md_alp, Md_sig
        worker, dt, param_dict = method_2_worker(Setup, version, *pm_RGD)

        plt_TargetErr_RunT([worker], [dt], ['RGD'])

    elif Show_Case == 1:
        Null_List = [0]
        mu_List   = [0.75]
        note      = '{}'.format(StateName)

        #wk_list, dt_list, label_list = Get_wk_RGD_MiFGD(InitX, Md_tr, Md_alp, Md_sig, mu_List, Data_pm) 
        wk_list, dt_list, label_list = Get_wk_diff_parameters(Null_List, note, Data_pm, mu_List, RG_info, Option=1)
        
        plt_TargetErr_RunT(wk_list, dt_list, label_list)

    elif Show_Case == 2:                #  different mea (shot number)

        #mea_List = [1200, 2400, 3600, 4800]
        mea_List = [2400, 3600]
        mu_List  = [0.75]
        note      = '{}'.format(StateName)

        #wk_list, dt_list, label_list = Get_wk_diff_shots(mea_List, InitX, Md_tr, Md_alp, Md_sig, mu_List, Data_pm)        
        wk_list, dt_list, label_list = Get_wk_diff_parameters(mea_List, note, Data_pm, mu_List, RG_info, Option=2)

        tit = '{}-{}, Nr = {}, m = {}'.format(StateName, Nk, Nr, m)
        plt_TargetErr_RunT(wk_list, dt_list, label_list, tit)

    elif Show_Case == 4:                # different m  -->  Time & Err

        #m_List  = [409, 819, 1638]      # for Nk = 6
        #m_List  = [3276, 6553, 13107, 19660]    #  for Nk = 8
        m_List  = [3276, 6553, 13107]    #  for Nk = 8

        mu_List = [0.75]
        note      = '{}'.format(StateName)


        T_sum_dict, T_rec_dict, Final_Err_dict, note_dict = Get_Time_Err_diff_parm(m_List, Data_pm, RG_info, mu_List, note)

        Plt_convergence_Time_Err_diff_parm(m_List, T_sum_dict, T_rec_dict, Final_Err_dict, note_dict)


    elif Show_Case == 3:                #  different m

        #m_List  = [409, 819, 1638]             #  for Nk = 6
        #m_List  = [3276, 6553, 13107, 19660]    #  for Nk = 8
        #m_List  = [3276, 6553, 13107]    #  for Nk = 8
        #m_List  = [52428, 104857, 209715,314572]    #  for Nk = 10
        #m_List  = [52428, 104857, 209715]    #  for Nk = 10
        m_List  = [104857, 209715]    #  for Nk = 10

        #m_List  = [3276]    #  for Nk = 8
        #m_List  = [104857]    #  for Nk = 10
        #m_List  = [209715]    #  for Nk = 10
        #m_List  = [838860]

        #mu_List = [0.75]
        mu_List = [0.75, 0.5, 1/3]
        #mu_List = [0.75, 0.5, 1/3, 0.2]
        #mu_List = [0.5, 1/3, 0.2]
        #mu_List = [0.75, 0.5, 1/3, 0.2, 0.25]
        #mu_List = []

        method_List = [['RGD', 1], ['RGD', 0], ['MiFGD', 0]]               # [method, InitX]
        #method_List = [['RGD', 1], ['RGD', 0]]               # [method, InitX]
        #method_List = [['RGD', 1], ['MiFGD', 0]]               # [method, InitX]
        #method_List = []

        note    = '{}'.format(StateName)

        #wk_list, dt_list, label_list = Get_wk_diff_m(m_List, InitX, Md_tr, Md_alp, Md_sig, mu_List, Data_pm)            #  direct version (old)
        wk_list, dt_list, label_list, Case_list = Get_wk_diff_parameters(m_List, note, Data_pm, mu_List, RG_info, 3, method_List)

        tit = '{}-{}, Nr = {}, shot = {}'.format(StateName, Nk, Nr, mea)

        #fig, ax = plt.subplots(2,1, figsize=(8,6))
        #plt_TargetErr_RunT(wk_list, dt_list, label_list, tit, [ax[1], 1,1])
        plt_TargetErr_RunT(wk_list, dt_list, label_list, tit)

        plt.legend()
        #plt.xlim([-20, 2500])
        plt.plot()

        List_runT, List_T_RGD, List_T_MiF, Err_R0, Err_R1, Err_MiF, List_conv, List_Iter, List_note = \
            Compare_RunT_work(wk_list, label_list, Case_list, mu_List)

        #fig, ax11 = plt.subplots(1,1, figsize=(8,6))
        #Plt_scaling_compare(m_List, List_T_RGD, List_T_MiF, Err_R0, Err_R1, Err_MiF, ax11)
        Plt_scaling_compare(m_List, List_T_RGD, List_T_MiF, Err_R0, Err_R1, Err_MiF)


    elif Show_Case == 31:                #  different m

        mu_List = [0.75]
        method_List = [['RGD', 1], ['RGD', 0], ['MiFGD', 0]]               # [method, InitX]
        #method_List = [['RGD', 1], ['RGD', 0]]               # [method, InitX]
        #method_List = [['RGD', 1], ['MiFGD', 0]]               # [method, InitX]
        #method_List = []

        #m_List  = [3276, 6553, 13107, 19660]    #  for Nk = 8
        #m_List  = [838860]

        #note    = '{}'.format(StateName)
        note    = ''
        #tit = '{}-{}, Nr = {}, shot = {}'.format(StateName, Nk, Nr, mea)
        #tit0 = '{}-{}'.format(StateName, Nk)
        tit0 = '{}({})'.format(StateName, Nk)

        plt.rcParams.update({'font.size':14})
        fig,  ax  = plt.subplots(2,2, figsize=(8,6))
        
        method = 2
        if method == 1:
            #m_List  = [3276]    #  for Nk = 8
            m_List  = [52428]    #  for Nk = 10
            tit = '(a) {}, m = {}'.format(tit0, m_List[0])
            wk_list, dt_list, label_list, Case_list = \
                Get_wk_diff_parameters(m_List, note, Data_pm, mu_List, RG_info, 3, method_List)
            plt_TargetErr_RunT(wk_list, dt_list, label_list, tit, [ax[0,0], 0, 1])

            #m_List  = [6553]    #  for Nk = 8
            m_List  = [104857]    #  for Nk = 10
            tit = '(b) {}, m = {}'.format(tit0, m_List[0])
            wk_list, dt_list, label_list, Case_list = \
                Get_wk_diff_parameters(m_List, note, Data_pm, mu_List, RG_info, 3, method_List)
            plt_TargetErr_RunT(wk_list, dt_list, label_list, tit, [ax[0,1], 0, 0])

            #m_List  = [13107]    #  for Nk = 8
            m_List  = [209715]    #  for Nk = 10
            tit = '(c) {}, m = {}'.format(tit0, m_List[0])
            wk_list, dt_list, label_list, Case_list = \
                Get_wk_diff_parameters(m_List, note, Data_pm, mu_List, RG_info, 3, method_List)        
            plt_TargetErr_RunT(wk_list, dt_list, label_list, tit, [ax[1,0], 1, 1])

            #m_List  = [19660]    #  for Nk = 8
            m_List  = [314572]    #  for Nk = 10
            tit = '(d) {}, m = {}'.format(tit0, m_List[0])
            wk_list, dt_list, label_list, Case_list = \
                Get_wk_diff_parameters(m_List, note, Data_pm, mu_List, RG_info, 3, method_List)
            plt_TargetErr_RunT(wk_list, dt_list, label_list, tit, [ax[1,1], 1, 0])

        elif method == 2:
            figCap = ['(a)', '(b)', '(c)', '(d)']

            #m_List  = [3276, 6553, 13107, 19660]    #  for Nk = 8
            m_List  = [52428, 104857, 209715,314572]    #  for Nk = 10

            for ii, m in enumerate(m_List):
                px = ii%2
                py = int(ii/2)
                lbx = py
                lby = (px+1)%2      # to label y or not
                print('py = {}, px = {}, lby = {}, lbx = {}'.format(py, px, lby, lbx))
                tit = '{} {}, m = {}'.format(figCap[ii], tit0, m)
                wk_list, dt_list, label_list, Case_list = \
                    Get_wk_diff_parameters([m], note, Data_pm, mu_List, RG_info, 3, method_List)
                plt_TargetErr_RunT(wk_list, dt_list, label_list, tit, [ax[py,px], lbx, lby])


        #plt.legend()
        #plt.xlim([-20, 2500])
        
        FigName = './fig/Fig_{}_diff_m.pdf'.format(StateName)
        plt.savefig(FigName)
        #plt.plot()


    elif Show_Case == 52:                # different [Nk, m]

        mu_List   = [0.75]
        #mu_List   = [0.75, 0.5]

        method_List = [['RGD', 1], ['RGD', 0], ['MiFGD', 0]]               # [method, InitX]
        #method_List = [['RGD', 1]]               # [method, InitX]

        #method_List = []

        note      = '{}'.format(StateName)

        #Nk_m_List = [[6, 819], [8, 3276], [10, 52428], [12, 838860]]     #  [Nk, m] 
        #Nk_m_List = [[6, 819], [8, 13107], [10, 209715], [12, 838860]]     #  [Nk, m] 

        Nk_m_List = [[6, 1228], [8, 13107], [10, 209715], [12, 838860]]     #  [Nk, m] 
        #Nk_m_List = [[10, 209715]]     #  [Nk, m] 

        Nk_list   = [Nk for Nk, m in Nk_m_List]

        figCap = ['(a)', '(b)', '(c)', '(d)']

        plt.rcParams.update({'font.size':14})
        fig, ax = plt.subplots(1,3, figsize=(10,5))
        #fig, ax = plt.subplots(1,3, figsize=(8,6))

        fig.tight_layout()
        plt.subplots_adjust(left=0.1)



        Name_List = ['GHZ', 'Had', 'rand']
        #Name_List = ['GHZ']

        ALL_conv_Collect  = []
        ALL_Iter_Collect  = []
        #ALL_conv_Collect  = [[] for nn in Name_List]
        #ALL_Iter_Collect  = [[] for nn in Name_List]

        for ss, StateName in enumerate(Name_List):

            if StateName == 'rand':
                version = [1, 1, 0]   # [proj version s, measure version, noise model]
                StVer = [1, 0]
                Nr = 3
            else:
                version = [1, 1, -1]   # [proj version s, measure version, noise model]
                StVer = 0
                Nr = 1

            Data_pm = [Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method]

            # ------------------------------------- #
            #   collect data for this StateName     #
            # ------------------------------------- #

            List_T_RGD = [[], []]               #  Run Time for RGD   [InitX = 0, InitX = 1]
            List_T_MiF = [[] for ii in range(len(mu_List))]     #  each elem [] for each mu
        
            List_runT  = []
            List_Err   = []
            #List_Iter  = []
            List_conv  = []      #  converged or not
            List_meth  = []      #  RGD --> InitX = 0 or 1   |  MiFGD --> mu
            List_note  = []      #  elem = [RGD, 0], [RGD, 1], [MiFGD, mu]

            ALL_Iter_RGD = [[], []]  #  ALL_Iter = [[RGD0, RGD1], [mu1, mu2, ...]]
            ALL_conv_RGD = [[], []]  #  ALL_conv = [[RGD0, RGD1], [mu1, mu2, ...]]
            ALL_Iter_MiF = [[] for mu in mu_List]
            ALL_conv_MiF = [[] for mu in mu_List]


            Err_R0     = []
            Err_R1     = []
            Err_MiF    = [[] for mu in mu_List]

            for ii, Nk_m in enumerate(Nk_m_List):
                #Nk, m = Nk_m

                wk_list, dt_list, label_list, Case_list = \
                    Get_wk_diff_parameters([Nk_m], note, Data_pm, mu_List, RG_info, 5, method_List)

                # --------------------------------------- #
                #   collect all parm_List cases to plot   #
                # --------------------------------------- #

                Now_runT, Now_T_RGD, Now_T_MiF, e_R0, e_R1, e_MiF, Now_conv, Now_Iter, Now_note = \
                    Compare_RunT_work(wk_list, label_list, Case_list, mu_List)

                List_runT.append(Now_runT)

                List_T_RGD[0].append(Now_T_RGD[0][0])      #  RGD, InitX = 0
                List_T_RGD[1].append(Now_T_RGD[1][0])      #  RGD, InitX = 1
                for kk, mu in enumerate(mu_List):
                    #List_T_MiF[kk].append(Now_T_MiF[kk][0])
                    List_T_MiF[kk] += Now_T_MiF[kk]

                #Err_R0.append(e_R0)         #  target Error for RGD, InitX=0
                #Err_R1.append(e_R1)         #  target Error for RGD, InitX=1
                Err_R0.append(e_R0[0])         #  target Error for RGD, InitX=0
                Err_R1.append(e_R1[0])         #  target Error for RGD, InitX=1
                for kk, mu in enumerate(mu_List):
                    #Err_MiF[kk].append(e_MiF[kk])
                    Err_MiF[kk].append(e_MiF[kk][0])

                # ----------------------------------- #
                #   collect the iterations | conv     #
                # ----------------------------------- #
                ALL_Iter_RGD[0] += Now_Iter[0][0]
                ALL_Iter_RGD[1] += Now_Iter[0][1]

                ALL_conv_RGD[0] += Now_conv[0][0]
                ALL_conv_RGD[1] += Now_conv[0][1]

                for kk, mu in enumerate(mu_List):
                    ALL_Iter_MiF[kk] += Now_Iter[1][kk]
                    ALL_conv_MiF[kk] += Now_conv[1][kk]

            # ----------------------------------------- #
            #   done for this StateName  -->  to plot   #
            # ----------------------------------------- #
            ALL_Iter = [ALL_Iter_RGD, ALL_Iter_MiF]
            ALL_conv = [ALL_conv_RGD, ALL_conv_MiF]

            ALL_Iter_Collect.append(ALL_Iter)
            ALL_conv_Collect.append(ALL_conv)
            #ALL_Iter_Collect[ss].append(ALL_Iter)
            #ALL_conv_Collect[ss].append(ALL_conv)

            #Plt_scaling_compare(Nk_list, List_T_RGD, List_T_MiF, Err_R0, Err_R1, Err_MiF, \
            #                    ALL_Iter, ALL_conv)

            Plt_scaling_compare(Nk_list, List_T_RGD, List_T_MiF, Err_R0, Err_R1, Err_MiF, \
                                ALL_Iter, ALL_conv, ax[ss])

            tit0 = '{} {}'.format(figCap[ss], StateName)
            #ax[ss].set_ylim([-100, 12000])            
            #ax[ss].set_title(tit0, y = 0.8)
            ax[ss].text(x=0.1, y=0.9, s=tit0, transform=ax[ss].transAxes, 
                        fontsize=16)


            # -------------------------------- #
            #   labeling for each subplot      #
            # -------------------------------- #
            ax[ss].set_xlabel('qubit number $n$')
            
            if ss == 0:
                ax[ss].set_ylabel('computation time (sec)') 

            if ss == 0:
                ax[ss].legend(loc='center right',fontsize=12)

        #plt.xlim([-20, 4100])
        #plt.ylim([-50, 25000])

        #FigName = './fig/Fig_scaling.pdf'
        FigName = './fig/Fig_scaling.eps'
        plt.savefig(FigName)
        #plt.plot()


    elif Show_Case == 51:                # different [Nk, m]

        mu_List   = [0.75]
        method_List = [['RGD', 1], ['RGD', 0], ['MiFGD', 0]]               # [method, InitX]
        #method_List = []

        note      = '{}'.format(StateName)

        #Nk_m_List = [[6, 819], [8, 3276], [10, 52428], [12, 838860]]     #  [Nk, m] 
        #Nk_m_List = [[6, 819], [8, 13107], [10, 209715], [12, 838860]]     #  [Nk, m] 

        Nk_m_List = [[6, 1228], [8, 13107], [10, 209715], [12, 838860]]     #  [Nk, m] 

        Nk_list   = [Nk for Nk, m in Nk_m_List]

        figCap = ['(a)', '(b)', '(c)', '(d)']
        tit0 = '{}'.format(StateName)

        plt.rcParams.update({'font.size':14})
        #fig, ax = plt.subplots(2,2, figsize=(8,6))
        fig, ax = plt.subplots(2,2, figsize=(8,6))

        method = 2
        if method == 2:

            for ii, Nk_m in enumerate(Nk_m_List):
                Nk, m = Nk_m

                px = ii%2
                py = int(ii/2)
                lbx = py
                lby = (px+1)%2      # to label y or not
                print('py = {}, px = {}, lby = {}, lbx = {}'.format(py, px, lby, lbx))
                print('  Nk_m = {}, Nk = {}, m = {}'.format(Nk_m, Nk, m))

                #tit = ' {} {} $k$ = {}, $m$ = {}'.format(figCap[ii], tit0, Nk, m)
                tit = ' {} {}$({})$, $m$ = {}'.format(figCap[ii], tit0, Nk, m)


                wk_list, dt_list, label_list, Case_list = \
                    Get_wk_diff_parameters([Nk_m], note, Data_pm, mu_List, RG_info, 5, method_List)
                plt_TargetErr_RunT(wk_list, dt_list, label_list, tit, [ax[py,px], lbx, lby])


        #plt.xlim([-20, 4100])

        FigName = './fig/Fig_{}_Nk_m.pdf'.format(StateName)
        plt.savefig(FigName)
        plt.plot()

    elif Show_Case == 5:                # different [Nk, m]

        #Nk_m_List = [[4, 100], [6, 819]]           #  [Nk, m]
        #Nk_m_List = [[6, 1638], [6, 819]]          #  [Nk, m]
        #Nk_m_List = [[8, 3276], [10, 52428], [12, 838860]]     #  [Nk, m]
        #Nk_m_List = [[6, 819], [8, 3276], [10, 52428], [12, 838860]]     #  [Nk, m]
        #Nk_m_List = [[6, 1228], [8, 13107], [10, 209715], [12, 838860]]     #  [Nk, m] 
        Nk_m_List = [[6, 1228]]     #  [Nk, m] 

        mu_List   = [0.75]
        method_List = [['RGD', 1], ['RGD', 0], ['MiFGD', 0]]               # [method, InitX]
        #method_List = []

        note      = '{}'.format(StateName)

        #wk_list, dt_list, label_list =  Get_Wk_diff_Nk_m(Data_pm, parm_List, InitX, Md_tr, Md_alp, Md_sig, mu_List)       #  old direct version
        #wk_list, dt_list, label_list = Get_wk_diff_parameters(Nk_m_List, note, Data_pm, mu_List, RG_info, Option=5)
        wk_list, dt_list, label_list, Case_list = \
            Get_wk_diff_parameters(Nk_m_List, note, Data_pm, mu_List, RG_info, 5, method_List)

        tit = 'update [Nk, m]'
        plt_TargetErr_RunT(wk_list, dt_list, label_list, tit)

        #plt.xlim([-20, 4100])
        plt.plot()

    elif Show_Case == 6:                # different [Nk, m, version]

        mu_List   = [0.75]
        method_List = [['RGD', 1], ['RGD', 0], ['MiFGD', 0]]               # [method, InitX]
        #method_List = []

        #parm_List = [[6, 409, [2, 1, -1]], [6, 819, [1, 1, -1]], [6, 1638, [1, 1, -1]]]    #  [Nk, m, version]
        #parm_List = [[8, 3276, [1,1,-1]], [10, 52428,[1,1,-1]], [12, 838860,[3,1,-1]]]      #  [Nk, m, version]
        parm_List = [[8, 3276, [1,1,-1]], [10, 104857,[1,1,-1]], [12, 838860,[3,1,-1]]]      #  [Nk, m, version]

        tit       = 'update [Nk, m, version]'    

        note      = '{}'.format(StateName)

        #wk_list, dt_list, label_list = Get_wk_diff_parameters(parm_List, note, Data_pm, mu_List, RG_info, Option=6)
        wk_list, dt_list, label_list = Get_wk_diff_parameters(parm_List, note, Data_pm, mu_List, RG_info, 6, method_List)

        plt_TargetErr_RunT(wk_list, dt_list, label_list, tit)

        #plt.legend(loc='lower right')
        plt.legend()
        plt.xlim([-20, 2500])
        plt.plot()


    elif Show_Case == 10:         # different [StateName, m, version]

        mu_List   = [0.75]
        method_List = [['RGD', 1], ['RGD', 0], ['MiFGD', 0]]               # [method, InitX]
        #method_List = [['RGD', 1], ['RGD', 0]]               # [method, InitX]
        #method_List = [['RGD', 1], ['MiFGD', 0]]               # [method, InitX]
        #method_List = []

        case_Nk = 10
        if case_Nk == 6:
            #parm_List = [['GHZ', 819, [1,1,-1],0], ['Had', 1228,[1,1,-1],0], ['rand', 819, [1,1,0],[1,0]]]      #  [StateName, m, version]
            #parm_List = [['GHZ', 819, [1,1,-1],0], ['Had', 1228,[1,1,-1],0]]      #  [StateName, m, version]

            #parm_List = [['GHZ', 819, [1,1,-1],0]]      #  [StateName, m, version]
            #parm_List = [['GHZ', 1228, [1,1,-1],0]]      #  [StateName, m, version]
            #parm_List = [['Had', 1228,[1,1,-1],0]]      #  [StateName, m, version]
            #parm_List = [['rand', 819, [1,1,0],[1,0]]]      #  [StateName, m, version]
            parm_List = [['rand', 1228, [1,1,0],[1,0]]]      #  [StateName, m, version]

        elif case_Nk == 8:
            #parm_List = [['GHZ', 6553, [1,1,-1],0], ['Had', 6553,[1,1,-1],0], ['rand', 6553,[2,1,0], [1,0]]]      #  [StateName, m, version]
            #parm_List = [['GHZ', 6553, [1,1,-1],0], ['Had', 6553,[1,1,-1],0]]      #  [StateName, m, version]

            #parm_List = [['GHZ', 6553, [1,1,-1], 0]]      #  [StateName, m, version]
            #parm_List = [['GHZ', 13107, [1,1,-1], 0]]      #  [StateName, m, version]

            #parm_List = [['Had', 6553,[1,1,-1], 0]]      #  [StateName, m, version]
            #parm_List = [['Had', 13107,[1,1,-1], 0]]      #  [StateName, m, version]
            
            #parm_List = [['rand', 6553,[2,1,0], [1,0]]]      #  [StateName, m, version]
            parm_List = [['rand', 13107,[1,1,0], [1,0]]]      #  [StateName, m, version]

        elif case_Nk == 10:
            #parm_List = [['GHZ', 104857, [1,1,-1],0], ['Had', 104857,[1,1,-1],0], ['rand', 52428,[1,1,0], [1,0]]]      #  [StateName, m, version]

            #parm_List = [['GHZ', 104857, [1,1,-1],0]]      #  [StateName, m, version]
            #parm_List = [['Had', 104857,[1,1,-1],0]]      #  [StateName, m, version]
            #parm_List = [['rand', 52428,[1,1,0], [1,0]]]      #  [StateName, m, version]
            parm_List = [['rand', 104857,[1,1,0], [1,0]]]      #  [StateName, m, version]

        elif case_Nk == 12:
            #parm_List = [['GHZ', 838860, [3,1,-1],0], ['Had', 838860,[3,1,-1],0], ['rand', 838860,[1,1,0], [1,0]]]      #  [StateName, m, version]

            #parm_List = [['GHZ', 838860, [3,1,-1],0]]      #  [StateName, m, version]
            #parm_List = [['Had', 838860,[3,1,-1],0]]      #  [StateName, m, version]
            parm_List = [['rand', 838860,[1,1,0], [1,0]]]      #  [StateName, m, version]

        tit       = 'update [StateName, m, version]'    

        note      = '{}'.format(Nk)

        #wk_list, dt_list, label_list = Get_wk_diff_parameters(parm_List, note, Data_pm, mu_List, RG_info, Option=6)
        wk_list, dt_list, label_list = Get_wk_diff_parameters(parm_List, note, Data_pm, mu_List, RG_info, 10, method_List)

        plt_TargetErr_RunT(wk_list, dt_list, label_list, tit)

        plt.legend()
        #plt.xlim([-20, 2500])
        plt.plot()

    elif Show_Case == 71:      # different [RGD_info, InitX]   &   diff [r]

        mu_List   = [0.75]
        method_List = [['RGD', 1], ['RGD', 0], ['MiFGD', 0]]               # [method, InitX]
        #method_List = []


        note    = '{}'.format(StateName)
        #tit0    = '{}: $k$ = {}'.format(StateName, Nk)    
        tit0    = '{}({})'.format(StateName, Nk)    

        fig, ax = plt.subplots(2,2, figsize=(8,6))
        figCap = ['(a)', '(b)', '(c)', '(d)']

        method = 2
        if method == 2:
            Nr_List = [3, 6, 8, 10]

            for ii, Nr in enumerate(Nr_List):
                px = ii%2
                py = int(ii/2)
                lbx = py
                lby = (px+1)%2      # to label y or not
                print('py = {}, px = {}, lby = {}, lbx = {}'.format(py, px, lby, lbx))
                tit = '{} {}, $r$ = {}'.format(figCap[ii], tit0, Nr)
                wk_list, dt_list, label_list, Case_list = \
                    Get_wk_diff_parameters([Nr], note, Data_pm, mu_List, RG_info, 7, method_List)
                plt_TargetErr_RunT(wk_list, dt_list, label_list, tit, [ax[py,px], lbx, lby])

        #plt.legend()
        #plt.xlim([-20, 2500])

        FigName = './fig/Fig_rand{}_diff_Nr.pdf'.format(Nk)
        plt.savefig(FigName)
        #plt.plot()

    elif Show_Case == 7:      # different [RGD_info, InitX]   &   diff [r]

        mu_List   = [0.75]
        method_List = [['RGD', 1], ['RGD', 0], ['MiFGD', 0]]               # [method, InitX]
        #method_List = []

        #Nr_List = [3, 6, 8, 10]
        Nr_List = [10]

        note      = '{}'.format(StateName)

        wk_list, dt_list, label_list = Get_wk_diff_parameters(Nr_List, note, Data_pm, mu_List, RG_info, 7, method_List)

        tit       = '{}: Nk = {}'.format(StateName, Nk)    
        plt_TargetErr_RunT(wk_list, dt_list, label_list, tit)

        #plt.legend()
        #plt.xlim([-20, 2500])
        plt.plot()

    elif Show_Case == 8:      # different [RGD_info, InitX]   &   diff [Nk]

        Nk_List = [10]
        mu_List   = [0.75]
        method_List = [['RGD', 1], ['RGD', 0], ['MiFGD', 0]]               # [method, InitX]
        #method_List = []

        note      = '{}'.format(StateName)

        wk_list, dt_list, label_list = Get_wk_diff_parameters(Nk_List, note, Data_pm, mu_List, RG_info, 8, method_List)

        tit       = '{}'.format(StateName)    
        plt_TargetErr_RunT(wk_list, dt_list, label_list, tit)


    elif Show_Case == 9:                            # different [RGD_info, InitX]

        Null_List = [0]
        mu_List   = [0.75, 0.5]
        method_List = [['RGD', 1], ['RGD', 0], ['MiFGD', 0]]               # [method, InitX]

        note      = '{}'.format(StateName)

        wk_list, dt_list, label_list = Get_wk_diff_parameters(Null_List, note, Data_pm, mu_List, RG_info, 1, method_List)

        tit       = '{}: Nk = {}'.format(StateName, Nk)    
        plt_TargetErr_RunT(wk_list, dt_list, label_list, tit)

    #  to implement
    #   updatet  mea (shots)



def template():

    Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method = Data_pm
 
    wk_list = []
    dt_list = []
    pm_list = []
    label_list = []

    Data_pm = [Nk, StateName, m, mea, Nr, version, New_Pj_shot, StVer, Pj_method, mea_method, measure_method]
    #print('Data_pm = {}'.format(Data_pm))

    Setup = Data_Info(Data_pm)

    # -----------------     RGD     --------------------- #
    pm_RGD = 'RGD', InitX, Md_tr, Md_alp, Md_sig
    wk_list, dt_list, pm_list = Get_wk_list(Setup, version, wk_list, dt_list, pm_list, pm_RGD)




# %%

