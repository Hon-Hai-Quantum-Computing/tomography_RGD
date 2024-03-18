# Overview

1. The program uses the optimization method of Riemannian Gradient Descent (RGD) algorithm to do the tomography problem.
   The implementation follows from the arXiv:2210.04717

2. Compare the optimization method of Momentum-inspired Factorized Gradient Descent (MiFGD).

# Package dependence
1. The package depedence is the same as those from MiFGD.



# Usage

1. must create a directory ./calc and data stored in ./calc

2. python main_QST.py

----
# program structures

## necessary parameters:

the basic parameters (Nk, iterations, ...), version control (for different data or samples), path for saving/loading projectors or measurements, target density matrix, and so on, are stored in the dictionary variable params_dict.

The path processing the data is defined in 
Utility.py/State_Naming/Dir0.
The default value of Dir0 is ./calc.





1. system parameters

   Nk    =  number of qubits
   m     =  number of sampled Pauli matrices
   mea   =  number of shot measurements

         [usage]   either directly define them
                  or     call  basic_sys_setting(Choose) @  Input_setting.py

   StateName  =  'GHZ',  'Had',  'rand'   (where 'GHZ' and 'Had' are pure states)

   Nr      =   number of rank  [  fixed to 1   for  'GHZ', 'Had'   ]
   
   Noise   =   noise model  [only needed for mixed sates]
         [  only implementing exact results now, i.e. not really implemented in the code]
         different types for states
         i.e.  pure states: = 0            (GHZ or Had)
               'rand':      = [0, 0, 0.1]  (actually no effect now)


2.  running cases record, i.e.  labelling the run cases 

	StVer   =  control of the 'rand' state generation   
   ('pure states': 0;  
    'rand': [which version of generated data, to generate or not), 
            (eg)  [1, 0] = [version 1, do not generate] 
                  [1, 1] = [version 1, generate a new one]         
   )

	version =  [version of generated states, version of generated measurement,  Noise]



the optimization program is to run directly
   - MiFGD_optRun.py   for MiFGD method
   - RGD_optRun.py     for  RGD  method

In MiFGD_otpRun.py, execute 
    worker = methodsMiFGD.BasicWorker(params_dict, input_S)
    worker.compute(InitX_MiFGD, Ld_MiFGD)


as long as all the data are prepared in the directories.
The path for the saving/loading the projectos and measured y 





main_QST is the example to do the optimization over prepared measured y from the target_density_matrix. Feel free to do the same process by directly calling each necessary ingredient, including

1. 


## producing sampled Pauli matrices 

Input_Utility/Gen_Projectors

## producing target_density_matrix, measured y, 

Input_measurements/Get_measurement_by_labels
