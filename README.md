# Overview

1. The program uses the optimization method of Riemannian Gradient Descent (RGD) algorithm to do the tomography problem.
   The implementation follows from the arXiv:2210.04717

2. Compare the optimization method of Momentum-inspired Factorized Gradient Descent (MiFGD). 

# Code implementation and Installation
1. The code development is based on the MiFGD code from
https://github.com/gidiko/MiFGD. 

2. The repository can call either the MiFGD or RGD method for calculating the same input data.

3.  Create a conda environment `QST` by the following command:
```
conda env create -f environment.yml 
```


# Usage

1. must create a directory ./calc and data stored in ./calc

2. Run the exampleRun.ipynb in the jupyternote environment.

----
# program structures

## necessary parameters:

The basic parameters (Nk, iterations, ...), version control (for different data or samples), path for saving/loading projectors or measurements, target density matrix, and so on, are stored in the dictionary variable params_dict.

The path processing the data is defined in 
Utility.py/State_Naming/Dir0.
The default value of Dir0 is ./calc.


1. system parameters

   -  Nk    =  number of qubits
   -  m     =  number of sampled Pauli matrices
   -  mea   =  number of shot measurements

         [usage]   either directly define them
                  or     call  basic_sys_setting(Choose) in  Input_setting.py

   - StateName  =  'GHZ',  'Had',  'rand'   (where 'GHZ' and 'Had' are pure states)

   - Nr  =   number of rank  [  fixed to 1   for  'GHZ', 'Had'   ]
   
   - Noise   =   noise model  [only needed for mixed sates]
      - only implementing exact results now, i.e. not really implemented in the code
      -  different types for different states, i.e.
         -  pure states: = 0            (GHZ or Had)
         -  'rand':      = [0, 0, 0.1]  (actually no effect now) 

2.  running cases record, i.e.  labelling the run cases 

	- StVer   =  control of the 'rand' state generation   
      - 'pure states': 0, i.e.  no need of this parameter
      - 'rand': list of [which version of generated data, to generate or not], i.e.
         - [1, 0] = [version 1, do not generate] 
         - [1, 1] = [version 1, generate a new one]         

	- version =  [version of generated states, version of generated measurement,  Noise]


## Calling the optimization method

1. the optimization program is to run directly
   - MiFGD_optRun.py   for MiFGD method
   - RGD_optRun.py     for  RGD  method

2. In MiFGD_otpRun.py, execute 
   - worker = methodsMiFGD.BasicWorker(params_dict, input_S)
   - worker.compute(InitX_MiFGD, Ld_MiFGD)

3. main_QST is the example to do the optimization over prepared measured y from the target_density_matrix. Feel free to do the same process by directly calling each necessary ingredient, 
   - as long as all the data are prepared in the specified directories.
   - The path for the saving/loading the projectors and measured y can be different


## producing sampled Pauli matrices 

Input_Utility/Gen_Projectors

## producing target_density_matrix, measured y, 

Input_measurements/Get_measurement_by_labels
