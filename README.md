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

   - Nk    =  number of qubits
   - m     =  number of sampled Pauli matrices
   - mea   =  number of shot measurements
```
         [usage]   either directly define them
                  or     call  basic_sys_setting(Choose) @  Input_setting.py
```

   - StateName  =  'GHZ',  'Had',  'rand'   (where 'GHZ' and 'Had' are pure states)

   - Nr  =   number of rank  [  fixed to 1   for  'GHZ', 'Had'   ]
   
   - Noise   =   noise model  [only needed for mixed sates]
      - only implementing exact results now, i.e. not really implemented in the code
      -  different types for different states, i.e.
         -  pure states: = 0            (GHZ or Had)
         -  'rand':      = [0, 0, 0.1]  (actually no effect now) 

2.  controlling parameters for each cases, i.e.  labelling each running case
   - Notes for loading or generating sampled Pauli matrices & measurements 
      - projectors = labels = Pauli matrices (obs = observables) 
      - measurement: 
         - using qiskit package: shot measurement for pure states       
         - using qutip package: direct exact calculation for mixed states  
         - Note the measurements for pure states & mixed states are different for the treatment now                           

	- StVer   =  control of the 'rand' state generation   
      - 'pure states': 0, i.e.  no need of this parameter
      - 'rand': list of [which version of generated data, to generate or not], i.e.
         - [1, 0] = [version 1, do not generate] 
         - [1, 1] = [version 1, generate a new one]         

	- version =  [version of sampled Pauli operators, version of generated measurement,  Noise]
      - version[0] =  version of sampled Projectors             
      - version[1] =  version of measurement                          
      - version[2] =  specification of Noise                                
      - each version will have separate directory
         - if loading Proj/measurement  --> will check if directory exists 
         - if generating new Pj/measure --> will create new directory     

3.  program control parameters for projectors | measurement

   - Pj_method      = which method to generate | load projectors   (default = 1)
      - = 0  (each Pauli matrix is stored separately)          
      - = 1  (all Pauli matrices saved in chunks)    
         - if: not too many --> saved in one Pj_List.pickle file 
         - else: will save in several chunks of Pj_list_xx.pickle 
 
   - mea_method     = which method to generate | load measurements (default = 1)
      - = 0  (measurement for each Pauli obs is saved separately)  
      - = 1  (measurement result is stored in chunks)              

	- measure_method = to use parallel generation of measurement 
      - = 1   directly calculate the whole label_list, i.e. no parallel calculation
      - = 3   parallel measurement for each label part, i.e. parallel calculation [default]

4.  Data_In = [DirRho, Dir_proj, Dir_meas]
   - Note
      - (default)  = []   (i.e. not specified)
      - can also be [DirRho, Dir_proj]  (not specified Dir_meas)
      - each directory can be empty '' (must check the consistency)
   - DirRho = directory of storing Rho (only necessary for 'rand')
   - Dir_proj = directory of storing projectors
   - Dir_meas = directory of storing measurement results

5.  New_Pj_shot = [how to obtain sampled Pauli operator, how to have measurements]
   - default usage:
      - New_Pj_shot = [1, 1] : creating new Projector & measurement  
      - Caveat: other choice of New_Pj_shot must be consistent with Data_In, existing directory specified (Directory of Rho, projectors, measurements), etc.

   - New_Pj_shot[0] = specify how to obtain the sampled Pauli operators
      - = 0  loading sampled projectors              
      - = 1  creating new sampled projectors    
            (if the version already exists, then the version New_Pj_shot[0] will +1 automatically)
      - = -1  Pauli operators generated from single & nearest-neighbor site Pauli matrices 
      - =  2  combining each Pj_list_xx.pickle into single file first 
    
   - New_Pj_shot[1] = specify how to obtain the measured results
      - = 0  loading measurement result from file                   
      - = 1  doing measurement (qiskit shot or qutip calculation)    
              (if the version already exists, then the version New_Pj_shot[1] will +1 automatically)
      - = 2  converting shot measured data_dict into the coef needed, 
            where measurement_list = the coef of each Pauli obs
            (no longer maintained -> may need check)


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
