

#import numpy as np
#from numpy import linalg as LA
#from functools import reduce

#from mpi4py import MPI
#from qiskit.quantum_info import state_fidelity

import projectors
import measurements

import os
import pickle
import time

# --------------------------------- #
#	the origianl RGD package		#
# --------------------------------- #
#from projectors import Sm_dxd_list
#import qutip as qu

#import scipy.sparse as sparse
#import time

############################################################
## Parallel projFGD base class
## XXX WIP
############################################################
class WorkerParm:
	def __init__(self,
				 process_idx,
				 num_processes,
				 params_dict, input_S):
#				 params_dict):

		projector_store_path     = params_dict.get('projector_store_path', None)
		num_iterations           = params_dict['num_iterations']

		beta                     = params_dict.get('beta', None)
		trace                    = params_dict.get('trace', 1.0)
		target_state             = params_dict.get('target_state', None)
		target_DM                = params_dict.get('target_DM', None)

		convergence_check_period = params_dict.get('convergence_check_period', 10)
		relative_error_tolerance = params_dict.get('relative_error_tolerance', 0.0001)		#  original setting
		#relative_error_tolerance = params_dict.get('relative_error_tolerance', 0.001)

		parity_flavor            = params_dict.get('parity_flavor', 'effective')

		pauli_correlation_measurements_fpath = params_dict.get('pauli_correlation_measurements_fpath',
															   None)
		pauli_basis_measurements_fpath       = params_dict.get('pauli_basis_measurements_fpath',
															   None)

		# alternatives: 'pauli_basis', 'pauli_correlation'
		measurements_type                    = params_dict.get('measurements_type', 'pauli_correlation') 

		measurement_store_path   = params_dict.get('measurement_store_path', None)

		tomography_labels            = params_dict.get('tomography_labels', None)
		tomography_bitvectors_list   = params_dict.get('tomography_bitvectors_list', None)
		
		density_matrix           = params_dict.get('density_matrix', None)

		label_format             = params_dict.get('label_format', 'big_endian')

		# XXX deactivated
		# store_load_batch_size    = params_dict.get('store_load_batch_size', 10)
		debug                    = params_dict.get('debug', True)  

		seed                     = params_dict.get('seed', 0)

		Nr                       = params_dict.get('Nr', 1)		#  RGD: rank of the target matrix
		self.Nr                  = Nr

		self.StateName           = params_dict.get('StateName', None)
		self.Noise               = params_dict.get('Noise', None)

		Pj_method                = params_dict.get('Pj_method', 0)
		mea_method               = params_dict.get('mea_method', 0)
		measure_method           = params_dict.get('measure_method', 1)
		version                  = params_dict.get('version', None)

		self.Pj_method           = Pj_method
		self.meas_method         = mea_method
		
		self.meas_path           = measurement_store_path

		#print('  tomography_labels = {}'.format(tomography_labels))
		print('  self.Noise = {}'.format(self.Noise))
		print('  Pj_method  = {}'.format(Pj_method))
		print('  mea_method = {}'.format(mea_method))

		#################################################################################
		#		loading label_list from pickle											#
		#		(eg)	tomography_labels   = ['XII', 'YIY', 'YZY', 'ZXY', 'ZXZ']		#
		#				label_list          = ['XII', 'YIY', 'YZY', 'ZXY', 'ZXZ']		#
		#################################################################################

		if self.Noise != None:		
			if Pj_method == 0:			#  original method
				tomography_labels = [fname.split('.')[0] for fname in os.listdir(projector_store_path)]
			elif Pj_method == 1:		#  new method
				tomography_labels = projectors.ProjectorStore.load_labels_from_file(projector_store_path)
			#print('  -->  loading from Projectors: tomography_labels = {}'.format(tomography_labels))

		# tomography_labels (and tomography_bitvectors)
		if tomography_labels == None:
			if mea_method == 0:			#  the original method
				tomography_labels  = measurements.MeasurementStore.load_labels(measurement_store_path)
				
			elif mea_method == 1:		#  the new method
				if measure_method == 1:			
					tomography_labels  = measurements.MeasurementStore.load_labels_from_file(measurement_store_path)
				elif measure_method == 3:
					measurement_output = measurements.MeasurementStore.Load_measurement_list(measurement_store_path, 'ALL')
					tomography_labels = measurement_output[0]

					#print(' len(measure_output)  = {}'.format(len(measurement_output)))
					#print(' measurement_output   = {}'.format(measurement_output))
					#print(' measurement_output[0]   = {}'.format(measurement_output[0]))
					#print(' measurement_output[1]   = {}'.format(measurement_output[1]))

		#print(' tomography_labels   = {}'.format(tomography_labels))

		n = len(tomography_labels[0])
		d = 2 ** n

		num_tomography_labels = len(tomography_labels) 
		label_list            = split_list(tomography_labels, num_processes)[process_idx]
		num_labels            = len(label_list)

		#print(' tomography_labels   = {}'.format(tomography_labels))
		#print(' num_processes       = {}'.format(num_processes))
		#print(' label_list          = {}'.format(label_list))
		#print(' split = {}'.format(split_list(tomography_labels, 2)))
		#print(' process_idx         = {}'.format(process_idx))
		#print(' measurements_type   = {}'.format(measurements_type))	## here = pauli_correlation
		#self.label_list = label_list

		################################################################################
		# Loading projectors
		########Loading projectors########################################################################

		#print('          projector_store_path = {}'.format(projector_store_path))
		#print('          num_labels           = {}'.format(num_labels))
		#print('          debug                = {}'.format(debug))
		tp1 = time.time()

		if projector_store_path is not None:
			# load projectors in batches
			projector_dict = {}
			start = 0
			# end   = min(num_labels,  store_load_batch_size)
			end = num_labels
			if debug:
				print('Loading projectors')                
			while True:
				if measurements_type == 'pauli_correlation':
					if Pj_method == 0:		#  original method: saving each Pj separately
						data_dict = projectors.ProjectorStore.load(projector_store_path, label_list[start:end])
					elif Pj_method == 1:	#  new method: saving all together
						#data_dict = projectors.ProjectorStore.load_PoolMap(projector_store_path)
						data_dict = projectors.ProjectorStore.load_PoolMap(projector_store_path, label_list[start:end])

				elif measurements_type == 'pauli_basis':
					data_dict = projectors.PauliStringProjectorStore.load(projector_store_path,
																		  label_list[start:end],
																		  bitvector_list[start:end])
				projector_dict.update(data_dict)
				start = end
				# end   = min(num_labels, start  + store_load_batch_size)
				end = num_labels

				#print(' -------------------------------------------------------- ')
				#print('projector_dict  = {}'.format(projector_dict))
				#print('data_dict       = {}'.format(data_dict))
				#print('label_list      = {}'.format(label_list))
				#print('label_list      = {}'.format(label_list[0:num_labels]))
				#print(' start   =  {}'.format(start))
				#print(' end     =  {}'.format(end))
				#print(' debug   =  {}'.format(debug))

				if debug:
					print('   %d projectors loaded' % len(projector_dict))
				if start == num_labels:
					break

				#print(' ----------------- (not entering here)  -------------------------- ')
				#print('projector_dict  = {}'.format(projector_dict))
				#print(' start   =  {}'.format(start))
				#print(' end     =  {}'.format(end))

			if measurements_type == 'pauli_correlation':
				projector_list = [projector_dict[label] for label in label_list]
			elif measurements_type == 'pauli_basis':
				projector_list = []
				for label, bitvector in zip(*[label_list, bitvector_list]):
					projector_list.append(projector_dict[label][bitvector])
				#projector_list = [projector_dict[label][bitvector] for (label, bitvector) in
				#zip(*[label_list, bitvector_list])]
			#print('projector_list = {}'.format(projector_list))
			del projector_dict
		else: 
			if debug:
				print('Computing projectors')
			projector_list   = [projectors.build_projector_naive(label, label_format) for label in label_list]
			projector_list   = [projectors.build_projection_vector(label, bitvector) for label, bitvector
								in zip(*[label_list, bitvector_list])]
		if debug:
			print('   Projectors ready to compute with\n')

		#print(' ********************************************** ')
		#print('  projector_list   = {}'.format(projector_list))
		#print('projector_dict  = {}'.format(projector_dict))	# UnboundLocalError: local variable 'projector_dict' referenced before assignment
		#self.projector_list = projector_list

		tp2 = time.time()
		print('\n  ******  projectors  -->  time = {}   ****** \n'.format(tp2 - tp1))

		################################################################################
		# Loading measurements
		################################################################################

		#print('pauli_correlation_measurements_fpath = {}'.format(pauli_correlation_measurements_fpath))
		#print('pauli_basis_measurements_fpath       = {}'.format(pauli_basis_measurements_fpath))
		#print('measurement_store_path               = {}'.format(measurement_store_path))		
		#print('density_matrix                       = {}'.format(density_matrix))

		tm1 = time.time()

		if pauli_correlation_measurements_fpath is not None:			#  not entering
			if debug:
				print('Loading Pauli correlation measurements')                
			if pauli_correlation_measurements_fpath.endswith('.json'):
				with open(pauli_correlation_measurements_fpath, 'r') as f:
					data_dict = json.load(f)
			measurement_list = [data_dict[label] for label in label_list]
			del data_dict

		if pauli_basis_measurements_fpath is not None:					#  not entering
			if debug:
				print('Loading Pauli basis measurements')                
			if pauli_basis_measurements_fpath.endswith('.json'):
				with open(pauli_basis_measurements_fpath, 'r') as f:
					data_dict = json.load(f)
			measurement_list = [data_dict[label][bitvector] for label, bitvector
								in zip(*[label_list, bitvector_list])]
			del data_dict

		elif self.Noise != None:							#  entering noise model
			print(' *******  using noise model as measured results ********')
			print('measurement_store_path = {}'.format(measurement_store_path))
			
			Model    = self.Noise[0]
			ver_meas = version[1]
			F_noiseMea = '{}/zN{}_v{}_measurements'.format(measurement_store_path, Model, ver_meas)

			#with open(measurement_store_path, 'rb') as f:
			with open(F_noiseMea, 'rb') as f:
				#yProj_Exact, zN, measurement_list, Noise, s_label, params_setup,\
				#	  target_density_matrix, rho = pickle.load(f)
				s_label, measurement_list, zN = pickle.load(f)	

			print('  --->  loading measured results from noise model DONE  \n')
			#print('s_label = {}'.format(s_label))
			#print('yProj   = {}'.format(measurement_list))
			#print('zN      = {}'.format(zN))

		elif measurement_store_path is not None:			#  entering this (normal shot measurements)
			print('   ===>  to GET measurment_list  \n')

			# load measurements in batches
			measurement_dict = {}
			start = 0
			# end   = min(num_labels,  store_load_batch_size)
			end = num_labels
			if debug:
				print('Loading measurements')    

			# ----------------------------------------- #
			# 	to get    measurement_list				#	
			# ----------------------------------------- #
			
			if mea_method == 1:				#  the new method (already calc measurement_list)

				if measure_method == 1:			#   from direct qiskit measurement						
					#measurement_list = measurements.MeasurementStore.Load_measurement_list(measurement_store_path)
					#labs, measurement_list = measurements.MeasurementStore.Load_measurement_list(measurement_store_path)
					output_list = measurements.MeasurementStore.Load_measurement_list(measurement_store_path)

					if len(output_list) == 1:
						measurement_list = output_list
					elif len(output_list) == 2:					#   [labels,  measurement_list]
						measurement_list = output_list[1]
				elif measure_method == 3:		#   from parallel cpus to have measurement_list
					measurement_list = measurement_output[1]

				#print('  -->  loaded measurment_list = \n{}'.format(measurement_list))

			elif mea_method == 0:			#  the original method

				while True:
					if mea_method == 0:			#  the original method: each Pauli op in separate file
						data_dict = measurements.MeasurementStore.load(measurement_store_path, label_list[start:end])
					elif mea_method == 1:		#  the    new   method:  ALL Pauli op in ONE file   (will not enter)
						data_dict = measurements.MeasurementStore.load_OneDict(measurement_store_path, label_list[start:end])

					measurement_dict.update(data_dict)
					start = end
					# end   = min(num_labels, start + store_load_batch_size)
					end = num_labels
					if debug:
						print('   %d measurements loaded' % len(measurement_dict))
					if start == num_labels:
						break
				
				count_dict_list = [measurement_dict[label] for label in label_list]
				measurement_object_list = [measurements.Measurement(label, count_dict) for (label, count_dict) in 
										zip(*[label_list, count_dict_list])]

				#print(' ------------  will enter here ------------------ ')
				#print(' start = {}   |   end = {}'.format(start, end))
				#print(' debug            = {}'.format(debug))
				#print(' label_list       = {}'.format(label_list))
				#print(' count_dict_list  = {}'.format(count_dict_list))
				#print(' parity_flavor    = {}'.format(parity_flavor))
				#print(' beta             = {}'.format(beta))
				#print(' measurement_dict = {}\n'.format(measurement_dict))
				#print(' measurement_object_list = {}\n'.format(measurement_object_list))


				if measurements_type == 'pauli_correlation':
					measurement_list = [measurement_object.get_pauli_correlation_measurement(beta, parity_flavor)[label] for
										(label, measurement_object) in zip(*[label_list,
																			measurement_object_list])]

					#print('\n ****  the measurement_type = {} *** \n'.format(measurements_type))
					print(' measurement_list = {}'.format(measurement_list))				

				elif measurements_type == 'pauli_basis':
					measurement_list = [measurement_object.get_pauli_basis_measurement(beta)[label][bitvector] for
										(label, bitvector, measurement_object) in zip(*[label_list,
																						bitvector_list,
																						measurement_object_list])]
				del count_dict_list, measurement_object_list

				#self.measurement_list = measurement_list
				#print('measurement_list  = {}'.format(measurement_list))

		elif density_matrix is not None:
			if debug:
				print('Computing measurements')
			if measurements_type == 'pauli_correlation':
								measurement_list = get_measurements(density_matrix, projector_list)
		if debug:
			print('   Measurements ready to compute with')

		tm2 = time.time()
		print('\n   **** Loading measurements   -->   time = {}  **** \n'.format(tm2 - tm1))



		self.num_iterations = num_iterations    # number of iterations (gradient steps)
		
		self.trace                     = trace
		self.target_state              = target_state
		self.target_DM				   = target_DM			#  target state -->  density matrix

		self.convergence_check_period  = convergence_check_period   # how often to check convergence
		self.relative_error_tolerance  = relative_error_tolerance   # user-decided relative error 
		self.seed                      = seed
		
		self.label_list     = label_list
		# self.bitvector_list = bitvector_list
		self.measurement_list = measurement_list
		#self.zm               = np.ones(num_labels)

		self.projector_list   = projector_list

		self.num_tomography_labels = num_tomography_labels
		self.num_labels            = num_labels

		# self.num_bitvectors        = num_bitvectors
		
		n = len(label_list[0])
		d = 2 ** n
		
		self.n            = n   # number of qubits
		self.num_elements = d
		

		self.process_idx   = process_idx
		self.num_processes = num_processes

		#self.FrobErr_st                 = []	# 	Frobenius norm difference from State
		#self.FrobErr_Xk				 = [] 	#   Frobenius norm difference from matrix Xk
		self.Target_Err_st              = []	# 	Frobenius norm difference from State
		self.Target_Err_Xk				= [] 	#   Frobenius norm difference from matrix Xk

		self.target_error_list          = []
		self.target_relative_error_list = []


		self.fidelity_Xk_list			= []	#   Fidelity between (Xk, target)
		self.Err_relative_Xk            = []

		self.error_list                 = []
		#self.relative_error_list        = []
		self.Err_relative_st       		= []

		
		self.iteration                  = 0
		self.converged                  = False
		self.convergence_iteration      = 0

		self.fidelity_list              = []

		# additional attributes useful for experimentation
		self.circuit_name = params_dict.get('circuit_name', '')
		self.backend      = params_dict.get('backend', 'local_qiskit_simulator')

		if self.Noise == None:
			self.num_shots    = params_dict.get('num_shots', 8192)
		
		self.complete_measurements_percentage = params_dict.get('complete_measurements_percentage', 40)

	#def Nice(self):
	#	print(' Here you go!')


############################################################
## Utility functions
## XXX To modularize/package
############################################################

def split_list(x, num_parts):
	n = len(x)
	size = n // num_parts
	parts = [x[i * size: (i+1) * size] for i in range(num_parts - 1 )]
	parts.append(x[(num_parts - 1) * size:])
	return parts

