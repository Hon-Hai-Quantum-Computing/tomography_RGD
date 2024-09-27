
# coding: utf-8

import random
import numpy as np


import qiskit
#from qiskit.tools import outer
from qiskit.quantum_info import Statevector
from qiskit.providers.basic_provider import BasicSimulator

#import sys
#sys.path.append('../')
#from Utility import Generate_All_labels

class State:
    """ qiskit circuit generator of specified quantum state
    """
    def __init__(self, n):
        """Initializes State class

            - n: number of qubits
            - quantum register: object to hold quantum information
            - classical register: object to hold classical information
            - circuit_name: circuit name; defined in each subclass (GHZState, HadamardState, RandomState)

        Args:
            n (int): number of qubits
        """
        #def A(S, rho):
        #    # S random set
        #    # rho density op
        #    return np.sqrt(2 ** d) / np.sqrt(m) * np.asarray([qu.expect(op, rho) for op in S])

        quantum_register   = qiskit.QuantumRegister(n, name='qr')
        classical_register = qiskit.ClassicalRegister(n, name='cr')

        self.quantum_register   = quantum_register
        self.classical_register = classical_register
    
        self.n = n
        
        self.circuit_name = None
        self.circuit      = None
        self.measurement_circuit_names = []
        self.measurement_circuits      = []

        
    def create_circuit(self):
        raise NotImplemented

    
    def execute_circuit(self):
        # XXX not needed?
        pass

    
    def get_state_vector(self):
        """Executes circuit by connecting to Qiskit object, and obtain state vector

        Returns:
            (class Statevector): an ndarray specifying the state vector after the circuit has been created
        """
        if self.circuit is None:
            self.create_circuit()
        
        # XXX add probe?
        # ****** for old qiskit 0.x version, e.g. 0.13.0  ******
        #backend_engine = qiskit.Aer.get_backend('statevector_simulator')
        #job     = qiskit.execute(self.circuit, backend_engine)
        #result  = job.result()
        #state_vector = result.get_statevector(self.circuit_name)

        # *****  usage of new qiskit 1.x version  *****
        state_vector = Statevector(self.circuit)

        return state_vector        

    
    def get_state_matrix(self):
        """Obtain density matrix by taking an outer product of state vector

        Returns:
            (ndarray of complex128): the target denstiy matrix
        """
        state_vector = self.get_state_vector()
        #state_matrix = outer(state_vector)         # deprecated in qiskit 1.x

        state_conj = np.array(state_vector).conj()
        state_matrix = np.outer(state_vector, state_conj)
        return state_matrix
    
    
    def create_measurement_circuits(self, labels, label_format='big_endian'):
        """Prepares measurement circuits

        Args:
            labels (list of str): each element is a string of Pauli matrices (e.g. XYZXX) 
            label_format (str, optional): Specify the ordering of qubits. Defaults to 'big_endian'.
        """
        
        if self.circuit is None:
            self.create_circuit()
        
        qubits = range(self.n)
        
        for label in labels:

            # for aligning to the natural little_endian way of iterating through bits below
            if label_format == 'big_endian':
                effective_label = label[::-1]
            else:
                effective_label = label[:]
            probe_circuit = qiskit.QuantumCircuit(self.quantum_register, 
                                           self.classical_register, 
                                           name=label)

            for qubit, letter in zip(*[qubits, effective_label]): 
                probe_circuit.barrier(self.quantum_register[qubit])
                if letter == 'X':
                    #probe_circuit.u2(0.,       np.pi, self.quantum_register[qubit])     # H (qiskit 0.x)
                    probe_circuit.u(0.5*np.pi, 0., np.pi, self.quantum_register[qubit])  # H (qiskit 1.x)

                elif letter == 'Y':
                    #probe_circuit.u2(0., 0.5 * np.pi, self.quantum_register[qubit])           # H.S^* (qiskit 0.x)
                    probe_circuit.u(0.5*np.pi, 0., 0.5 * np.pi, self.quantum_register[qubit])  # H.S^* (qiskit 1.x)

                probe_circuit.measure(self.quantum_register[qubit], 
                                      self.classical_register[qubit])
                
            measurement_circuit_name = self.make_measurement_circuit_name(self.circuit_name, label)    
            #measurement_circuit      = self.circuit + probe_circuit             #  qiskit 0.x
            measurement_circuit      = self.circuit.compose(probe_circuit)       #  qiskit 1.x
            
            self.measurement_circuit_names.append(measurement_circuit_name)
            self.measurement_circuits.append(measurement_circuit)
        

    @staticmethod    
    def make_measurement_circuit_name(circuit_name, label):
        """ Measurement circuit naming convention

        Args:
            circuit_name (str): name of the circuit, usually the state name
            label (str): the Pauli string referring to the measurement

        Returns:
            (str): name of measurement circuit
        """
        name = '%s-%s' % (circuit_name, label)
        return name

    
    def execute_measurement_circuits(self, labels:list[str],
                                     backend   = BasicSimulator(),
                                     num_shots = 100,
                                     #num_shots = 10000,
                                     label_format='big_endian'):
        """ Executes measurement circuits

        Args:
            labels (list[str]):list of labels where each label is a string of Pauli matrices (e.g. XYZXX) 
            backend (class BasicSimulator, optional): a class defining a quantum simulator or a real device 
                            Defaults to BasicSimulator().
            num_shots (int, optional): number of shots measurement is taken to get empirical frequency through counts
                        Defaults to 100.
            label_format (str, optional): _description_. Defaults to 'big_endian'.

        Returns:
            (list): each elment is a dictionary storing all information of each Pauli meausrement 
        """
        if self.measurement_circuit_names == []:
            self.create_measurement_circuits(labels, label_format)
        
        circuit_names = self.measurement_circuit_names

        #  ****  usage of old qiskit 0.x  version  ****
        #backend   = 'qasm_simulator'
        #backend_engine = qiskit.Aer.get_backend(backend)
        #job = qiskit.execute(self.measurement_circuits, backend_engine, shots=num_shots)

        #  **** usage of qiskit 1.x version ******
        job = backend.run(self.measurement_circuits, shots=num_shots)

        result = job.result()

        data_dict_list = []
        for i, label in enumerate(labels):
            measurement_circuit_name = self.make_measurement_circuit_name(self.circuit_name, label)
            data_dict = {'measurement_circuit_name' : measurement_circuit_name,
                         'circuit_name'             : self.circuit_name,
                         'label'                    : label,
                         'count_dict'               : result.get_counts(i),
                         'backend'                  : backend,
                         'num_shots'                : num_shots}
            data_dict_list.append(data_dict)
        return data_dict_list


class GHZState(State):
    """Constructor for GHZState class

    Args:
        State (class State): the constructed circuit
    """
    def __init__(self, n):
        """ initialization of GHZState class by calling class State

        Args:
            n (int): number of quibts of the state
        """
        State.__init__(self, n)
        self.circuit_name = 'GHZ'

        
    def create_circuit(self):    
        """
        Creates a quantum circuit for a GHZ state using qiskit
        """    
        circuit = qiskit.QuantumCircuit(self.quantum_register, 
                                 self.classical_register, 
                                 name=self.circuit_name)

        circuit.h(self.quantum_register[0])
        for i in range(1, self.n):
            circuit.cx(self.quantum_register[0], self.quantum_register[i])   
        
        self.circuit = circuit

        
class HadamardState(State):
    '''
    Constructor for HadamardState class
    '''
    def __init__(self, n):        
        """ initialization of Hadamard State class by calling class State

        Args:
            n (int): number of quibts of the state
        """
        State.__init__(self, n)
        self.circuit_name = 'Hadamard'
        
        
    def create_circuit(self):
        """
        Creates a quantum circuit for a Hadamard state using qiskit
        """
        circuit = qiskit.QuantumCircuit(self.quantum_register, 
                                 self.classical_register, 
                                 name=self.circuit_name)
        
        for i in range(self.n):
            circuit.h(self.quantum_register[i])   
        
        self.circuit = circuit

        
class RandomState(State):
    '''
    Constructor for RandomState class
    '''
    def __init__(self, n, seed=0, depth=40):
        State.__init__(self, n)
        self.circuit_name = 'Random-%d' % (self.n, )

        self.seed  = seed
        self.depth = depth

        
    def create_circuit(self):
        random.seed(a=self.seed)
        circuit = qiskit.QuantumCircuit(self.quantum_register, 
                                 self.classical_register, 
                                 name=self.circuit_name)

        for j in range(self.depth):
            if self.n == 1:
                op_ind = 0
            else:
                op_ind = random.randint(0, 1)
            if op_ind == 0: # U3
                qind = random.randint(0, self.n - 1)
                #circuit.u3(random.random(), random.random(), random.random(),
                #           self.quantum_register[qind])                            # qiskit 0.x
                circuit.u(random.random(), random.random(), random.random(),
                           self.quantum_register[qind])                             # qiskit 1.x

            elif op_ind == 1: # CX
                source, target = random.sample(range(self.n), 2)
                circuit.cx(self.quantum_register[source],
                           self.quantum_register[target])
        
        self.circuit = circuit


if __name__ == '__main__':

    ############################################################
    ### Example of creating and running an experiment
    ############################################################

    n = 3
    #labels = projectors.generate_random_label_list(20, n)
    labels = ['YXY', 'IXX', 'ZYI', 'XXX', 'YZZ']
    #labels = ['YZYX', 'ZZIX', 'XXIZ', 'XZIY', 'YXYI', 'ZYYX', 'YXXX', 'IIYY', 'ZIXZ', 'IXXI', 'YZXI', 'ZZYI', 'YZXY', 'XYZI', 'XZXI', 'XZYX', 'YIXI', 'IZYY', 'ZIZX', 'YXXY']
    #labels = ['IIIX', 'IYIY', 'YYXI', 'ZZYY', 'ZYIX', 'XIII', 'XXZI', 'YXZI', 'IZXX', 'YYIZ', 'XXIY', 'XXZY', 'ZZIY', 'YIYX', 'YYZZ', 'YZXZ', 'YZYZ', 'ZXYY', 'IXIZ', 'XZII']
    #labels = Generate_All_labels(n)

    #state   = GHZState(n)
    #state   = HadamardState(n)
    state   = RandomState(n)

    #
    # DO the shot measurement
    #
    state.create_circuit() 
    data_dict_list = state.execute_measurement_circuits(labels)


    target_density_matrix = state.get_state_matrix()
    target_state          = state.get_state_vector()            
