

import os
import pickle
import multiprocessing
import numpy as np
import time

# ----------------------------------------- #
#   class method to deal with each label    #
# ----------------------------------------- #

class Measurement:
    def __init__(self, label, count_dict):
        '''
        Initializes Measurement class
        - labels: string of Pauli matrices (e.g. XYZXX)
        - count_dict: counts of each possible binary output string (e.g. '000...00', '000...01', etc.)
        - num_shots: number of shots measurement is taken to get empirical frequency through counts
        '''
        akey = list(count_dict.keys())[0]
        assert len(label) == len(akey)

        self.label            = label
        self.count_dict       = self.zfill(count_dict)
        self.num_shots        = sum(count_dict.values())

    @staticmethod
    def zfill(count_dict):
        n = len(list(count_dict.keys())[0])
        d = 2 ** n
        result = {}
        for idx in range(d):
            key = bin(idx)[2:]
            key = key.zfill(n)
            result[key] = count_dict.get(key, 0)
        return result



    @staticmethod
    def naive_parity(key):
        return key.count('1')

    
    def effective_parity(self, key):
        indices    = [i for i, symbol in enumerate(self.label) if symbol == 'I']
        digit_list = list(key)
        for i in indices:
            digit_list[i] = '0'
        effective_key = ''.join(digit_list)

        #print(' *************  parity  ********************')
        #print(' label           = {}'.format(self.label))
        #print(' key             = {}'.format(key))
        #print(' indices         = {}'.format(indices))
        #print('digit_list       = {}'.format(digit_list))
        #print('effective_key    = {}'.format(effective_key))
        #print('effective parity = {}'.format(effective_key.count('1')))

        return effective_key.count('1')

    def parity(self, key, parity_flavor='effective'):
        if parity_flavor == 'effective':
            return self.effective_parity(key)
        else:
            return self.naive_parity(key)
            
    
    def get_pauli_correlation_measurement(self, beta=None, parity_flavor='effective'):
        '''
        Generate Pauli correlation measurement (expectation value of Pauli monomials).
        Note that summation of d=2^n Pauli basis measurement corresponds to one Pauli correlation measurement.
        '''
        if beta == None:
        #    beta = 0.50922         #  original MiFGD usage
            beta = 0.0              #  this one looks more exact
        num_shots          = 1.0 * self.num_shots
        num_items          = len(self.count_dict)

        #frequencies        = {k : (v + beta) / (num_shots + num_items * beta) for k, v in self.count_dict.items()}
        #parity_frequencies = {k : (-1) ** self.parity(k, parity_flavor) * v for k, v in frequencies.items()}
        #correlation        = sum(parity_frequencies.values())
        #data = {self.label : correlation}

        freq2          = {k : (v ) / (num_shots ) for k, v in self.count_dict.items()}
        parity_freq2   = {k : (-1) ** self.parity(k, parity_flavor) * v for k, v in freq2.items()}
        correlation2   = sum(parity_freq2.values())
        data2 = {self.label : correlation2}

        #print(' ------------------------------  ')
        #print('freq2                   = {}'.format(freq2))
        #print('parity_freq2            = {}'.format(parity_freq2))
        #print('parity_freq2.values()   = {}'.format(parity_freq2.values()))
        #print('correlation2            = {}'.format(correlation2))
        #print(' ------------------------------  ')
        #print('parity_flavor           = {}'.format(parity_flavor))
        #print('beta                    = {}'.format(beta))
        #print('label                   = {}'.format(self.label))
        #print('self.count_dict         = {}'.format(self.count_dict))
        #print('self.count_dict.items() = {}'.format(self.count_dict.items()))
        #print('num_shots               = {}'.format(num_shots))
        #print('num_items               = {}'.format(num_items))
        #print('frequencies             = {}'.format(frequencies))
        #print('frequencies.items()     = {}'.format(frequencies.items()))

        #print('parity_frequencies      = {}'.format(parity_frequencies))
        #print('data                    = {}'.format(data))

        return data2
        #return data

    def get_pauli_basis_measurement(self, beta=None):
        '''
        Generate Pauli basis measurement. 
        Note that summation of d=2^n Pauli basis measurement corresponds to one Pauli correlation measurement.
        '''
        if beta == None:
            beta = 0.50922
        num_shots   = 1.0 * self.num_shots
        num_items   = len(self.count_dict)
        frequencies = {k : (v + beta) / (num_shots + num_items * beta) for k, v in self.count_dict.items()}
        data = {self.label: frequencies}
        return data



    def _pickle_save(self, fpath):
        with open(fpath, 'wb') as f:
            pickle.dump(self.count_dict, f)


    def _hdf5_save(self, fpath):
        f = h5py.File(fpath, 'a')
        group = f.create_group(self.label)
        
        items = [[key, value] for key, value in self.count_dict.items()]
        keys   = np.array([item[0] for item in items], dtype='S')
        values = np.array([item[1] for item in items], dtype='int32')
        
        dataset = group.create_dataset('keys',   data = keys)
        dataset = group.create_dataset('values', data = values) 
        f.close()

    # Save the measurement to disk
    def save(self, path):
        if os.path.isdir(path):
            fpath = os.path.join(path, '%s.pickle' % self.label)
            self._pickle_save(fpath)
        elif path.endswith('.hdf5'):
            fpath = path
            self._hdf5_save(fpath)

            
    @classmethod
    def _pickle_load(cls, fpath):
        with open(fpath, 'rb') as f:
            data = pickle.load(f)
        return data


    @classmethod
    def _hdf5_load(cls, fpath, label):

        f = h5py.File(fpath, 'r')
        group = f[label]
        keys   = group['keys'][:]
        values = group['values'][:]

        data = {k: v for k, v in  zip(*[keys, values])}
        return data


    # Load a measurement from disk
    @classmethod
    def load(cls, path, label, num_leading_symbols=0):
        if os.path.isdir(path):
            if num_leading_symbols == 0:
                fpath = os.path.join(path, '%s.pickle' % label)
                count_dict  = cls._pickle_load(fpath)
            else:
                fragment_name = label[:num_leading_symbols]
                fpath         = os.path.join(path, fragment_name, '%s.pickle' % label)
                count_dict    = cls._pickle_load(fpath) 
        elif path.endswith('.hdf5'):
            fpath = path
            count_dict = cls._hdf5_load(fpath, label)
        measurement = cls(label, count_dict)
        return measurement

# --------------------------------------------- #
#    functions to calculate measurement_list    #
#       -->  good for parallel calc             #
# --------------------------------------------- #
def measurement_list_calc(label_list, count_dict_list):
    """         label_list = self.label_list
                data_dict  = self.data_dict
    """
    
    parity_flavor='effective'
    beta = None
    measurement_object_list = [Measurement(label, count_dict) for (label, count_dict) in 
                                    zip(*[label_list, count_dict_list])]

    measurement_list = [measurement_object.get_pauli_correlation_measurement(beta, parity_flavor)[label] for
                            (label, measurement_object) in zip(*[label_list,
                                                                measurement_object_list])]
    #print('label_list       = {}'.format(label_list))
    #print('measurement_list = {}'.format(measurement_list))
    return measurement_list


def measurement_list_calc_wID_wrap(argv):
    #print('       argv = {}'.format(argv))

    #ID, measurement_list = measurement_list_calc_wID(*argv)
    #return [ID, measurement_list]

    out = measurement_list_calc_wID(*argv)
    return out


def measurement_list_calc_wID(ID, label_list, count_dict_list):

    print('      start to calc    {}-th   label_list'.format(ID))
    #print('      list calc  ID = {}  -->   label_list = {}'.format(ID, label_list))
    #print('                 data_dict = {}\n'.format(count_dict_list))

    parity_flavor='effective'
    beta = None
    measurement_object_list = [Measurement(label, count_dict) for (label, count_dict) in 
                                    zip(*[label_list, count_dict_list])]

    measurement_list = [measurement_object.get_pauli_correlation_measurement(beta, parity_flavor)[label] for
                            (label, measurement_object) in zip(*[label_list,
                                                                measurement_object_list])]
    #print('label_list       = {}'.format(label_list))
    #print('measurement_list = {}'.format(measurement_list))
    #return measurement_list
    #return [ID, measurement_list]
    return {ID: measurement_list}


# ----------------------------- #
#       for parallel CPU        #
# ----------------------------- #
def split_list(x, num_parts):
	n = len(x)
	size = n // num_parts
	parts = [x[i * size: (i+1) * size] for i in range(num_parts - 1 )]
	parts.append(x[(num_parts - 1) * size:])
	return parts


# --------------------------------- #
#   class to deal with all labels   #
# --------------------------------- #

class MeasurementStore:
    def __init__(self,
                 measurement_dict):
        self.measurement_dict = measurement_dict
        
        self.labels = list(measurement_dict.keys())
        self.size   = len(self.labels)

    @classmethod
    def calc_measurement_list(cls, label_list, data_dict, method=0):

        print('    -------------------         start calc measurment_list        ------------------  \n')
        #print(' label_list = {}\n'.format(label_list))

        count_dict_list = [data_dict[label] for label in label_list]
        del data_dict

        #method = 0

        if method == 0:
            print('       ****   directly calc & save measurement_list  ')
            print('       ****   len(label_list) = {},  len(data_dict) = {}\n'.format(len(label_list), len(count_dict_list)))
            measurement_list = measurement_list_calc(label_list, count_dict_list)

        elif method == 1:           #  parallel CPU

            num_CPU = multiprocessing.cpu_count()
            if num_CPU > len(label_list):
                num_CPU = 3
            print('       ****  use parallel #CPU = {} to calc measurement_list \n'.format(num_CPU))
                        
            ID_list         = np.arange(num_CPU)
            label_part      = split_list(label_list, num_CPU)            
            count_dict_part = split_list(count_dict_list, num_CPU)
            #print('label_part = {}'.format(label_part))
            #print('count_dict_part = {}'.format(count_dict_part))

            del count_dict_list

            pool = multiprocessing.Pool(num_CPU)

            Run_parallel = 2

            if Run_parallel == 1:
                L_pt = pool.starmap(measurement_list_calc_wID, zip(ID_list, label_part, count_dict_part))
                pool.close()
                pool.join()
                #print('L_pt = {}\n'.format(L_pt))

                #ml_dict = {ID: ml for ID, ml in L_pt}      #  if return  [ID, measurement_list]            
                ml_dict = {}
                {ml_dict.update(xx) for xx in L_pt}         #  if return  {ID, measurement_list}

            elif Run_parallel == 2:

                ml_dict = {}
                mp_collect = []
                for ii, labels in enumerate(label_part):
                    print('       **************  {}-th label_part to be parallelly calc   *******'.format(ii))
                    #out = measurement_list_calc_wID(ii, labels, count_dict_part[ii])
                    #out = measurement_list_calc_wID_wrap([ii, labels, count_dict_part[ii]])
                    #print('     out = {}\n'.format(out))
                    #ml_dict.update(out)

                    mp_out = pool.apply_async(measurement_list_calc_wID_wrap, ([ii, labels, count_dict_part[ii]], ))
                    mp_collect.append(mp_out)
                pool.close()
                pool.join()
                del count_dict_part

                for xx in mp_collect:
                    out = xx.get()
                    #print('     out = {}\n'.format(out))

                    ml_dict.update(out)
                #print(' ml_dict = {}\n'.format(ml_dict))


            measurement_list = []
            for xx in ID_list:
                measurement_list += ml_dict[xx]
            #print(' ml_dict  = {}\n'.format(ml_dict))
            #print(' -->  measurement_list  = {}'.format(measurement_list))

        #measurement_list2 = measurement_list_calc(label_list, count_dict_list)
        #print(' -->  measurement_list2 = {}'.format(measurement_list2))
        
        #print(' calc method = {}'.format(method))
        #print('label_list       = {}'.format(label_list))
        #print('  calc  -->  measurement_list = {}'.format(measurement_list))
        print('    -------------         DONE of calculating measurment_list       ------------------  \n')

        return measurement_list

    @classmethod
    def Save_measurement_by_data_dict(cls, path, label_list, data_dict, method=0, Name=None, ToSave=1):

        #print('   *****  Having data_dict  -->  to calc & save measurement_list  ****')
        #print('   label_list = {}'.format(label_list))
        #print('   data_dict  = {}\n'.format(data_dict))


        tt1 = time.time()
        measurement_list = cls.calc_measurement_list(label_list, data_dict, method)
        tt2 = time.time()
        print('    ******   cls.calc_measurement_list   -->  time = {}   *****\n'.format(tt2-tt1))

        if ToSave == 1:
            tt3 = time.time()
            cls.Save_measurement_list_file(path, label_list, measurement_list, Name)
            tt4 = time.time()
            print('    ******   cls.Save_measurement_list_file   -->  time = {}   ****'.format(tt4-tt3))

        return label_list, measurement_list


    def Save_measurement_list(self, path, method=0, Name=None):
        print(' -------------  to calc & save measurement_list  -------------- \n')

        label_list = self.label_list
        data_dict  = self.data_dict

        measurement_list = self.calc_measurement_list(label_list, data_dict, method)

        self.Save_measurement_list_file(path, label_list, measurement_list, Name)

        return label_list, measurement_list

    @classmethod
    def Load_measurement_list(cls, path, Name=None):
        
        if Name == None:
            ml_file = '{}/measurement_list.pickle'.format(path)
            print(' ml_file  = {}'.format(ml_file))

            with open(ml_file, 'rb') as f:
                label_list, measurement_list = pickle.load(f)

        else:
            ml_file = '{}/measureL_{}.pickle'.format(path, Name)

            with open(ml_file, 'rb') as f:
                Name, label_list, measurement_list = pickle.load(f)

        print('  --> loading measurement_list DONE from  ml_file = {}\n'.format(ml_file))
        return label_list, measurement_list

    @classmethod
    def Save_measurement_list_file(cls, path, label_list, measurement_list, Name=None):

        print('\n')
        print('          path = {} \n'.format(path))
        print('          Name = {} \n'.format(Name))
        if Name == None:
            #print(' None')
            ml_file = '{}/measurement_list.pickle'.format(path)

            with open(ml_file, 'wb') as f:
                pickle.dump([label_list, measurement_list], f)

        else:
            ml_file = '{}/measureL_{}.pickle'.format(path, Name)

            #print('Name = {},  label_list = {}'.format(Name, label_list))
            #print('  measurement_list = {}'.format(measurement_list))

            with open(ml_file, 'wb') as f:
                pickle.dump([Name, label_list, measurement_list], f)

        print('          ***  measurement_list --> ml_file = {} is saved\n'.format(ml_file))


    def saveInOne(self, path, Name=None):
        format = 'pickle'
        if not os.path.exists(path):
            os.mkdir(path)

        label_list = sorted(self.labels)

        data_dict = {}
        for label in label_list:
            data_dict[label] = self.measurement_dict[label]

        self.label_list = label_list
        self.data_dict  = data_dict
        #print('  saveInOne  -->  labels = {}'.format(self.labels))
        #print('label_list = {}'.format(label_list))
        #print('measurement_dict = {}'.format(self.measurement_dict))
        #print('data_dict = {}'.format(data_dict))


        # --------  save files ---------------- #   
        if Name == None:
            # --------  saved file for labels  ---------------- #
            label_file = '{}/labels.pickle'.format(path)

            with open(label_file, 'wb') as f:
                pickle.dump(label_list, f)
            print('      ***   label list is stored in {} \n'.format(label_file))

            # --------  saved file for count_dict  ------------- #
            Fname = '{}/count_dict.pickle'.format(path)

            with open(Fname, 'wb') as f:
                #pickle.dump(self.measurement_dict, f)
                pickle.dump(data_dict, f)
            print('      ***   count_dict file = {} is dumped (i.e. saved) \n'.format(Fname))

        else:
            label_count_file = '{}/labs_data_{}.pickle'.format(path, Name)

            with open(label_count_file, 'wb') as f:
                pickle.dump([label_list, data_dict], f)
            print('      ***   label_list  &  count_dict_list  is stored in {} \n'.format(label_count_file))
    
    @classmethod
    def load_labels_data(cls, path, Name):
        label_count_file = '{}/labs_data_{}.pickle'.format(path, Name)
        print('    label_count_file = {}\n'.format(label_count_file))

        with open(label_count_file, 'rb') as f:
            label_list, data_dict = pickle.load(f)

        return label_list, data_dict


    # Load measurements previously saved under a disk folder
    @classmethod
    def load_OneDict(cls, path, labels=None):

        if labels == None:
            labels = cls.load_labels_from_file(path)
            #print(' loaded labels = {}'.format(labels))

        Fname = '{}/count_dict.pickle'.format(path)
        print('  --->   Loading count_dict from {}\n'.format(Fname))
        with open(Fname, 'rb') as f:
            measurement_dict = pickle.load(f)
        print('  --->   count_dict is loaded \n')

        return measurement_dict

    @classmethod
    def load_labels_from_file(cls, path):

        label_file = '{}/labels.pickle'.format(path)
        with open(label_file, 'rb') as f:
            labels = pickle.load(f)
        print('  --->   label_file = {} is loaded \n'.format(label_file))

        return labels


    def save(self, path):
        format = 'hdf5'
        if not path.endswith('.hdf5'):
            format = 'pickle'
            if not os.path.exists(path):
                os.mkdir(path)
        for label, count_dict in self.measurement_dict.items():
            Measurement(label, count_dict).save(path)

    @classmethod
    def load_labels(cls, path):
        if path.endswith('.hdf5'):
            with h5py.File(path, 'r') as f:
                labels = f.keys()
        else:
            try:
                labels = [fname.split('.')[0] for fname in os.listdir(path)]
            except:
                fragment_paths = [os.path.join(path, fragment) for fragment in os.listdir(path)]
                labels = []
                for fragment_path in fragment_paths:
                    fragment_labels = [fname.split('.')[0] for fname in os.listdir(fragment_path)]
                    labels.extend(fragment_labels)
        return labels

    
    
    # Load measurements previously saved under a disk folder
    @classmethod
    def load(cls, path, labels=None):

        if labels == None:
            labels = cls.load_labels(path)

        # checking if the store is fragmented and compute num_leading_symbols
        names = os.listdir(path)
        aname = names[0]
        apath = os.path.join(path, aname)
        if os.path.isdir(apath):
            num_leading_symbols = len(aname)
        else:
            num_leading_symbols = 0

        #print('num_leading_symbols = {}'.format(num_leading_symbols))
        #print(' apath  =  {}'.format(apath))
        #print(' is a directory ?? {}'.format(os.path.isdir(apath)))

        # load the store
        measurements = [Measurement.load(path, label, num_leading_symbols) for label in labels]
        measurement_dict = {}
        for label, measurement in zip(*[labels, measurements]):
            measurement_dict[label] = measurement.count_dict
        return measurement_dict



