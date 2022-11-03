
import numpy as np


def dict2array(counts, n_qubits, n_obs ):
    """
    Transform the dictionary of count in to a list.

        Args:
            counts:   The counts of the measurements performed for each QuantumCircuit, if auto=True not needed.
            n_qubits: Number of qubits or players (i.e Alice and Bob).
            n_obs:    Number of observables per qubit or player (i.e Alice and Bob)

        Output:
            Array of count elements
    """
    arr = np.zeros( n_obs**n_qubits )
    
    for idx in counts :    
        idx = idx.replace(' ','')
        arr[ int(idx[::-1],n_obs) ] = counts[idx]

    return arr



def verify( result, S ):
    """ 
    For each element of result verify if violate or not the bell inequality

        Args:
            result: The value of the bell inequality given a specific sets of measurements.
            S:      The upper bound of the bell inequality.

        Output:
            Shows what QuantumCircuit violate the bell inequality and what is the specific value of the violation.
    """

    cc = 0     
    for j, k in enumerate(result):
        if abs(k) > S:
            cc +=1
            print("The inequality has been broken in the circuit number {} with value {}!".format(j, k))
    
    if cc == 0:
        print("This state and observable combination don't violate the bell inequalitie!")
            

                
def digit2key( n_obs, n_qubits ):
    """
    Transform a number to a number base representation (i.e binary, trinary, ...) with the same length

        Args:
            n_obs:    Number of observables per qubit or player (i.e Alice and Bob)
            n_qubits: Number of qubits or players (i.e Alice and Bob).

        Output:
            Return a list of equal length in the base number of representation, of the given decimal base number.
    """
    obs_conf = []
    for j in range(  n_obs**n_qubits ):
        obs_conf.append( np.base_repr( j, base=n_obs ).zfill( n_qubits ))

    return obs_conf

