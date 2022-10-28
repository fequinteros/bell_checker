
import numpy as np


def dict2array(counts, n_qubits, n_obs ):
    """

    """
    arr = np.zeros( n_obs**n_qubits )
    
    for idx in counts :    
        idx = idx.replace(' ','')
        arr[ int(idx[::-1],n_obs) ] = counts[idx]

    return arr


def verify( result, S ):
    """ 
    
    
    """
    for j, k in enumerate(result):
            if k > S:
                print("The inequality has been broken in the circuit number {} with value {}!".format(j, k))

                
def digit2key( n_obs, n_qubits ):
    """ 
    
    
    """    
    obs_conf = []
    for j in range(  n_obs**n_qubits ):
        obs_conf.append( np.base_repr( j, base=n_obs ).zfill( n_qubits ))

    return obs_conf
