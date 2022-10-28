# import qiskit tools
from itertools import count
import qiskit
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, transpile, Aer, IBMQ
from qiskit.tools.visualization import circuit_drawer
from qiskit.tools.monitor import job_monitor, backend_monitor, backend_overview

# import python stuff
import numpy as np

class inequality:
    def __init__( self, init_state, n_systems=None , bases=None, bases_share=False ):
        self.qudit       = 2           # for now only qubit states, expect to add qutrits with module threerra
        self.init_sate   = init_state  # input state to evaluate
        self.n_systems   = n_systems   # number of players (i.e  Alice, Bob, Charlie, ...)
        self.bases       = bases       # array containing all the bases for each player
        self.bases_share = bases_share # all players share or have the same bases or not
        
        assert n_systems != None            , 'must define the number of players'
        assert isinstance( n_systems, int ) , 'the number of players must be nonzero positive integer: ' + str(n_systems)
        assert n_systems > 0                , 'the number of players must be nonzero positive'
        assert bases != None                , 'must gave the observables (bases of medition)'
        
    def get_observable( self, channel_index, base_index ):
        """
        Returns the observables of choice in channel_index and base_index i.e (0,3) == Alice, base #3 
            Input:
                channel_index , shape (1, N <= n_systems)
                base_index    , shape (1, N <= n_systems, qudit**2)

            Output:
                bases_list, shape(N, N, qudit**2), N <= n_systems
        """
        bases     = self.bases

        if self.bases_share is False:
            assert channel_index.shape[0] is base_index.shape[0], 'channels and list of bases must be the same length'             
            
            zipped = sorted(zip(channel_index, base_index), key = lambda x: x[0])
            base_eval = [bases[channel][base] for (channel, base) in zipped]

        else:
            temp_bindx = base_index         
            zip(channel_index, temp_bindx).sort(axis=[0])
            base_eval  = [bases[base] for base in base_index]
    
        return base_eval  

    def obserbable( self, channel_index, base_index ):
        n_systems  = self.n_systems
        init_state = self.init_state

        assert n_systems >= channel_index.shape[0]   ,  'the lenght of the list of observables must be less or equal the numbers of players (i.e Alice, Bob, Charlie, ...)'
        assert n_systems >= np.max( channel_index )  ,  'the max index number must be less or equal the numbers of players (i.e Alice, Bob, Charlie, ...)'
        assert 0         <= np.min( channel_index )  ,  'the min index number must be a positive number'

        base_eval = self.get_observable( channel_index, base_index )
        
        qr = QuantumCircuit( n_systems, 'q' )
        cr = ClassicalRegister( n_systems, 'c' )
        qc = QuantumCircuit( qr, cr )
        
        qc.compose( init_state, qubits=range(n_systems), inplace=True )
        qc.barrier( )
        
        for i, obs in enumerate( base_eval ):
            qc.unitary( obs, i, label='$M_{}$'.format(i) )
        
        qc.measure( range(qr), range(cr) )
        return qc

            
class checker:
    def __init__(self):
        pass
    
    def witness(self, counts, w_quant):
        """Computes expectation values for the inequality, for each
            angle (theta) between measurement axis.

                Args: counts (list[dict]): dict of counts for each experiment

                Returns:
                    Tuple(List, List): Tuple of lists with the witnesses
        """    
        witness = [[] for k in range(w_quant)]
        

        for i in range(0, len(counts), w_quant):
            dict_con = counts[i: i+ w_quant]
            wit_list = np.zeros(w_quant)
            for gates in dict_con:
                for element in gates:
                    for j, function in enumerate(w_quant):
                        wit_list[j] += function(element)
            
            for l in range(w_quant):
                witness[l].apppend(wit_list[l])
        
        return witness
    
#print(['{}'.format(i) for i in range(10)])