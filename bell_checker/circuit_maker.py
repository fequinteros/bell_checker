# Qiskit package and modules
import qiskit
from qiskit import QuantumCircuit, transpile, Aer
from qiskit.tools.visualization import circuit_drawer
from qiskit.tools.monitor import job_monitor, backend_monitor, backend_overview

# Python packages
import numpy as np

#Our model
import bell_checker.tools as tools

class CircuitMaker:
    """ Create the circuits needed to evaluate the given inequality bell. """
    
    def __init__( self, auto=False, n_shots=1024, backend='aer_simulator' ):
        """
        Initialize the class
    
            Args:
                auto:    If the backend, run, result and get_counts are aumatized or not.
                n_shots: Number of shots when running the circuits.
                backend: The backend on the circuits are excecuted.
        """

        supported_backend = ['aer_simulator', 'qasm_simulator']

        assert type( auto )    == bool, 'The auto parameter must be a boolean (default value is False)'
        assert type( n_shots ) == int, 'The n_shots parameter must be a integer (default value is 1024)'
        assert type( backend ) == str, 'The backend parameter must be a string (default value is \'aer_simulator\')'
        assert backend in supported_backend, 'The allowed backend parameters are {}'.format(supported_backend)

        self.n_shots  = n_shots
        self.backend  = backend
        self.auto     = auto
        
        

    def constructor( self, init_state, obs, theta_matr=[0] ):
        """
        Construct the circuits needed to evaluate the bell inequality

            Args:
                init_state:  The state to evaluate in QuantumCircuit, array or list form.
                obs:         Observables to measure for each qubit.
                theta_matr:  The array or list of angles that rotate each observable, the combination of angles are given in sub-arrays (i.e [[pi]] --> rotates Alice obs stay fix and Bob obs is rotated in pi ).
        
            Output:
                The set of QuantumCircuits needed to evaluate the bell inequality.
        """
        obs = np.array( obs )
        
        n_qubits = obs.shape[0]
        n_obs    = obs.shape[1]

        self.n_qubits = n_qubits        
        self.n_obs    = n_obs

        obs_conf = tools.digit2key( n_obs, n_qubits )
        
        qcs = []
        for index in range(  n_obs**n_qubits ):
            qc = self.circuit_blueprint( init_state, obs, obs_conf, theta_matr, index )
            qcs.append(qc)
        

        self.circuits = qcs

        if self.auto == True:
            self.backend_run( )


        return qcs
 


    def circuit_blueprint( self, init_state, obs, obs_conf, theta_matr, index ):
        """
        The general form of each circuit, evaluate different cases and make the correspond circuit.
            
            Args:
                init_state:  The state to evaluate in QuantumCircuit, array or list form.
                obs:         Observables to measure for each qubit.
                obs_conf:    
                theta_matr:  The array or list of angles that rotate each observable, the combination of angles are given in sub-arrays (i.e [[pi]] --> rotates Alice obs stay fix and Bob obs is rotated in pi ).
                index:       The index of the observables to use in the QuantumCircuit
            
            Output:
                Each QuantumCircuit needed to evaluate the bell inequality.
        """
        n_qubits = self.n_qubits

        qc = QuantumCircuit( n_qubits, n_qubits )
        qc = self.init_state_circuit(qc, init_state)
        qc.barrier( )
            
        
        for k, theta in enumerate( theta_matr ):
            qc.ry( theta, k+1 )

        for j in range( n_qubits ):
            qc.unitary( obs[ j ][ int(obs_conf[index][j]) ], j, label='{}_OBS:{}'.format( j, obs_conf[index][j] ))

                       
        qc.barrier( )    
        qc.measure( range( n_qubits ), range( n_qubits ) )

        return qc
    
    def init_state_circuit( self, qc, init_state ):
        """
        Initialize the initial state with vector or QuantumCircuit form
            
            Args:
                qc:         The QuantumCircuit to append the initial state.
                init_state: The initial state to evaluate and append.

            Output:
                The inititial state in the QuantumCircuit form.
        """
        n_qubits = self.n_qubits
        if type( init_state ) == QuantumCircuit:
            qc.compose( init_state, qubits=range( n_qubits ), inplace=True )
        else:
            qc.initialize( init_state, qc.qubits )

        return qc


    def backend_run(self):
        """
        If the backend, run, result and get_counts are aumatized excecute this function

            Args:
                self: Options initialized in the class.
        """
        bell_circuits = self.circuits
        n_shots       = self.n_shots
        backend       = self.backend

        simulator = Aer.get_backend( backend )   ## externalizar y verificar si es el string propuesto o un backen typo qiskit // IBMqbackend
        job = simulator.run( bell_circuits, shots=n_shots )
        results = job.result( )
        counts = results.get_counts( )

        self.counts = counts
    


    def witness( self, Sabxy, S, counts=None, shutdown=False):
        """
        Compute the specific bell inequality to evaluate.

            Args:
                Sabxy:  The probabilitie form of the bell inequality.
                S:      The upper bound of the bell inequality.
                counts: The counts of the measurements performed for each QuantumCircuit, if auto=True not needed.

            Output:
                The values of the bell inequality for each QuantumCircuit
        """
        n_obs    = self.n_obs
        n_qubits = self.n_qubits
        
        if self.auto == True:
            counts= self.counts

        Sabxy  = np.array( Sabxy ) 
        result = []

        for i in range( 0, len(counts), n_obs**n_qubits ):
            counts_temp = counts[ i:i + n_obs**n_qubits ]
            bell_val = 0

            for j, count in enumerate( counts_temp ):
                prob     = tools.dict2array( count, n_qubits, n_obs )
                no_shots = np.sum( prob )  # total size sample
                bell_val += np.sum( prob*Sabxy[j] )/no_shots  # normalization
                
            result.append( bell_val )
        if shutdown==False:
            tools.verify( result, S )
        
        return result