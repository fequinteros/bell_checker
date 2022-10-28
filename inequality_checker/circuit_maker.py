# Qiskit package and modules
import qiskit
from qiskit import QuantumCircuit, transpile, Aer
from qiskit.tools.visualization import circuit_drawer
from qiskit.tools.monitor import job_monitor, backend_monitor, backend_overview

# Python packages
import itertools as it
from itertools import permutations
import numpy as np

#Our model
import tools

class CircuitMaker:
    def __init__( self, n_qubits=None, n_obs=None ):
        """


        """
        self.n_qubits    = n_qubits
        self.n_obs       = n_obs       

        
    def constructor( self, init_state, obs, theta_vec=[0] ):
        """


        """
        obs = np.array( obs )
        if self.n_qubits  == None:
            n_qubits = obs.shape[0]
            self.n_qubits = n_qubits
        if self.n_obs == None:
            n_obs      = obs.shape[1]
            self.n_obs = n_obs
        
        qcs = self.circuit_blueprint( init_state, obs, theta_vec )
  
        return qcs
    


    def circuit_blueprint(self, init_state, obs, theta_vec ):
        """


        """
        n_obs    = self.n_obs
        n_qubits = self.n_qubits
      
        obs_conf = tools.digit2key( n_obs, n_qubits )
        
        qcs = []
        for theta in theta_vec:
            for i in range(  n_obs**n_qubits ):
                qc = QuantumCircuit( n_qubits, n_qubits )
                qc.compose( init_state, qubits=range( n_qubits ), inplace=True )
                qc.barrier( )

                if theta != 0:
                    qc.ry( theta, 0 )
                
                for j in range( n_qubits ):
                    qc.unitary( obs[ j ][ int(obs_conf[i][j]) ], j, label='{}_OBS:{}'.format( j, obs_conf[i][j] ))
                            
                qc.barrier()    
                qc.measure( range( n_qubits ), range( n_qubits ) )
                qcs.append( qc )

        return qcs
    

    
    def witness( self, counts, Sabxy, S ):
        """



        """
        n_obs    = self.n_obs
        n_qubits = self.n_qubits
        
        
        
        # Divide the list of dictionaries in sets of the number of combinations// ( number of bases )**(num de qubits)
        Sabxy  = np.array(Sabxy)
        result = []
        for i in range(0, len(counts), n_obs**n_qubits):
            counts_temp = counts[i:i + n_obs**n_qubits]
            
            if Sabxy.ndim == 1: # If the inequality is given in a expected value equation   
                no_shots = sum(counts_temp[-1][y] for y in counts_temp[-1])
                    
                bell_val = 0
                for j in range(len(counts_temp)):
                    for element in counts_temp[j]:
                        parity = (-1)**(int(element[0])+int(element[1])) # (int(element[k]) for k in range(len(ement)))
                        bell_val += parity*counts_temp[j][element]*(Sabxy[j])

                result.append(bell_val/no_shots)

                for j, count in enumerate(counts_temp):
                    print(count)
                    parity = (-1)**(int(count[0])+int(count[1])) # (int(element[k]) for k in range(len(ement)))
                    bell_val += parity*counts_temp[j][count]*(Sabxy[j])
                

            elif Sabxy.ndim == 2: # If the inequality is given in a probability equation
                bell_val    = 0
                for j, count in enumerate(counts_temp):
                    prob     = tools.dict2array(count, n_qubits, n_obs ) # 2 == > Es n_qubits ; numero de qubits
                    no_shots = np.sum(prob)
                    bell_val += np.sum(prob*Sabxy[j])/no_shots

                result.append( bell_val )

        #tools.verify( result, S )
        
        return result


