import numpy as np
import cvxpy as cp
import scipy
from math import sqrt

def multiport(bms): # Solo 4 permitido por ahora
    if bms == 4:
        U = np.array([[0.499 + 0.0*1j, 0.501 + 0.0*1j, 0.499 + 0.0*1j, 0.499 + 0.0*1j],
                [0.501 + 0.0*1j, 0.491 + 0.08*1j, -0.496 - 0.06*1j, -0.498 - 0.01*1j],
                [0.499 + 0.0*1j,  -0.495 - 0.06*1j, 0.498 + 0.03*1j, -0.499 + 0.03*1j],
                [0.499  + 0.0*1j, -0.499 - 0.01*1j, -0.499 + 0.03*1j, 0.499 - 0.01*1j]])
    U = U@np.linalg.inv(scipy.linalg.sqrtm(np.dot(U.conj().T, U)))

    return U

def choi_matrix(bms, entrada):
    """
    Calcula la matriz de Choi de la unitaria del multiport.
    IN
        U. array di x di. Matriz del multiport.
        entrada. array n_entradas. Fibras en las que codificamos el qubit.
    OUT
        dens_out. array di*n_entradas x di*n_entradas. Matriz de Choi del proceso.
    """
    U = multiport(bms)
    n_entradas = entrada.shape[0]
    dim = U.shape[0]
    diag_r = np.eye(n_entradas, n_entradas)
    diag_a = np.eye(dim, dim)
    psi = np.zeros((n_entradas*dim, 1))
    for k in range(n_entradas):
        psi = psi + np.kron(diag_r[:, k], diag_a[:, entrada[k]]).reshape(-1, 1)
    psi = psi/np.linalg.norm(psi)
    dens = np.kron(psi.conj().T, psi)
    print(psi)
    print(dens)
    U_comp = np.kron(diag_r, U)
    dens_out = U_comp@dens@U_comp.conj().T

    return dens_out, dens

#informacion a actualizar cada vez
dim_bms, qudit_dim = 4, 2
di = dim_bms*qudit_dim

entrada = np.array([0, 3]) # estado
dens_out, dens = choi_matrix(dim_bms, entrada) # (IxE)(|phi> <phi|) = |psi> <psi|, |phi> <phi|

#########################################
file = np.load('complete_data_nexp_13_100.npy', allow_pickle=True)
#print(len(file)) #num bases
#print(len(file[0])) #num pot exp
#print(len(file[0][0])) # 0 = fidelidad, 1= matrixes, 2 = gm_matrices
#print(len(file[0][0][0])) #elements of montecarlo

matr = file[19][0][1][0]
matr = np.round(matr,3)



def decompose_information(matr):
    temp_bs   = np.round(matr[0]/matr[0][0],3)/qudit_dim

    partial_bs_T = [[temp_bs[l] for l in range(len(np.real(temp_bs))//qudit_dim)],
               [temp_bs[l+dim_bms] for l in range(len(np.real(temp_bs))//qudit_dim)]]

    partial_bs   = [[partial_bs_T[f][k] for f in range(len(partial_bs_T))] for k in range(len(partial_bs_T[0]))]

    
    psi_a = np.trace(matr.reshape(qudit_dim,dim_bms,qudit_dim,dim_bms), axis1=0, axis2=2)  
    psi_dens = np.zeros(psi_a.shape[0]) + 0.0j
    for l in range(psi_a.shape[0]):        
        psi_dens[l] += psi_a[l][l]/sqrt(psi_a[l][l])
    psi_dens = psi_dens

    input_state = (np.abs(partial_bs_T)@psi_dens)/np.sum(np.abs(partial_bs_T)@psi_dens)
    input_state = np.round(input_state,3)

    return psi_a, psi_dens, input_state, partial_bs, partial_bs_T


psi_a, psi_dens, input_state, partial_bs, partial_bs_T = decompose_information(matr)


"""
print("part_1")
print(psi_a)
print('\n')
print("psi_dens")
print(psi_dens)
print('\n')
print("partial_bs")
print(partial_bs)
print(partial_bs_T)
print('\n')
print("input_state")
print(input_state)
"""