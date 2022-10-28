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
    U_comp = np.kron(diag_r, U)
    dens_out = U_comp@dens@U_comp.conj().T

    return dens_out




#informacion a actualizar cada vez
dim_bms = 4
qudit_dim=2
di = dim_bms*qudit_dim
entrada = np.array([0, 1]) # estado

dens = choi_matrix(dim_bms, entrada)

#########################################
file = np.load('complete_data_nexp_13_100.npy', allow_pickle=True)
#print(len(file)) #num bases
#print(len(file[0])) #num pot exp
#print(len(file[0][0])) # 0 = fidelidad, 1= matrixes, 2 = gm_matrices
#print(len(file[0][0][0])) #elements of montecarlo

#inf = file[6][0][0][0]
matr = file[19][0][1][0]
matr= np.array(matr)
np.save("estimacion_densidad_20_bases_100_pick_random.npy", matr)

"""
with open("matrices_densidad_estimadas.txt", mode='w') as f:
    f.write('Matriz densidad original, e-5.\n')
    f.write('\n')
    for i in range(len(dens)):
        for j in range(len(dens[i])):
            if j ==0:
                line = str(np.round(dens[i][j],5)) +'   '
            elif j == (len(dens)-1):
                line += str(np.round(dens[i][j],5))
            else:
                line += str(np.round(dens[i][j],5)) +'   '

        f.write(line + '\n')
    for k in range(3):
        f.write('\n')
    f.write('Matriz densidad estimada, diferentes decimales \{3,4,5,6\}. \n')
    f.write('\n')

    for l in range(3,7):
        for i in range(len(matr)):
            for j in range(len(matr[i])):
                if j ==0:
                    line = str(np.round(matr[i][j], l))  +'   '
                elif j == (len(matr)-1):
                    line += str(np.round(matr[i][j], l))
                else:
                    line += str(np.round(matr[i][j], l)) +'   '
            f.write(line + '\n')
        for k in range(3):
            f.write('\n')

f.close()
"""
