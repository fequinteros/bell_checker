import cvxpy as cp
import numpy as np
from math import sqrt
import math_methods as mtm

def compressed_sensing(prob, obs, epsilon=0.001):
    """
    Obtiene la matriz densidad a partir de la matriz con los observables en
    forma matricial y las probabilidades utilizando compressed_sensing.
    IN
        prob:array di*N. probabilidades obtenidas al medir los proyectores
        A: array di*N x di^2. Proyectores en forma vectorial
    OUT
        dens: array di x di. Matriz densidad reconstruida
    """
    # optimizacion convexa
    di2 = obs.shape[0] #base de medicion
    di = int(np.sqrt(di2))
    density = cp.Variable((di2, 1), complex=True)
    identity = np.eye(di, di) + 0*1j
    identity = mtm.matrix2vec(identity)
    objective = cp.Minimize(cp.normNuc(cp.reshape(density, (di, di), "F")))
    constraints = [(identity.conj().T@density)[0, 0] == 1 + 0*1j,
                    cp.norm2(obs.conj().T@density - prob.reshape(-1, 1)) <= epsilon]

    problema = cp.Problem(objective, constraints)
    result = problema.solve()
    dens_est = mtm.vect2matrix(density.value)
    return dens_est/np.trace(dens_est) ### normalizar matriz densidad


def decompose_information(matr, qudit_dim, bms_dim):
    temp_bs   = np.round(matr[0]/matr[0][0],3)/qudit_dim

    partial_bs_T = [[temp_bs[l] for l in range(len(np.real(temp_bs))//qudit_dim)],
               [temp_bs[l+bms_dim] for l in range(len(np.real(temp_bs))//qudit_dim)]]

    partial_bs   = [[partial_bs_T[f][k] for f in range(len(partial_bs_T))] for k in range(len(partial_bs_T[0]))]

    
    psi_a = np.trace(matr.reshape(qudit_dim,bms_dim,qudit_dim,bms_dim), axis1=0, axis2=2)  
    psi_dens = np.zeros(psi_a.shape[0]) + 0.0j
    for l in range(psi_a.shape[0]):        
        psi_dens[l] += np.real(psi_a[l][l])/sqrt(np.real(psi_a[l][l]))
    psi_dens = psi_dens

    input_state = (np.abs(partial_bs_T)@psi_dens)/np.sum(np.abs(partial_bs_T)@psi_dens)
    input_state = np.round(input_state,3)

    return psi_a, psi_dens, input_state, partial_bs, partial_bs_T


def comparacion_bases_gm(file, index):
    bases_x = []
    fid_matr, arr_matr, gm_matr = [], [], []
    for bas in range(len(file)):
        bases_x.append(int(bas)+2)
        temp_1, temp_2, temp_3 = [], [], []
        for fid in range(len(file[bas][0][0])):
            temp_1.append(file[bas][0][0][fid])
            temp_2.append(file[bas][0][1][fid])
            temp_3.append(file[bas][0][2][fid])

        fid_matr.append(temp_1)
        arr_matr.append(temp_2)
        gm_matr.append(temp_3)

    zipped_final, total_orden = [], []
    for i in range(len(gm_matr)):
        indexs = []
        gm_temp, fid_temp = [], []
        for j in range(len(gm_matr[i])):
            if j not in indexs:
                indexs.append(j)
                gm_temp.append(gm_matr[i][j])
                temp_fid = fid_matr[i][j]
                cc = 0
                if j < len(gm_matr[i]):
                    for f in range(j+1, len(gm_matr[i])):
                        if np.array_equal(gm_matr[i][j], gm_matr[i][f]) == True:
                            indexs.append(f)
                            temp_fid += fid_matr[i][f]
                            cc +=1
                if cc > 0:
                    temp_fid = temp_fid/cc
                fid_temp.append(temp_fid)

        zipped = list(zip(gm_temp, fid_temp))
        sort_zip = sorted(zipped, key = lambda x: x[1])
        total_orden.append(sort_zip)
        zipped_final.append(sort_zip[0])

    np.save('datos_bases_minimos_20_alphas_' + str(index) + '.npy', zipped_final)
    return total_orden, zipped_final

