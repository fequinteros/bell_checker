import numpy as np
import math_methods as mt_m
import compressed_information as cs_i

def tomography(dens, n_bases, potencia, iteraciones, gm, ind):
    total_matrix = []

    for bases in range(n_bases.shape[0]):
        iters_exp = []
        for exper in potencia:
            all_matr, inf_matr, arr_matr, gm_matr = [], [], [], []
            nu_exp = 2**exper
            for iter in range(iteraciones):                
                # mediciones
                obs = mt_m.setmatrix2vec(gm[bases][iter]) # vectores: array di*di x N. Agrupa todos los observables
                
                # medimos los valores esperados de los observables
                evals = mt_m.setexpval(dens, gm[bases][iter], nu_exp = nu_exp)
                
                # tomografia por compressed_sensing
                dens_est = cs_i.compressed_sensing(evals, obs, epsilon=0.01)
                
                # reconstruimos la matriz unitaria a partir de la matriz de Choi evaluamos la
                # fidelidad entra la salida de compressed_sensing y la matriz densidad original.
                infidelidad = abs(1 - mt_m.fid_dens(dens_est, dens))

                inf_matr.append(infidelidad)
                arr_matr.append(dens_est)
                gm_matr.append(ind[bases][iter])

                np.random.seed()

            all_matr.append(inf_matr)
            all_matr.append(arr_matr)
            all_matr.append(gm_matr)
            iters_exp.append(all_matr)
        total_matrix.append(iters_exp)
        print('Se termino la base: '+ str(bases+1))

    total_matrix = np.array(total_matrix, dtype=object)
    
    

    return total_matrix
