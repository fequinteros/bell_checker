from math import sqrt
import numpy as np

import pipeline_tomog as pipe
import initial_generation as ingen
import math_methods as mt_m
import compressed_information as cs_i


# parametros del problema a simular
###################################
iteraciones_por_base = 750 # iteraciones/montecarlo
base_baja, base_alta = 20, 21 # 2, di*di+1
#pot_baja, pot_alta = 20, 21 # 4, 13
cantidad_bases = np.arange(base_baja, base_alta) 
potencia_mediciones = [13] #np.arange(pot_baja, pot_alta)  2**N


fineza_coef_alph_bet = 50
entrada_canales = np.array([2, 3]) # canales de entrada [0, 1, 2, 3]
qudit_dim, bms_dim = 2, 4 # Solo bms_dim = 4 permitido por ahora
dim_hil = qudit_dim*bms_dim


# Creacion de las bases de medicion y definicion coef estados de entrada
#######################################################################

#dim, iteracion, base_baja, base_alta = 8, 1, 20, 31
GellM, gms, indx = ingen.gell_mann(dim_hil, iteraciones_por_base, base_baja, base_alta) # Crear matrices de gellman
np.save('datos_GMM_' +str(base_baja)+ '_' +str(base_baja)+ '.npy', GellM)

partition = np.arange(0, fineza_coef_alph_bet)
alphas = np.round(np.sqrt(partition/(partition.shape[0]-1)), 4)
betas  = np.round(np.sqrt(1 - np.power(alphas, 2)), 4)
############################################################################


alphas_matr, recuperation_matr = [], []
for index, (alphas_ind, betas_ind) in enumerate(zip(alphas, betas)):
    entrada_coefs = [float(alphas_ind), float(betas_ind)]

    dens = ingen.choi_matrix(bms_dim, entrada_canales, entrada_coefs) # Solo bms_dim = 4 permitido por ahora
    total_matr = pipe.tomography(dens, cantidad_bases, potencia_mediciones, iteraciones_por_base, gms, indx)
    np.save('datos_bases_20_alphas_' + str(index) + '.npy', total_matr)
    #total_orden, minimos = cs_i.comparacion_bases_gm(total_matr, qudit_dim, bms_dim)
    
    #total_matr {[alphas y betas distintos][num bases]
    #[num pot exp][0 = fidelidad, 1= matrixes, 2 = gm_matrices][elements of iteracion/montecarlo]}
    
    alphas_matr.append(total_matr)
    for nbases in range(total_matr.shape[0]):
        for nexp in range(total_matr.shape[1]):
            for iter in range(total_matr.shape[3]):
                inspect_matr = total_matr[0][0][1][0]
                psi_a, phi_dens, est_input_state, partial_bs, partial_bs_T = cs_i.decompose_information(inspect_matr)


    bms = ingen.multiport(bms_dim)
    bms_slide = [bms[:, k] for k in entrada_canales]
    fidelidad_multiport = [np.dot(partial_bs_T[l], bms_slide[l]) for l in range(len(entrada_canales))] 

    fidelidad_estado = (np.absolute(np.dot(entrada_coefs, est_input_state))**2)
    recuperation_matr.append([fidelidad_estado, fidelidad_multiport, entrada_coefs, est_input_state, partial_bs_T])

    print('Se termino coeficinte: ' +str(index+1)+' de ' + str(partition.shape[0]))

np.save('recuperacion_datos_bases_20_alphas_' + str(partition) + '.npy', recuperation_matr)

