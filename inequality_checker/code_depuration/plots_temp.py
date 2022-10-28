from math import comb, exp, factorial, sqrt

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import curve_fit
from scipy.stats import poisson


def multiport(bms): # Solo 4 permitido por ahora
    if bms == 4:
        U = np.array([[0.499 + 0.0*1j, 0.501 + 0.0*1j, 0.499 + 0.0*1j, 0.499 + 0.0*1j],
                [0.501 + 0.0*1j, 0.491 + 0.08*1j, -0.496 - 0.06*1j, -0.498 - 0.01*1j],
                [0.499 + 0.0*1j,  -0.495 - 0.06*1j, 0.498 + 0.03*1j, -0.499 + 0.03*1j],
                [0.499  + 0.0*1j, -0.499 - 0.01*1j, -0.499 + 0.03*1j, 0.499 - 0.01*1j]])
    U = U@np.linalg.inv(scipy.linalg.sqrtm(np.dot(U.conj().T, U)))

    return U

def image_1():
    file = np.load('recuperacion_datos_bases_20_alphas_50.npy', allow_pickle=True)
    for k in range(len(file)):
        datos = file[k]
        print(datos)

#fidelidad_estado, fidelidad_multiport, entrada_coefs, est_input_state, partial_bs_T
#image_1()

def decompose_information(matr):
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


# dim de mi problema
bms_dim, qudit_dim = 4, 2 # Solo bms_dim = 4 permitido por ahora
di = bms_dim*qudit_dim
iteraciones = 750 # iteraciones/montecarlo
potencia = [13] # potencia = np.arange(5,26)
n_bases = np.arange(20,21) # n_bases = np.arange(2,di*di+1)
entrada_canales = np.array([2, 3]) # estado

partition = 50
alphas = np.array([float(k/partition) for k in range(partition+1)])
betas  = np.array([sqrt(float(1 - alphas[k]**2)) for k in range(partition+1)]) 

bms = multiport(bms_dim)
bms_slide = [bms[:, k] for k in entrada_canales]



alphas_matr, recuperation_matr = [], []
fidelity_data, alpha_data = [], []
prom_fid, minum_fid, maxum_fid, diffs2_fid = [], [], [], []
prom_alpha, minum_alpha, maxum_alpha, diffs2_alpha = [], [], [], []
count_0_index, count_1_index = [], []

for index in range(alphas.shape[0]):
    # 20 , 35
    file_tomog = np.load('datos_bases_35_alphas_' + str(index) + '.npy', allow_pickle=True)
    entrada_coefs = [float(alphas[index]), float(betas[index])]

    temp_fid = []
    fid_temp_alpha, fid_temp_alpha_T = [], []
    #print(len(total_matr)) #num bases
    #print(len(total_matr[0])) #num pot exp
    #print(len(total_matr[0][0])) # 0 = fidelidad, 1= matrixes, 2 = gm_matrices
    #print(len(total_matr[0][0][0])) #elements of iteracion/montecarlo
    count_0, count_1 = 0,0
    for element in range(len(file_tomog[0][0][1])):
        fidelit_matr = file_tomog[0][0][0][element]
        inspect_matr = file_tomog[0][0][1][element]
        temp_fid.append(fidelit_matr)
        
        psi_a, phi_dens, est_input_state, partial_bs, partial_bs_T = decompose_information(inspect_matr)
        est_input_state_T = [est_input_state[1], est_input_state[0]]
        
        fidelidad_estado = np.dot(np.real(entrada_coefs), np.sqrt(np.real(est_input_state)))
        fidelidad_estado_T = np.dot(np.real(entrada_coefs), np.sqrt(np.real(est_input_state_T)))
        
        if index in [34,35,36]:
            print(index)
            print(entrada_coefs)
            print(est_input_state)
            print('\n')
        fid_temp_alpha.append(fidelidad_estado)
        fid_temp_alpha_T.append(fidelidad_estado_T)
        if fidelidad_estado > fidelidad_estado_T:
            count_0 += 1
        else:
            count_1 += 1
    
    if count_0 > count_1:
        count_0_index.append(index)
        prom_alpha.append(np.average(fid_temp_alpha))
        minum_alpha.append(np.max(fid_temp_alpha))
        maxum_alpha.append(np.min(fid_temp_alpha))
        diffs2_alpha.append(np.std(fid_temp_alpha))
    elif count_0 <= count_1:
        count_1_index.append(index)
        prom_alpha.append(np.average(fid_temp_alpha_T))
        minum_alpha.append(np.max(fid_temp_alpha_T))
        maxum_alpha.append(np.min(fid_temp_alpha_T))
        diffs2_alpha.append(np.std(fid_temp_alpha_T))
    else:
        pass
        #fidelidad_multiport = [np.dot(bms_slide[l], partial_bs_T[l]) for l in range(len(entrada_canales))] 
        #fidelidad_estado = np.dot(np.real(entrada_coefs), np.real(est_input_state))

        #recuperation_matr.append([fidelidad_estado, fidelidad_multiport, entrada_coefs, est_input_state, partial_bs_T])
        
    prom_fid.append(np.average(temp_fid))
    minum_fid.append(np.max(temp_fid))
    maxum_fid.append(np.min(temp_fid))
    diffs2_fid.append(np.std(temp_fid))

    fidelity_data.append(temp_fid)
    alpha_data.append(entrada_coefs)

# alphas_matr, recuperation_matr
bases_x = [x for x in range(len(alpha_data))]

fig = plt.figure(figsize=(8, 6))
gs = gridspec.GridSpec(nrows = 2,ncols=1, figure=fig, height_ratios=(10, 1),
                left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.05, hspace=0)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])

ax1.set_title("Fidelidad vs Alpha_Beta_Conf (2**13, 35b)")
ax1.set_xlabel("Fidelidad alpha-betha")
ax1.set_ylabel("Fidelidad")
ax1.set_xticks([])
ax1.set_yticks(ticks = [0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1])
ax1.set_xlim(-1,partition+1)
ax1.set_ylim(0.4, 1.1)


ax1.scatter(bases_x, prom_alpha, color='blue', linestyle='dashed', s=13, label = "prom alpha")
ax1.plot(bases_x, maxum_alpha, color='grey', lw=0.6, label= 'min/max MonC alpha')
ax1.plot(bases_x, minum_alpha, color='grey', lw=0.6)
ax1.errorbar(x=bases_x, y=prom_alpha, yerr=diffs2_alpha, fmt=' ' ,lw = 0.8, capsize = 1,
            capthick = 1, color='black')
ax1.vlines([25,35,43], ymin=0.4, ymax=1.1, lw=0.4, color='black', linestyles='dashed')

ax1.scatter(bases_x, prom_fid, color='red', linestyle='dashed', s=13, label = "prom fid")
ax1.plot(bases_x, maxum_fid, color='green', lw=0.6, label= 'min/max MonC fid')
ax1.plot(bases_x, minum_fid, color='green', lw=0.6)
ax1.errorbar(x=bases_x, y=prom_fid, yerr=diffs2_fid, fmt=' ' ,lw = 0.8, capsize = 1,
            capthick = 1, color='black')

alpha_2, beta_2 = [alphas[i]**2 for i in range(len(alphas))], [betas[i]**2 for i in range(len(betas))]
ax2.bar(bases_x, alpha_2, width=0.5)
ax2.bar(bases_x, beta_2, bottom=alpha_2, width=0.5)
ax2.hlines(0.25, xmin=-1, xmax=51, lw=0.7, color='grey', linestyles='dashed')
ax2.hlines(0.5, xmin=-1, xmax=51, lw=0.7, color='grey', linestyles='dashed')
ax2.hlines(0.75, xmin=-1, xmax=51, lw=0.7, color='grey', linestyles='dashed')

ax2.vlines([25,35,43], ymin=0, ymax=1.05, lw=0.6, color='black', linestyles='dashed')

ax2.set_yticks([])
ax2.set_xlim(-1,partition+1)
ax2.set_ylim(0, 1.05)

ax1.legend(loc ='upper right')
plt.xticks(ticks = bases_x)
plt.show()
