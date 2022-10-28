# %% codecell
import numpy as np
import pandas as pd
import cvxpy as cp
import multiprocessing as mp
import scipy
import time
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

def gell_mann(dim):
    """
    Construye las matrices generalizadas de Gell-Mann en dimension arbitraria
    IN
        di: int. Dimension de la matriz
    OUT
        GellMatrices: array di x di x 3. Matrices de Gell-Mann apiladas en
                    la tercera dimension.
    """
    # Gell Mann simetricas
    A = np.eye(dim-1, dim)
    B = np.eye(dim, dim)
    B = B[0:dim, 1:dim]
    ii = 0
    SGM = np.zeros((dim, dim, int(dim*(dim-1)/2))) + 0*1j
    for k in range(dim-1):
        for jj in range(k, dim-1):
            SGM_aux = (np.kron(A[k,:].reshape(1, -1), B[:, jj].reshape(-1, 1))
                    + np.kron(A[k, :].conj().reshape(-1, 1), B[:,jj].conj().reshape(1, -1)))
            SGM[:, :, ii] = SGM_aux
            ii = ii + 1
    # # Gell Mann anti-simetricas
    ii = 0
    AGM = np.zeros((dim, dim, int(dim*(dim-1)/2))) + 0*1j
    for k in range(dim-1):
        for jj in range(k, dim-1):
            AGM_aux = (1j*np.kron(A[k,:].reshape(1, -1), B[:, jj].reshape(-1, 1))
                    - 1j*np.kron(A[k, :].conj().reshape(-1, 1), B[:,jj].conj().reshape(1, -1)))
            AGM[:, :, ii] = AGM_aux
            ii = ii + 1
    C = np.eye(dim, dim) + 0*1j
    DGM = np.zeros((dim, dim, dim - 1)) + 0*1j
    for l in range(dim-1):
        norm_const = (2/((l+1)*(l+2)))**(1/2)
        DGM_aux = - (l+1)*np.kron(C[l+1,:].reshape(1, -1), C[l+1,:].conj().reshape(-1, 1))
        for jj in range(l+1):
            DGM_aux = DGM_aux + np.kron(A[jj,:].reshape(1, -1), A[jj,:].conj().reshape(-1, 1))
        DGM[:,:,l] = norm_const*DGM_aux
    GellMatrices = np.concatenate([SGM, AGM, DGM, C.reshape(dim, dim, 1)], 2)

    return GellMatrices

def setmatrix2vec(matrix):
    """
    Dada una lista de observables construye sus representaciones vectoriales
    IN
        matrix: array di x di x N. Observables.
    OUT
        vectores: array di x N. Agrupa todos los observables
    """
    N = matrix.shape[2]
    di = matrix.shape[0]
    A = np.zeros((int(di**2), int(N)), dtype="complex")
    for k in range(N):
        A[:, k] = matrix2vec(matrix[:, :, k])[:, 0]
    return A

def matrix2vec(matrix):
    """
    Construye la representación vectorial de una matriz
    IN
        matrix: array di x di. Matriz de entrada
    OUT
        matrix_vec: array di^2 x 1. Vector de la matriz de entrada
    """

    di = matrix.shape[0]
    matrix_vec = np.reshape(matrix, (di**2, 1), "F")
    return matrix_vec

def expval(psi_sis, obs, nu_exp = None, seed=None):
    """
    Calcula el valor esperado del observable en el sistema con estado psi_sis
    IN
        psi_sis: dim x 1 or dim x dim array.
        obs: dim x dim array. Observable a medir
        nu_exp: int. Numero de mediciones.
    OUT
        eval: float. Valor esperado
    """

    # construct a base of dim x dim projectors.
    np.random.seed(seed)
    valores, vectores = np.linalg.eig(obs)
    # calculates the fidelity
    fid = fidelidad(psi_sis, vectores, nu_exp = nu_exp, seed=seed)
    eval = np.real(np.dot(valores, fid))
    return eval

def setexpval(psi_sis, obs, nu_exp = None, seed=None):
    np.random.seed(seed)
    N = obs.shape[2]
    evals = np.zeros((N))
    for k in range(N):
        evals[k] = expval(psi_sis, obs[:, :, k], nu_exp = nu_exp, seed=nu_exp)
    return evals

def fid_dens(rho, sigma):
    """
    fidelidad teorica entre dos matrices densidad
    rho = compressed sensing
    sigma = real density matrix
    """
    raiz = scipy.linalg.sqrtm(scipy.linalg.sqrtm(rho)@sigma@scipy.linalg.sqrtm(rho))
    fid = np.real(np.trace(raiz))**2 #### aveces era mayot que 1

    return fid

def gram_schmidt(matriz):
    """
    from a matrix of n vectors of dimension n calculate a base using
    Gram-Schmidt algorithm.
    IN
    matriz: d x d matrix
    OUT
    Q: d x d orthogonal matrix
    """
    R = np.zeros((matriz.shape[1], matriz.shape[1]),dtype=complex)
    Q = np.zeros(matriz.shape,dtype=complex)
    for k in range(0, matriz.shape[1]):
	    R[k, k] = np.sqrt(np.dot(matriz[:, k].conj(), matriz[:, k]))
	    Q[:, k] = matriz[:, k]/R[k, k]
	    for j in range(k+1, matriz.shape[1]):
		    R[k, j] = np.dot(Q[:, k].conj(), matriz[:, j])
		    matriz[:, j] = matriz[:, j] - R[k, j]*Q[:, k]
    return Q

def fidelidad(psi_sis, base, nu_exp = None, seed=None):
    """
    Calcula la fidelidad entre el estado del sistema psi_sis y autovectores de la base
    IN
        psi_sis: dim x 1 or dim x dim array.
        base: dim x 1 or dim x dim array. It is the base.
        nu_exp: int. Numero de mediciones.
    OUT
        prob: d  array with the projections between the state of
            the system psi_sis and the set of psi_est
    """

    # construct a base of dim x dim projectors.
    np.random.seed(seed)
    base = base/np.linalg.norm(base, axis = 0) #verifies normalization
    dim = psi_sis.shape[0]                     #system dimension
    n_base = int(base.size/base.shape[0])      #states in the base
    if n_base < dim:
        base = np.c_[base, np.random.rand(dim,dim-n_base)*1j]
        base = gram_schmidt(base)
    # calculates the fidelity
    if nu_exp == None:
        return fid_teo(psi_sis, base)
    else:
        return fid_sim(psi_sis, base, nu_exp, seed=seed)

def fid_teo(psi_sis, base):
    """
    Calcula la fidelidad teorica entre el estado del sistema psi_sis
    y la base
    IN
        psi_sis: d x 1 or d x d complex vector.
        base: d x d, a 2-d array. A base asociated to a estimator.
    OUT
        prob: par real vector. These are the projections
                    of psi_sis on the base that contains psi_est.
    """
    dim = psi_sis.shape[0] # dimension of the state
    dim_2 = int(psi_sis.size/dim) # second dimension of the state
    if dim_2 == 1:
        prob = np.absolute(np.dot(psi_sis.conj().T, base))**2
        prob = (prob/np.sum(prob)).real
        return prob[0]
    else:
        prob = np.sum(np.multiply(base.conj(), np.dot(psi_sis, base)), 0)
        prob = (prob/np.sum(prob)).real
        return prob

def fid_sim(psi_sis, base, nu_exp, seed=None):
    """
    Calcula la fidelidad teorica entre el estado del sistema psi_sis
    y la base
    IN
        psi_sis: d x 1 complex vector.
        base: d x d, a 2-d array. A base asociated to a estimator.
        nu_exp: int. Number of counts asociated to a experiment.
    OUT
        prob: par real vector. These are the projections
                    of psi_sis on the base that contains psi_est.
    """
    np.random.seed(seed)
    dim = psi_sis.shape[0] # dimension of the state
    dim_2 = int(psi_sis.size/dim) # second dimension of the state
    if dim_2 == 1:
        prob = np.absolute(np.dot(base.conj().T, psi_sis))**2
        prob = (prob/np.sum(prob)).real
        prob = np.random.multinomial(nu_exp, prob[:, 0])/nu_exp
        prob = prob/np.sum(prob)
        return prob
    else:
        prob =  np.sum(np.multiply(base.conj(), np.dot(psi_sis, base)), 0).real
        prob = np.random.multinomial(nu_exp, prob)
        prob = prob/nu_exp #normalizado sobre num experimentos
        prob = prob/np.sum(prob) #normalizado sobre probabilidades
        return prob

def vect2matrix(matrix_vec):
    """
    Construye la representación matricial de un vector
    IN
        matrix_vec: array di^2 x 1. Vector de la matriz de entrada
    OUT
        matrix: array di x di. Matriz de entrada
    """

    di = matrix_vec.shape[0]
    matrix = np.reshape(matrix_vec, (int(np.sqrt(di)), int(np.sqrt(di))),  "F")
    return matrix

def compressed_sensing(prob, A, epsilon=0.001):
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
    di2 = A.shape[0]
    print(di2)
    di = int(np.sqrt(di2))
    density = cp.Variable((di2, 1), complex=True)
    identity = np.eye(di, di) + 0*1j
    identity = matrix2vec(identity)
    objective = cp.Minimize(cp.normNuc(cp.reshape(density, (di, di), "F")))
    sa= (identity.conj().T@density)[0, 0] == 1
    print(sa)
    constraints = [(identity.conj().T@density)[0, 0] == 1 + 0*1j,
                    cp.norm2(A.conj().T@density - prob.reshape(-1, 1)) <= epsilon]
    # constraints = [cp.norm2(A.conj()@dens - prob.reshape(-1, 1)) <= epsilon]
    problema = cp.Problem(objective, constraints)
    result = problema.solve()
    #print("valor optimo: ", result)
    dens_est = vect2matrix(density.value)
    return dens_est/np.trace(dens_est) ### normalizar matriz densidad

def estado(dim, n_par, dens=False):
    """
    Matriz que genera n_par estados cuanticos arbitrarios
    Args:
        dim: int. Dimension del sistema
        n_par: int. Numero de estados que queremos
    Returns
        psi: array (dim, n_par). Matriz donde cada columna es un estado.

    """
    if dens:
        if n_par==0:
            psi = (np.random.normal(loc=0.0, scale=1.0,
                size=(dim, dim))
                + np.random.normal(loc=0.0, scale=1.0,
                size=(dim, dim))*1j)
            psi = np.dot(psi, psi.conj().T)
            psi = psi/np.trace(psi)
        else:
            psi = (np.random.normal(loc=0.0, scale=1.0,
                size=(dim, dim, n_par))
                + np.random.normal(loc=0.0, scale=1.0,
                size=(dim, dim, n_par))*1j)
            for k in range(n_par):
                psi[:, :, k] = np.dot(psi[:, :, k], psi[:, :, k].conj().T)
                psi[:, :, k] = psi[:, :, k]/np.trace(psi[:, :, k])
    else:
        if n_par==0:
            psi = (np.random.normal(loc=0.0, scale=1.0, size=(dim,))
                + np.random.normal(loc=0.0, scale=1.0, size=(dim,))*1j)
            psi = psi/np.linalg.norm(psi, axis=0)
        else:
            psi = (np.random.normal(loc=0.0, scale=1.0,
                size=(dim, n_par))
                + np.random.normal(loc=0.0, scale=1.0,
                size=(dim, n_par))*1j)
            psi = psi/np.linalg.norm(psi, axis=0)
    return psi

def tomography(dens, GellM, n_bases, potencia, iteraciones, ind_arr):
    total_matrix = []

    for bases in range(len(n_bases)):
        iters_exp = []
        for n_exper in potencia:
            all_matr, inf_matr, arr_matr = [], [], []
            gm_matr, alph_matr = [], []
            exper_dim = 2**n_exper
            bajo = exper_dim - int(sqrt(exper_dim))
            alto = exper_dim + int(sqrt(exper_dim))
            for iter in range(iteraciones):
                nu_exp = np.random.randint(bajo, high = (alto+1))
                gm = GellM[:, :, ind_arr[bases]]
                # mediciones
                obs = setmatrix2vec(gm)
                # medimos los valores esperados de los observables
                evals = setexpval(dens, gm, nu_exp = nu_exp)  #######Variar entre bases el nu_exp !!!
                # tomografia por compressed_sensing
                dens_est = compressed_sensing(evals, obs, epsilon=0.01)
                # reconstruimos la matriz unitaria a partir de la matriz de Choi
                # evaluamos la fidelidad entra la salida de compressed_sensing y la matriz
                # densidad original
                infidelidad = abs(1 - fid_dens(dens_est, dens))
                inf_matr.append(infidelidad)
                arr_matr.append(dens_est)
                gm_matr.append(obs)
                alph_matr.append(evals)
                np.random.seed()

            all_matr.append(inf_matr)
            all_matr.append(arr_matr)
            all_matr.append(gm_matr)
            all_matr.append(alph_matr)

            iters_exp.append(all_matr)
        total_matrix.append(iters_exp)
        print('Se termino la base: '+ str(bases+1))

    total_matrix = np.array(total_matrix, dtype=object)
    #np.save('complete_data_nexp_13_100.npy', total_matrix)

    return total_matrix

init_time = time.time()

#dim de mi problema
bms_dim, qudit_dim = 4, 2 #Solo bms_dim = 4 permitido por ahora
di = bms_dim*qudit_dim
Mon_iteraciones = 100
potencia = [13]
#potencia = np.arange([12])
n_bases  = np.arange(2,di*di+1)
#n_bases  = np.array([2,64])
entrada  = np.array([0, 1]) #estado

dens  = choi_matrix(bms_dim, entrada) #Solo bms_dim = 4 permitido por ahora

file  = np.load('index.npy', allow_pickle=True)
GellM = np.load('datos_GMM.npy', allow_pickle=True)

ind_arr  = []
for i in range(len(file)):
    ind  = file[i].split(", ")
    temp = []
    for j in range(len(ind)):
        temp.append(int(ind[j]))
    ind_arr.append(temp)

total_matr = tomography(dens, GellM, n_bases, potencia, Mon_iteraciones, ind_arr)

print("tiempo de computo: "+ str(time.time()-init_time))
"""
            bajo = exper_dim - int(sqrt(exper_dim))
            alto = exper_dim + int(sqrt(exper_dim))
            for monte in range(iteraciones):
                nu_exp = np.random.randint(bajo, high = (alto+1))
"""
