import numpy as np
import scipy


def setmatrix2vec(matrix):
    """
    Dada una lista de observables (bases de medicion) construye sus representaciones vectoriales
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

def setexpval(psi_sis, obs, nu_exp = None, seed=None):
    np.random.seed(seed)
    N = obs.shape[2]
    evals = np.zeros((N))
    for k in range(N):
        evals[k] = expval(psi_sis, obs[:, :, k], nu_exp = nu_exp, seed=nu_exp)
    return evals

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

def fid_dens(rho, sigma):
    """
    fidelidad teorica entre dos matrices densidad
    rho = compressed sensing
    sigma = real density matrix
    """
    raiz = scipy.linalg.sqrtm(scipy.linalg.sqrtm(rho)@sigma@scipy.linalg.sqrtm(rho))
    fid = np.real(np.trace(raiz))**2 #### aveces era mayot que 1

    return fid

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


