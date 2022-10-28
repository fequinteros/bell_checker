import numpy as np
import scipy


def multiport(bms): # Solo 4 permitido por ahora
    if bms == 4:
        U = np.array([[0.499 + 0.0*1j, 0.501 + 0.0*1j, 0.499 + 0.0*1j, 0.499 + 0.0*1j],
                [0.501 + 0.0*1j, 0.491 + 0.08*1j, -0.496 - 0.06*1j, -0.498 - 0.01*1j],
                [0.499 + 0.0*1j,  -0.495 - 0.06*1j, 0.498 + 0.03*1j, -0.499 + 0.03*1j],
                [0.499  + 0.0*1j, -0.499 - 0.01*1j, -0.499 + 0.03*1j, 0.499 - 0.01*1j]])
    U = U@np.linalg.inv(scipy.linalg.sqrtm(np.dot(U.conj().T, U)))

    return U

def choi_matrix(bms, entrada_canales, entrada_coefs):
    """
    Calcula la matriz de Choi de la unitaria del multiport.
    IN
        U. array di x di. Matriz del multiport.
        entrada. array n_entradas. Fibras en las que codificamos el qubit.
    OUT
        dens_out. array di*n_entradas x di*n_entradas. Matriz de Choi del proceso.
    """
    U = multiport(bms)
    n_entradas = entrada_canales.shape[0]
    dim = U.shape[0]
    diag_r = np.eye(n_entradas, n_entradas)
    diag_a = np.eye(dim, dim)
    psi = np.zeros((n_entradas*dim, 1))
    for k in range(n_entradas):
        psi = psi + entrada_coefs[k]*np.kron(diag_r[:, k], diag_a[:, int(entrada_canales[k])]).reshape(-1, 1)
    psi = psi/np.linalg.norm(psi)

    dens = np.kron(psi.conj().T, psi)
    U_comp = np.kron(diag_r, U)
    dens_out = U_comp@dens@U_comp.conj().T

    return dens_out

def random_gell_mann(GllM, dim, iteracion, base_baja, base_alta):
    """
    Entrega un conjunto aleatorio de n_bases elegidas entre las matrices de
    gellman de dimension di
    
    m = np.random.permutation(N)[:n_bases]
    """
    # se debe mejorar codigo reemplazar appends por metodos numpy
    #############################################################
    N = GllM.shape[2]


    assert base_alta <= N, "base_alta debe ser menor a N, dim del espacio de hilbert al cuadrado (AxB)**2"
    assert 0< base_baja < base_alta, "base_baja debe ser mayor que 1 y menor que base_alta"

    
    gms, indx = [], []
    for k in range(base_baja, base_alta+1):
        temp_gm, temp_indx = [], []
        for j in range(iteracion):
            m = np.arange(N, dtype=int)
            np.random.shuffle(m)
            m = m[:k]
            temp_gm.append(GllM[:, :, m])
            m.sort()
            temp_indx.append(m)

        gms.append(temp_gm)
        indx.append(temp_indx)

    gms, indx = np.array(gms), np.array(indx)
    return gms, indx

def gell_mann(dim, iteracion, base_baja, base_alta):
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
    
    gms, indx = random_gell_mann(GellMatrices, dim, iteracion, base_baja, base_alta)

    return GellMatrices, gms, indx

#dim, iteracion, base_baja, base_alta = 8, 1, 20, 31
#gell_mann(dim, iteracion, base_baja, base_alta)