#!/usr/bin/env python
# coding=utf-8
import numpy as np


def generate_A(param):
    N = param.N
    d = param.d
    P = param.P
    mode=param.mode
    if mode=='polynomial':
        A=np.random.randn(P,N, d)
    elif mode == 'bandlimited':
        A=[]
        j = (-1) ** (1/2)
        for p in range(P):
            Ar = np.random.randn(N,d)
            Aq = np.random.randn(N,d)
            A.append(Ar+j*Aq)
        A.append(np.random.randn(N,d))
        for p in range(P):
            A.append(np.conj(A[P-p-1]))
        A=np.array(A)
    return A


def generate_trajectory(param, A):
    N = param.N
    d = param.d
    P = param.P
    mode=param.mode
    omega=param.omega
    Nt= param.Nt()
    j = (-1) ** (1/2)
    X = np.zeros((Nt,N,d))+0*j
    if mode =='polynomial':
        for i,t in enumerate(param.time_list()):
            for p in range(P):
                X[i] += t**p * A[p]
    elif mode == 'bandlimited':
        for i,t in enumerate(param.time_list()):
            for p in range(P+1):
                ci = np.exp(j*p*omega*t)
                X[i] += ci*A[p+P]+np.conj(ci)*A[P-p]
        X -= A[P]
    X = np.real(X)
    return X

def gram_matrix(X):
    return np.dot(X, X.T)

def matrix_distance(G):
    dG = np.diag(G)[:,None]
    D = dG-2*G+dG.T
    return D

def distance_trajectory(X):
    Nt, N, d = X.shape
    G = np.zeros((Nt, N, N))
    D = np.zeros((Nt, N, N))
    for i in range(Nt):
        G[i]= gram_matrix(X[i])
        D[i] = matrix_distance(G[i])
    return D

def W(param, t):
    #Renvoie les w_k à l'instant t
    P = param.P
    K= param.K()

    M = np.zeros((K,K))
    T = np.zeros((K,1))
    time= param.time_list()
    for i in range(K):
        M[i,:] = time**i
        T[i,0] = t**i
    w_t = np.matmul(np.linalg.inv(M), T)
    return w_t

def G_t(G, w_t):
    #renvoie G à l'instant t
    K = np.shape(w_t)[0]
    G_tot = w_t[0]*G[0]
    for k in range(K-1):
        G_tot += w_t[k+1]*G[k+1]
    return G_tot

def rankProj(G,param):
    '''
    Projette sur la dimension problème d
    On garde que les d premières valeurs propres'''
    N = param.N
    d = param.d
    K = param.K()
    Gproj = [None for k in range(K)]
    try:
        for k in range(K):
            [U,S,V] = np.linalg.svd(G[k].value,full_matrices=True)
            S[d:N] = 0
            Gproj[k] = np.matmul(np.matmul(U,np.diag(S)),V)
    except:
        for k in range(K):
            [U,S,V] = np.linalg.svd(G[k],full_matrices=True)
            S[d:N] = 0
            Gproj[k] = np.matmul(np.matmul(U,np.diag(S)),V)
    return Gproj

def rotation(X,Y):
    '''
    Return the rotation matrix with Y anchor points'''
    M =Y.shape[0]
    XY = np.dot(Y.T, X[:M])
    U,_,V = np.linalg.svd(XY,full_matrices = True)  ## XJY' = UV' = R
    R = np.dot(U,V)
    return R

def gramtoX(G,param):
    X=[]
    for t in param.time_list():
        w_t = W(param, t)
        Gt = G_t(G, w_t)
        N = Gt.shape[0]
        [U,S,V] = np.linalg.svd(Gt,full_matrices=True)
        S = S ** (1/2)
        S = S[0:param.d]
        X.append(np.matmul(np.diag(S),V[0:param.d]).T)
    return np.array(X)

def aligned(X):
    "Met le barycentre à 0"
    (Nt, N, d)=X.shape
    print(Nt, N, d)
    J = np.eye(N) - np.ones((N,N))/N
    Xaligned = np.array([np.matmul(J, X[i]) for i in range(Nt)])
    return Xaligned

def mask_t(D, param):
    M= param.number_of_connections()
    N =D.shape[0]
    mask = np.zeros((N,N))
    for i in range(N):
        index=np.argsort(D[i])
        for j in range(M):
            mask[i, index[j]]=1
    mask =np.ceil((mask+mask.T)/2)
    return mask
    
def is_connected_t(D, W):
    Connect = np.zeros((W.shape))
    N =D.shape[0]
    infty=1e30
    for i in range(N):
        for j in range(N):
            Connect[i,j]=infty
    for i in range(N):
        for j in range(i, N):
            if W[i,j]==1:
                Connect[i,j]=D[i,j]*W[i,j]
                Connect[j,i]=D[i,j]*W[i,j]
    
    for i in range(N):
        for j in range(N):
            for k in range(N):
                test= Connect[i,k]+Connect[k,j]
                if Connect[i,j]>test:
                    Connect[i,j]=test
    drapeau =0
    for i in range(N):
        for j in range(N):
            if Connect[i,j] >= infty:
                drapeau =1
    
    if drapeau ==1:
        return False
    return True