#!/usr/bin/env python
# coding=utf-8
import numpy as np
import cvxpy as cp
import annexe


def optim_kdm(D, param):
    Nt= param.Nt
    N = param.N
    P= param.P
    K = param.K()
    G=[]
    constraint=[]
    for k in range(K):
        G.append(cp.Variable((N,N), PSD=True)) #On résout dans les SDP 1ere contrainte
        constraint.append(cp.sum(G[k],axis = 0) == 0) #Barycentre Nul 2e contrainte
    
    objective=0
    alpha=0
    print(D.shape)
    for i, t in enumerate(param.time_list()):
        alpha+=1
        w_t = annexe.W(param, t)
        Gt = annexe.G_t(G, w_t)
        constraint.append(Gt>>0) #Constraint SDP
        # Dt = matrix_distance(Gt)
        dG = cp.vstack( cp.diag(Gt) )
        Dt = cp.matmul(dG ,np.ones((1,N)))-2*Gt + cp.matmul(np.ones((N,1)),dG.T)

        Wi = annexe.mask_t(D[i], param)
        if not annexe.is_connected_t(D[i], Wi):
            print("ATTENTION LE GRAPHE N'EST PAS CONNECTE A t=".format(t))
            return
        
        # print(Wi*D[i]==D[i])
        W_vec = np.diag(Wi.flatten())
        alpha = (np.linalg.norm( np.matmul(W_vec, D[i].flatten()) ) )**(-2)
        objective += alpha*cp.norm( cp.matmul(W_vec, cp.vec(D[i]-Dt) ) )**2
    
    
    # objective/=alpha
    obj = cp.Minimize(objective)
    prob = cp.Problem(obj,constraint)

    try:
        prob.solve(solver=cp.SCS, verbose=True,normalize = True, max_iters=1000)
    except Exception as message:
        print(message)
        
    return G


def reconstruction_kedm(G, Y, param):
    Gproj = annexe.rankProj(G, param)
    print(np.linalg.norm(np.array(G)-np.array(Gproj)))
    X= annexe.gramtoX(Gproj, param)
    print(X.shape)
    Xaligned=annexe.aligned(Y)
    Xfin=[]
    for i_t, t in enumerate(param.time_list()):
        R=annexe.rotation(X[i_t], Xaligned[i_t, :, :param.anchor])
        Xfin.append(np.dot(R, X[i_t]))

    errorX = np.linalg.norm(Xaligned-Xfin) / np.linalg.norm(Xaligned)
    
    return X, Xaligned, np.array(Xfin), errorX

def reconstruction_basis(G, Y, param):
    Gproj = annexe.rankProj(G, param)
    print(np.linalg.norm(np.array(G)-np.array(Gproj)))
    X= annexe.gramtoXbasis(Gproj, param)
    print(X.shape)
    Xaligned=annexe.aligned(Y)
    Xfin=[]
    for i_t, t in enumerate(param.time_list()):
        R=annexe.rotation(X[i_t], Xaligned[i_t, :, :param.anchor])
        Xfin.append(np.dot(R, X[i_t]))

    errorX = np.linalg.norm(Xaligned-Xfin) / np.linalg.norm(Xaligned)
    
    return X, Xaligned, np.array(Xfin), errorX