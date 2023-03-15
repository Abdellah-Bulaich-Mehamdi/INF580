#!/usr/bin/env python
# coding=utf-8
import numpy as np
from matplotlib import pyplot as plt

def affichage(X, param):
    colors = np.random.rand(3,param.N)
    if param.d==2:
        for i in range(len(X)):
            for n in range(param.N):
                plt.plot(X[i][:,0,n], X[i][:,1,n], c = colors[:,n], label = 'i= {} '.format(n), linewidth=2)
            plt.legend()
            plt.show()
    elif param.d==3:
        for i in range(len(X)):
            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(1, 2, 1, projection='3d')
            for n in range(param.N):
                ax.plot(X[i][:, 0,n], X[i][:,1,n],X[i][:,2,n],c = colors[:,n],linewidth=2)
                ax.scatter(X[i][0,0, n], X[i][0, 1, n],X[i][0,2, n],c = colors[0,n],marker='x')