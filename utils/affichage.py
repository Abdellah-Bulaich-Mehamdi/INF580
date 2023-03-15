#!/usr/bin/env python
# coding=utf-8
import numpy as np
from matplotlib import pyplot as plt

def affichage(X, param):
    colors = np.random.rand(3,param.N)
    if param.d==2:
        for i in range(len(X)):
            for n in range(param.N):
                plt.plot(X[i][:,n,0], X[i][:,n,1], c = colors[:,n], linewidth=2)
            plt.show()
    elif param.d==3:
        for i in range(len(X)):
            fig = plt.figure(figsize=plt.figaspect(0.5))
            ax = fig.add_subplot(1, 2, 1, projection='3d')
            for n in range(param.N):
                ax.plot(X[i][:,n, 0], X[i][:,n, 1],X[i][:,n, 2],c = colors[:,n],linewidth=2)
                ax.scatter(X[i][0,n, 0], X[i][0,n, 1],X[i][0,n, 2],c = colors[0,n],marker='x')