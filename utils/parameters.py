#!/usr/bin/env python
# coding=utf-8
import numpy as np

class parameters:
    def __init__(self, N, d, P, std,  omega, start, end, Nt, taux, anchor, mode=None):
        self.N = N
        self.d = d
        self.P = P
        self.Nt= Nt
        self.mode= mode
        self.omega=omega
        self.start = start
        self.end = end
        self.taux = taux
        self.anchor = anchor
        self.std=std
        self.path = '../../../results/kedm/python3/'
    
    def define_mode(self, mode):
        self.mode=mode    
    
    def time_list(self):
        return np.linspace(self.start, self.end, self.Nt)
    
    def optim_time_list(self):
        return np.linspace(self.start, self.end, self.K())
    
    def number_of_connections(self):
        return int(self.taux*self.N)

    def K(self):
        if self.mode == 'polynomial':
            return self.P*2+1
        elif self.mode == 'bandlimited':
            return self.P*4+1