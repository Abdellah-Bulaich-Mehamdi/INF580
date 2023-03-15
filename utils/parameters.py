#!/usr/bin/env python
# coding=utf-8
import numpy as np

class parameters:
    def __init__(self, N, d, P, omega, start, end, taux, anchor, mode=None):
        self.N = N
        self.d = d
        self.P = P
        self.mode= mode
        self.omega=omega
        self.start = start
        self.end = end
        self.taux = taux
        self.anchor = anchor
        self.std=1
        self.path = '../../../results/kedm/python3/'
    
    def define_mode(self, mode):
        self.mode=mode
    
    def Nt(self):
        return self.K()
    
    
    def time_list(self):
        return np.linspace(self.start, self.end, self.Nt())
    
    def number_of_connections(self):
        return int(self.taux*self.N)

    def K(self):
        if self.mode == 'polynomial':
            return self.P*2+1
        elif self.mode == 'bandlimited':
            return self.P*4+1