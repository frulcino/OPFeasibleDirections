#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 11:01:49 2024

@author: frulcino
"""

from gurobi_optimods import datasets, opf
import gurobipy as gp
import numpy as np

#%% helper fuctions

class Network():
    
    def __init__(self,case):
        "loads case in matpower format"
        self.network = case
        
    def create_admittance_matrix(self, branch_i):
        """
        Creates admittance matrix of n, since it's sparse it's saved as a dictionary of tuples (indexes)
        
        return: G, B s.t. Y = G + Bi
        """
        
        case = self.network
        G = {}
        B = {}
        for branch in case["branches"]:
            fbus = branch["fbus"]
            tbus = branch["tbus"]
            r = branch["r"] # real part of impendence
            x = branch["x"] # complex part of impendence
            b = branch["b"] # line charging susceptance
            t = branch["ratio"]
            v = branch["angle"]
            y = r + x*1j # impendence
            A = (1 / y + (b/2)*1j) #no trasformer admittance
            Y22 = A 
            if t == 0: 
                Y11 = A
                Y21 = - 1 / y
                Y12 = - 1 / y
            else:
                Y11 = A / t**2
                Y12 = - 1 / (y * t * np.exp(-v*1j))
                Y21 = - 1 / (y * t * np.exp(v*1j))
            
            G[(fbus,fbus)] = Y11.real
            G[(fbus,tbus)] = Y12.real
            G[(tbus,fbus)] = Y21.real
            G[(tbus,tbus)] = Y22.real
            
            B[(fbus,fbus)] = Y11.imag
            B[(fbus,tbus)] = Y12.imag
            B[(tbus,fbus)] = Y21.imag
            B[(tbus,tbus)] = Y22.imag
            
            self.G = G
            self.B = B
            
        return G,B



#i'm a potato
#%%load case

case = datasets.load_opf_example("case9")
#result = opf.solve_opf(case, opftype="AC")

#%% test
n = Network(case)

Y = n.admittance_matrix(1)
        





#%%opf


m = gp.Model()

Pg = m.addVars(len(case["gen"]), name=["Pg_{}".format(i) for i in np.arange(len(case["gen"]))])
Qg = m.addVars(len(case["gen"]), name=["Qg_{}".format(i) for i in np.arange(len(case["gen"]))])
v = m.addVars(len(case["bus"]), name=["v_{}".format(bus["bus_i"]) for bus in case["bus"]])
c = m.addVars(len(case["branch"]), name=["c_{}{}".format(branch["fbus"],branch["tbus"]) for branch in case["branch"]])
s = m.addVars(len(case["branch"]), name=["s_{}{}".format(branch["fbus"],branch["tbus"]) for branch in case["branch"]])
P = m.addVars([(branch["fbus"],branch["tbus"]) for branch in case["branch"]] + [(branch["tbus"], branch["fbus"]) for branch in case["branch"]])
Q = m.addVars([(branch["fbus"],branch["tbus"]) for branch in case["branch"]] + [(branch["tbus"], branch["fbus"]) for branch in case["branch"]])