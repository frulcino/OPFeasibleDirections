#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 19:44:54 2024

@author: frulcino
"""

import gurobipy as gp
from gurobi_optimods import datasets, opf
import numpy as np
from gurobipy import GRB
import networkx as nx
from itertools import chain, combinations


class Network:
    
    def __init__(self, case):
        """
        case: matpower case in optimod format
        """
        # dictionary generator: bus since one bus can have more than one generator
        generators = {i: case["gen"][i]['bus'] for i in range(len(case["gen"]))}
        # list of buses
        buses = [case["bus"][i]['bus_i'] for i in range(len(case["bus"]))]
        bus_to_index = dict(zip(buses, np.arange(len(buses)))) #maps bus_i to the index of the corresponding element in the list case["bus"]
        # list of branches
        branches = [(case["branch"][i]['fbus'], case["branch"][i]['tbus']) for i in range(len(case["branch"]))]
        # list of branches with both orientations
        double_branches = branches + [(b, a) for (a, b) in branches]
        # various limits of variables
        Pg_min = {i: case["gen"][i]['Pmin'] for i in generators.keys()}
        Pg_max = {i: case["gen"][i]['Pmax'] for i in generators.keys()}
        Qg_min = {i: case["gen"][i]['Qmin'] for i in generators.keys()}
        Qg_max = {i: case["gen"][i]['Qmax'] for i in generators.keys()}
        V_min = {i: case["bus"][bus_to_index[i]]['Vmin'] for i in buses}
        V_max = {i: case["bus"][bus_to_index[i]]['Vmax'] for i in buses}
        # demand
        Pd = {i: case["bus"][bus_to_index[i]]['Pd'] for i in buses}
        Qd = {i: case["bus"][bus_to_index[i]]['Qd'] for i in buses}
        self.generators = generators
        self.buses = buses
        self.bus_to_index = bus_to_index
        self.branches = branches
        self.double_branches = double_branches
        self.Pg_min = Pg_min
        self.Pg_max = Pg_max
        self.Qg_min = Qg_min
        self.Qg_max = Qg_max
        self.V_min = V_min
        self.V_max = V_max
        self.Pd = Pd
        self.Qd = Qd
        self.case = case
        
    
    def create_admittance_matrix(self, branch):
        """
        Creates admittance matrix of n, since it's sparse it's saved as a dictionary of tuples (indexes)
        
        return: G, B s.t. Y = G + Bi
        
        """
        case = self.case
        G = {}
        B = {}
        #for branch in case["branch"]:
        fbus = branch["fbus"]
        tbus = branch["tbus"]
        r = branch["r"] # real part of impendence
        x = branch["x"] # complex part of impendence
        b = branch["b"] # line charging susceptance
        t = branch["ratio"]
        v = branch["angle"]
        z = r + x*1j # impendence
        y = 1/z
        A = (y + (b/2)*1j) #no trasformer admittance
        
        if t == 0: 
            t = 1 # optimod lo fa quindi io lo faccio
        
        Y22 = A 
        Y11 = A / t**2
        Y12 = - y / (t * np.exp(-v*1j))
        Y21 = - y / (t * np.exp(v*1j))
        
        G[(fbus,fbus)] = Y11.real
        G[(fbus,tbus)] = Y12.real
        G[(tbus,fbus)] = Y21.real
        G[(tbus,tbus)] = Y22.real
        
        B[(fbus,fbus)] = Y11.imag
        B[(fbus,tbus)] = Y12.imag
        B[(tbus,fbus)] = Y21.imag
        B[(tbus,tbus)] = Y22.imag
            
        return G,B
    
    class Model(self):
        """
        Initializes the OPF (Optimal Power Flow) variables for the given network.

        This function creates and initializes the variables required for the OPF problem.
        The variables include:
        - Pg: The active power generation for each generator.
        - Qg: The reactive power generation for each generator.
        - v: The voltage squared for each bus.
        - c: The current magnitude for each branch.
        - s: The current angle for each branch.
        - P: The active power flow for each double branch.
        - Q: The reactive power flow for each double branch.
        - s_abs: The absolute value of the current for each branch.
        - u: The binary variable indicating the status of each branch.
        - z_C: The binary variable indicating the status of each loop.
        - z_A: The binary variable indicating the status of each set of branches.
        - lambda_A: The binary variable indicating the status of each set of branches.
        - m_A: The integer variable indicating the number of branches in each set.

        Parameters:
            self (object): The instance of the class.

        Returns:
            m (gp.Model): The Gurobi model containing the OPF variables.
        """
        
        m = gp.Model()
        #%% Jabr variables
        Pg = m.addVars(generators.keys(), lb=Pg_min.values(), ub=Pg_max.values(),
                    vtype=GRB.CONTINUOUS, name=["Pg_{}".format(i) for i in generators])
        Qg = m.addVars(generators.keys(), lb=Qg_min.values(), ub=Qg_max.values(),
                    vtype=GRB.CONTINUOUS, name=["Qg_{}".format(i) for i in generators])
        v = m.addVars(buses, lb=[V_min[i]**2 for i in buses], ub=[V_max[i]**2 for i in buses],
                    vtype=GRB.CONTINUOUS, name=["v_{}".format(i) for i in buses])
        c = m.addVars(branches, lb=0, ub=[V_max[i]*V_max[j] for (i, j) in branches],
                    vtype=GRB.CONTINUOUS, name=["c_{}{}".format(i, j) for (i, j) in branches])
        s = m.addVars(branches, lb=[-V_max[i]*V_max[j] for (i, j) in branches], ub=[V_max[i]*V_max[j] for (i, j) in branches],
                    vtype=GRB.CONTINUOUS, name=["s_{}{}".format(i, j) for (i, j) in branches])
        P = m.addVars(double_branches, lb=-GRB.INFINITY, ub=GRB.INFINITY,
                    vtype=GRB.CONTINUOUS, name=["P_{}{}".format(i, j) for (i, j) in double_branches])
        Q = m.addVars(double_branches, lb=-GRB.INFINITY, ub=GRB.INFINITY,
                    vtype=GRB.CONTINUOUS, name=["Q_{}{}".format(i, j) for (i, j) in double_branches])
        #%% Loop constraint linearization variables

        s_abs = m.addVars(branches, lb = 0, ub = 1,
                        vtype=GRB.CONTINUOUS, name=["s_abs_{}{}".format(i,j) for (i,j) in branches])
        u = m.addVars(branches, vtype=GRB.BINARY, name=["u_{}{}".format(i,j) for (i,j) in branches])

        self.find_loops()
        loops = self.loops
        loops_b = self.loops_b
        z_C = m.addVars(np.arange(len(loops)), lb=0, ub=1, vtype=GRB.CONTINUOUS, name = ["z_cycle_{}".format(i) for i in np.arange(len(loops))])

        As_list = [[A for A in chain.from_iterable(combinations(loop, r) for r in range(0, len(loop)+1, 2))] for loop in loops_b]
        #Aset = set(sum([A for A in As] for As in As_list))
        #Aset = list(Aset)
        z_A_index = [(ii, this_A) for ii in range(len(As_list)) for this_A in As_list[ii]]
        z_A = m.addVars(z_A_index, lb=0, ub=1, vtype = GRB.CONTINUOUS, name = ["z_{}".format(A) for A in z_A_index])
        lambda_A = m.addVars(z_A_index, vtype=GRB.BINARY, name = ["lambda_{}".format(A) for A in z_A_index])
        m_A = m.addVars(z_A_index, vtype=GRB.INTEGER, lb = 0, ub = [len(A[1]) // 2 for A in z_A_index], name = ["m_{}".format(A) for A in z_A_index])
        #r_A = m.addVars(As, vtype=GRB.BINARY, name = ["r_{}".format(A) for A in As])
        return m  
    
    def find_loops(self):
        """
        Finds the loops in the network using the minimum cycle basis algorithm.

        Returns:
            loops (list): A list of ordered loops in the network. Each loop is represented as a list of nodes.
        
        """
        branches = self.branches
        double_branches = self.double_branches
        G = nx.Graph()  # building the support graph of the network
        for (i, j) in branches:
            G.add_edge(i, j)
        loops_unord = nx.minimum_cycle_basis(G)  # evaluate minimum (in terms of length if unweighted graph) cycle basis of network
        loops_ordered = []
        for loop_unord in loops_unord:
            loop_ord = [loop_unord[0]]
            loop_unord.remove(loop_unord[0])
            cont = 0
            while cont < len(loop_unord):
                if (loop_ord[-1], loop_unord[cont]) in double_branches:
                    loop_ord.append(loop_unord[cont])
                    loop_unord.remove(loop_unord[cont])
                    cont = 0
                else:
                    cont += 1
            loops_ordered.append(loop_ord)
        loops = loops_ordered

        self.graph = G
        #nx.draw(G, with_labels = True)
        # loop_v = loops[0]
        loops_b = [[(loop_v[i-1], loop_v[i]) for i in np.arange(len(loop_v)) if (loop_v[i-1], loop_v[i])  in branches] + [(loop_v[i], loop_v[i-1]) for i in np.arange(len(loop_v)) if (loop_v[i-1], loop_v[i])  not in branches] for loop_v in loops]
        self.loops = loops
        self.loops_b = loops_b

    
