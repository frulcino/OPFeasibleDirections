#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 11:01:49 2024

@author: frulcino
"""

import gurobipy as gp
from gurobi_optimods import datasets, opf
import numpy as np
from gurobipy import GRB
import networkx as nx
from itertools import chain, combinations


def adj(i, j):
    return [[0, 0], [0, 0]], [[0, 0], [0, 0]]


case = datasets.load_opf_example("case9")
# result = opf.solve_opf(case, opftype="AC")

# dictionary generator: bus since one bus can have more than one generator
generators = {i: case["gen"][i]['bus'] for i in range(len(case["gen"]))}
# list of buses
buses = [case["bus"][i]['bus_i'] for i in range(len(case["bus"]))]
# list of branches
branches = [(case["branch"][i]['fbus'], case["branch"][i]['tbus']) for i in range(len(case["branch"]))]
# list of branches with both orientations
double_branches = branches + [(b, a) for (a, b) in branches]
# various limits of variables
Pg_min = {i: case["gen"][i]['Pmin'] for i in generators.keys()}
Pg_max = {i: case["gen"][i]['Pmax'] for i in generators.keys()}
Qg_min = {i: case["gen"][i]['Qmin'] for i in generators.keys()}
Qg_max = {i: case["gen"][i]['Qmax'] for i in generators.keys()}
V_min = {i: case["bus"][i]['Vmin'] for i in buses}
V_max = {i: case["bus"][i]['Vmax'] for i in buses}
# demand
Pd = {i: case["bus"][i]['Pd'] for i in buses}
Qd = {i: case["bus"][i]['Qd'] for i in buses}

# Writing the model
m = gp.Model()

# variables
Pg = m.addVars(generators.keys(), lb=Pg_min.values(), ub=Pg_max.values(),
               vtype=GRB.CONTINUOUS, name=["Pg_{}".format(i) for i in generators])
Qg = m.addVars(generators.keys(), lb=Qg_min.values(), ub=Qg_max.values(),
               vtype=GRB.CONTINUOUS, name=["Qg_{}".format(i) for i in generators])
v = m.addVars(buses, lb=[V_min[i]**2 for i in buses], ub=[V_max[i]**2 for i in buses],
              vtype=GRB.CONTINUOUS, name=["v_{}".format(i) for i in buses])
c = m.addVars(branches, lb=[-V_max[k]*V_max[m] for (k, m) in branches], ub=[V_max[k]*V_max[m] for (k, m) in branches],
              vtype=GRB.CONTINUOUS, name=["c_{}{}".format(i, j) for (i, j) in branches])
s = m.addVars(branches, lb=[-V_max[k]*V_max[m] for (k, m) in branches], ub=[V_max[k]*V_max[m] for (k, m) in branches],
              vtype=GRB.CONTINUOUS, name=["s_{}{}".format(i, j) for (i, j) in branches])
P = m.addVars(double_branches, lb=-GRB.INFINITY, ub=GRB.INFINITY,
              vtype=GRB.CONTINUOUS, name=["P_{}{}".format(i, j) for (i, j) in double_branches])
Q = m.addVars(double_branches, lb=-GRB.INFINITY, ub=GRB.INFINITY,
              vtype=GRB.CONTINUOUS, name=["Q_{}{}".format(i, j) for (i, j) in double_branches])

# constraints
for (k, m) in branches:
    m.addConstr(c[(k, m)] ** 2 + s[(m, k)] ** 2 <= v[k] * v[m])  # (2)


G, B = adj(k, m)
for (k, m) in branches:
    m.addLConstr(P[(k, m)] == G[(k, k)] * v[k] + G[(k, m)] * c[(k, m)] + B[(k, m)] * s[(k, m)])   # (3)
    m.addLConstr(Q[(k, m)] == -B[(k, k)] * v[k] - B[(k, m)] * c[(k, m)] + G[(k, m)] * s[(k, m)])  # (4)
    m.addLConstr(P[(m, k)] == G[(m, m)] * v[m] + G[(m, k)] * c[(k, m)] - B[(m, k)] * s[(k, m)])   # (3)
    m.addLConstr(Q[(m, k)] == -B[(m, m)] * v[m] - B[(m, k)] * c[(k, m)] - G[(m, k)] * s[(k, m)])  # (4)

for k in buses:
    gen_at_bus = [i for i in generators.keys() if generators[i] == k]
    linked_buses = [i for i in buses if (i, k) in double_branches]
    m.addLConstr(sum(Pg[j] for j in gen_at_bus) - Pd[k] == sum(P[(k, m)] for m in linked_buses))  # (5) + (6) real
    m.addLConstr(sum(Qg[j] for j in gen_at_bus) - Qd[k] == sum(Q[(k, m)] for m in linked_buses))  # (5) + (6) imaginary

# TODO objective function


# Loop constraints and extended flower inequalities
G = nx.Graph()  # building the support graph of the network
for (i, j) in branches:
    G.add_edge(i, j)
loops = nx.minimum_cycle_basis(G)  # evaluate minimum (in terms of length if unweighted graph) cycle basis of network

for loop in loops:
    As = chain.from_iterable(combinations(loop, r) for r in range(0, len(loop)+1, 2))


