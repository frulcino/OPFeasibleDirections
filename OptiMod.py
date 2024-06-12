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

#TODO che radiants stuff
# c'è da fiss

def create_admittance_matrix(case, branch):
    """
    Creates admittance matrix of n, since it's sparse it's saved as a dictionary of tuples (indexes)
    
    return: G, B s.t. Y = G + Bi
    
    """
    
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

case = datasets.load_opf_example("case9")
# result = opf.solve_opf(case, opftype="AC")

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

# Writing the model
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


# Loop constraints and extended flower inequalities
G = nx.Graph()  # building the support graph of the network
for (i, j) in branches:
    G.add_edge(i, j)
loops = nx.minimum_cycle_basis(G)  # evaluate minimum (in terms of length if unweighted graph) cycle basis of network
loop_v = loops[0]
loop = [(loop_v[i-1], loop_v[i]) for i in np.arange(len(loop_v)) if (loop_v[i-1], loop_v[i])  in branches] + [(loop_v[i], loop_v[i-1]) for i in np.arange(len(loop_v)) if (loop_v[i-1], loop_v[i])  not in branches]
# assummiamo ci sia un solo ciclo eheh
z_C = m.addVars(np.arange(len(loops)), lb=0, ub=1, vtype=GRB.CONTINUOUS, name = ["z_cycle_{}".format(i) for i in np.arange(len(loops))])
#for loop in loops:
As = [A for A in chain.from_iterable(combinations(loop, r) for r in range(0, len(loop)+1, 2))]

z_A = m.addVars(As, lb=0, ub=1, vtype = GRB.CONTINUOUS, name = ["z_{}".format(A) for A in As])
lambda_A = m.addVars(As, vtype=GRB.BINARY, name = ["lambda_{}".format(A) for A in As])
m_A = m.addVars(As, vtype=GRB.INTEGER, lb = 0, ub = [len(A) // 2 for A in As], name = ["m_{}".format(A) for A in As])
#r_A = m.addVars(As, vtype=GRB.BINARY, name = ["r_{}".format(A) for A in As])

#%% constraints
for (i, j) in branches:
    m.addConstr(c[(i, j)] ** 2 + s[(i, j)] ** 2 <= v[i] * v[j])  # (2)
    m.addConstr(s[(i, j)] == V_max[i] * V_max[j] * (2 * u[(i,j)] - 1) * s_abs[(i,j)]) # 

def StdFormRelx(z, x, indices):
    for index in indices:
        m.addConstr(z <= x[index])
    m.addConstr(z + sum(1 - x[index] for index in indices) >= 1)
    
    
for A in As:
    m.addLConstr(lambda_A[A] + 2*m_A[A] == sum(u[h] for h in A))
    monomials = [s_abs[h] for h in A] + [c[k] / (V_max[k[0]] * V_max[[k[1]]]) for k in loop if k not in A]
    StdFormRelx(z_A[A], monomials , np.arange(len(monomials)))
    

#U_C = np.prod(V_max[i]**2 for i in loop_v)
m.addConstr(sum((-1)**(len(A) // 2) * (1 - 2*lambda_A[A]) * z_A[A] for A in As) == z_C[0])
     


z = z_C
x = v
for index in loop_v:
    m.addConstr(z[0] <= x[index] / (V_max[index]**2))
m.addConstr(z[0] + sum(1 - x[index] / (V_max[index]**2) for index in loop_v) >= 1)


    


#m.addLConstr(v[1] == 1)
baseMVA = case["baseMVA"]
for branch in case["branch"]:
    i = branch["fbus"]
    j = branch["tbus"]
    G, B = create_admittance_matrix(case, branch)
    print("MIAO")
    m.addLConstr(P[(i, j)] / (baseMVA) == G[(i, i)] * v[i] + G[(i, j)] * c[(i, j)] + B[(i, j)] * s[(i, j)])   # (3)
    m.addLConstr(Q[(i, j)] / (baseMVA) == -B[(i, i)] * v[i] - B[(i, j)] * c[(i, j)] + G[(i, j)] * s[(i, j)])  # (4)
    m.addLConstr(P[(j, i)] / (baseMVA) == G[(j, j)] * v[j] + G[(j, i)] * c[(i, j)] - B[(j, i)] * s[(i, j)])   # (3)
    m.addLConstr(Q[(j, i)] / (baseMVA) == -B[(j, j)] * v[j] - B[(j, i)] * c[(i, j)] - G[(j, i)] * s[(i, j)])  # (4)

for i in buses:
    gen_at_bus = [g for g in generators.keys() if generators[g] == i]
    linked_buses = [j for j in buses if (i, j) in double_branches]
    m.addLConstr(sum(Pg[g] for g in gen_at_bus) - Pd[i] == sum(P[(i, j)] for j in linked_buses))  # (5) + (6) real
    m.addLConstr(sum(Qg[g] for g in gen_at_bus) - Qd[i] == sum(Q[(i, j)] for j in linked_buses))  # (5) + (6) imaginary

# TODO objective function
totalcost = 0
for index, gen in enumerate(case["gen"]):
    gencost = case["gencost"][index]
    if gencost["costtype"] == 2: #polynomial cost
        c = gencost["costvector"]
        if len(c) == 3:
            totalcost += c[0] * Pg[index] ** 2 + c[1] * Pg[index] + c[2]
            
m.setObjective(totalcost)

m.optimize()


#%%
branch = case["branch"][0]
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



#eheh infeasible eheh



