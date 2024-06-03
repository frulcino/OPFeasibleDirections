#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 11:01:49 2024

@author: frulcino
"""

from gurobi_optimods import datasets, opf
import gurobipy as gp
import numpy as np

case = datasets.load_opf_example("case9")
#result = opf.solve_opf(case, opftype="AC")

m = gp.Model()

Pg = m.addVars(len(case["gen"]), name=["Pg_{}".format(i) for i in np.arange(len(case["gen"]))])
Qg = m.addVars(len(case["gen"]), name=["Qg_{}".format(i) for i in np.arange(len(case["gen"]))])
v = m.addVars(len(case["bus"]), name=["v_{}".format(bus["bus_i"]) for bus in case["bus"]])
c = m.addVars(len(case["branch"]), name=["c_{}{}".format(branch["fbus"],branch["tbus"]) for branch in case["branch"]])
s = m.addVars(len(case["branch"]), name=["s_{}{}".format(branch["fbus"],branch["tbus"]) for branch in case["branch"]])
P = m.addVars([(branch["fbus"],branch["tbus"]) for branch in case["branch"]] + [(branch["tbus"], branch["fbus"]) for branch in case["branch"]])
Q = m.addVars([(branch["fbus"],branch["tbus"]) for branch in case["branch"]] + [(branch["tbus"], branch["fbus"]) for branch in case["branch"]])