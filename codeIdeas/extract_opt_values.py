import conversion
import sys
import os
import time
import networkx as nx
import numpy as np
import math

from gurobi_optimods import opf


def create_admittance_matrix(branch):
    G = {}
    B = {}
    fbus = branch["fbus"]
    tbus = branch["tbus"]
    r = branch["r"]  # real part of impedance
    x = branch["x"]  # complex part of impedance
    b = branch["b"]  # line charging susceptance
    t = branch["ratio"]
    v = branch["angle"]
    z = r + x * 1j  # impedance
    y = 1 / z
    A = (y + (b / 2) * 1j)  # no transformer admittance

    if t == 0:
        t = 1  # optimod handles this, so I'm manually setting it here

    Y22 = A
    Y11 = A / (t ** 2)
    Y12 = -y / (t * np.exp(-v * 1j))
    Y21 = -y / (t * np.exp(v * 1j))

    G[(fbus, fbus)] = Y11.real
    G[(fbus, tbus)] = Y12.real
    G[(tbus, fbus)] = Y21.real
    G[(tbus, tbus)] = Y22.real

    B[(fbus, fbus)] = Y11.imag
    B[(fbus, tbus)] = Y12.imag
    B[(tbus, fbus)] = Y21.imag
    B[(tbus, tbus)] = Y22.imag

    if fbus == tbus:
        raise ValueError("Branch from bus to itself")

    return G, B


local_path = "C:/Users/ytsha/AppData/Local/Programs/Python/Python310/Lib/site-packages/gurobi_optimods/data/opf/"
name = "case30"
path = local_path + name + ".m"

# Load the case (Case 30) from file
case = conversion.conversion(path)

# Set the voltage angle of bus 0 to 0 degrees
case["bus"][0][8] = 0  # Bus 0, voltage angle (8th column) = 0

# Solve the OPF problem with the fixed voltage angle at bus 0
solution = opf.solve_opf(case, opftype='AC', time_limit=30)

baseMVA = case["baseMVA"]

# If there are less than or equal to 500 buses, build the support graph and process loops
if len(case["bus"]) <= 500:
    G = nx.Graph()  # Build the support graph of the network
    for i in range(len(case["branch"])):
        fbus = case["branch"][i]["fbus"]
        tbus = case["branch"][i]["tbus"]
        G.add_edge(fbus, tbus)

    # Compute and print minimum cycle basis (loops)
    loops = nx.minimum_cycle_basis(G)

    for i in range(len(case["branch"])):
        fbus = case["branch"][i]["fbus"]
        tbus = case["branch"][i]["tbus"]
        for j in range(len(case["bus"])):
            if solution["bus"][j]["bus_i"] == fbus:
                fbusangle = solution["bus"][j]["Va"]
            elif solution["bus"][j]["bus_i"] == tbus:
                tbusangle = solution["bus"][j]["Va"]
        print(f"{fbus},{tbus} = {math.sin(fbusangle - tbusangle)}")

# Print the voltage angles and magnitudes at each bus
for i in range(len(case["bus"])):
    fbus = case["bus"][i]["bus_i"]
    val2 = solution["bus"][i]["Vm"]
    val = solution["bus"][i]["Va"]
    print(f"Theta {fbus} = {val:.4f}, Tension = {val2:.4f}")
