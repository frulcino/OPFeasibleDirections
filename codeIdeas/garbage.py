import gurobipy as gp
from gurobipy import GRB
from gurobi_optimods import opf
from gurobi_optimods import datasets
import numpy as np
import conversion
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import time
import itertools
import math
import matplotlib.pyplot as plt
import networkx as nx
import hyperplanes
import ObjBound
import matplotlib.cm as cm

def garbage():
    nb_step=10
    theta_lim=math.pi/10
    model = gp.Model("Garbag")
    theta=model.addVar(lb=-2*theta_lim,ub=2*theta_lim)
    c=model.addVar(lb=math.cos(2*theta_lim),ub=1)
    for k in range(0,nb_step+1):
        x_bar=-2*theta_lim+k/nb_step*4*theta_lim
        a=-math.sin(x_bar)
        b=math.cos(x_bar)-x_bar*a
        # plot_line(theta_lim,a,b,0,nb_step)
        model.addConstr(c<=a*(theta)+b)
    sum_A=0
    step=4*theta_lim/nb_step
    A={}
    for k in range(0,nb_step):
        A[k]=model.addVar(vtype=GRB.BINARY, name=f"A {k}")
        sum_A+=A[k]
        x_bar=-2*theta_lim+k/nb_step*4*theta_lim
        a=(math.cos(x_bar+step)-math.cos(x_bar))/step
        b=math.cos(x_bar)-x_bar*a
        print(a,b)
        plot_line(theta_lim,a,b,k,nb_step)
        model.addConstr(c>=a*(theta)+b-(1-A[k]),name=f"c_binary {k}")
    model.addConstr(sum_A==1)

    x = np.linspace(-2*theta_lim, 2*theta_lim, 100)
    y = np.cos(x)
    plt.plot(x,y,color="b")

    plt.show()
    model.optimize()
    if model.status == GRB.INFEASIBLE:
        print("Model is infeasible.")
        # You can try to find infeasibilities by calling the following:
        model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
        model.write("infeasible_model.ilp")

def plot_line(theta_lim,a, b,k,nb_step):
    # Generate values for x (from -10 to 10)
    x = np.linspace(-2*theta_lim, 2*theta_lim, 5)
    
    # Calculate the corresponding y values for y = ax + b
    y = a * x + b
    
    # Create the plot
    plt.plot(x, y, label=f'y = {a}x + {b}',color= cm.coolwarm(k/nb_step))

garbage()