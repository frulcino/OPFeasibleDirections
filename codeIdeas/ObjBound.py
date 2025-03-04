import gurobipy as gp
from gurobipy import GRB
from gurobi_optimods import opf
from gurobi_optimods import datasets
import numpy as np
import conversion
from tqdm import tqdm
import time
import itertools
import math
import matplotlib.pyplot as plt
import networkx as nx
import hyperplanes



def GetBounds(model,s,c,v2):
    s_up_bounds=[]
    s_low_bounds=[]
    c_up_bounds=[]
    c_low_bounds=[]
    
    for type in range(2):
        for direction in range(2):
            for bounds in range(len(s)):
                objective = gp.QuadExpr()
                if type==0:
                    objective=s[bounds]
                else:
                    objective=c[bounds]
                if direction==0:
                    model.setObjective(objective, GRB.MINIMIZE)
                else:
                    model.setObjective(objective, GRB.MAXIMIZE)        
                model.optimize()
                
                if model.status == GRB.INFEASIBLE:
                    print("Model is infeasible.")
                    model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
                    model.write("infeasible_model.ilp")
                    val=["Error"]
                else:
                    obj=model.getObjective()
                    val=[obj.getValue()]

                if type==0:
                    if direction==0:
                        s_low_bounds.append(val[-1])
                    else:
                        s_up_bounds.append(val[-1])
                else:
                    if direction==0:
                        c_low_bounds.append(val[-1])
                    else:
                        c_up_bounds.append(val[-1])
    return s_low_bounds,s_up_bounds,c_low_bounds,c_up_bounds          
            
def GetBoundsnew(model,s_new,c_new,v2):
    s_up_bounds=[]
    s_low_bounds=[]
    c_up_bounds=[]
    c_low_bounds=[]
    
    for type in range(2):
        for direction in range(2):
            for bounds in range(1,len(s_new)+1):
                T=0
                Tlim=0
                pjabr=1
                pi2=0.15
                plim=1
                Tage=5
                Tftol=5
                eps=10E-8
                eps_par=10E-5/2
                eps_ftol=10E-5
                eps_slack=10E-5
                r=0
                objective = gp.QuadExpr()
                if type==0:
                    objective=s_new[bounds]
                else:
                    objective=c_new[bounds]
                if direction==0:
                    model.setObjective(objective, GRB.MINIMIZE)
                else:
                    model.setObjective(objective, GRB.MAXIMIZE)        
                model.optimize()
                
                if model.status == GRB.INFEASIBLE:
                    print("Model is infeasible.")
                    model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
                    model.write("infeasible_model.ilp")
                    val=["Error"]
                else:
                    obj=model.getObjective()
                    val=[obj.getValue()]

                if type==0:
                    if direction==0:
                        s_low_bounds.append(val[-1])
                    else:
                        s_up_bounds.append(val[-1])
                else:
                    if direction==0:
                        c_low_bounds.append(val[-1])
                    else:
                        c_up_bounds.append(val[-1])
    return s_low_bounds,s_up_bounds,c_low_bounds,c_up_bounds          

def GetBoundsTheta(model,case,s,c,v=[]):
    s_up_bounds=[]
    s_low_bounds=[]
    c_up_bounds=[]
    c_low_bounds=[]
    v_up_bounds=[]
    v_low_bounds=[]

    list_var=[]
    if v!=[]:
        list_var=[s,c,v]
    else:
        list_var=[s,c]
    
    for type in range(len(list_var)):
        for direction in range(2):
            for bounds in tqdm(list_var[type],ncols=100,desc=f" {type},{direction}"):
                objective = gp.QuadExpr()
                if type==0:
                    objective=s[bounds]
                elif type==1:
                    objective=c[bounds]
                elif type==2:
                    objective=v[bounds]
                if direction==0:
                    model.setObjective(objective, GRB.MINIMIZE)
                else:
                    model.setObjective(objective, GRB.MAXIMIZE)        
                model.optimize()
                
                if model.status == GRB.INFEASIBLE:
                    print("Model is infeasible.")
                    model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
                    model.write("infeasible_model.ilp")
                    val=["Error"]
                else:
                    obj=model.getObjective()
                    val=[obj.getValue()]

                if type==0:
                    if direction==0:
                        s_low_bounds.append(val[-1])
                    else:
                        s_up_bounds.append(val[-1])
                elif type==1:
                    if direction==0:
                        c_low_bounds.append(val[-1])
                    else:
                        c_up_bounds.append(val[-1])
                elif type==2:
                    if direction==0:
                        v_low_bounds.append(val[-1])
                    else:
                        v_up_bounds.append(val[-1])
    s_b={}
    compteur=0
    for bounds in s:
        s_b[bounds]=[s_low_bounds[compteur],s_up_bounds[compteur]]
        compteur+=1
    c_b={}
    compteur=0
    for bounds in c:
        c_b[bounds]=[c_low_bounds[compteur],c_up_bounds[compteur]]
        compteur+=1
    if v!=[]:
        v_b={}
        compteur=0
        for bounds in v:
            v_b[bounds]=[v_low_bounds[compteur],v_up_bounds[compteur]]
            compteur+=1
        return s_b,c_b,v_b 
    return s_b,c_b

if __name__=="__main__":
    name="case9"
    path="C:/Users/ytsha/AppData/Local/Programs/Python/Python310/Lib/site-packages/gurobi_optimods/data/opf/"+str(name)+".m"
    case = conversion.conversion(path)

    # print(GetBounds(case))