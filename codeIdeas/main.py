import CuttingPlane
import CuttingPlane_L
import CP_3L
import CP_3L_b
import conversion
import sys
import os
import time
from gurobi_optimods import opf

original_stdout = sys.stdout

local_path="C:/Users/ytsha/AppData/Local/Programs/Python/Python310/Lib/site-packages/gurobi_optimods/data/opf/"

do_plot=False
do_print=False
tightbound=False

names_to_do=[
    # "case9",
    # "case14",
    # "case30",
    # "case57",
]

memory={}
durations={}
solutions={}
solutions["case9"]=5296.69
solutions["case14"]=8081.52474
solutions["case30"]=576.89
solutions["case57"]=41737.7867
solutions["case118"]=129660.694
solutions["case_ACTIVSg500"]=93627.7080

# name="case30"
# path=local_path+name+".m"
# case = conversion.conversion(path)
# solution = opf.solve_opf(case, opftype='AC')
# print(solution["f"])

duration=time.time()
for name in names_to_do:
    print(name)

    if name not in solutions:
        solutions[name]=1

    path=local_path+name+".m"
    case = conversion.conversion(path)

    if not do_print:
        sys.stdout = open(os.devnull, 'w')
    memory[name]=[]
    durations[name]=[]
    memory[name].append(CuttingPlane.Solve(case,do_plot))
    memory[name][-1].append(memory[name][-1][0]/solutions[name])
    durations[name].append((time.time()-duration))
    duration=time.time()
    # memory[name].append(CuttingPlane_L.Solve(case,do_plot,tightbound))
    # memory[name][-1].append(memory[name][-1][0]/solutions[name])
    # durations[name].append((time.time()-duration))
    # duration=time.time()
    memory[name].append(CP_3L.Solve(case,do_plot,tightbound))
    memory[name][-1].append(memory[name][-1][0]/solutions[name])
    durations[name].append((time.time()-duration))
    duration=time.time()
    memory[name].append(CP_3L_b.Solve(case,do_plot,tightbound))
    memory[name][-1].append(memory[name][-1][0]/solutions[name])
    durations[name].append((time.time()-duration))
    duration=time.time()
    sys.stdout = original_stdout
for name in names_to_do:
    print(name+' :')
    print(memory[name])
    print(durations[name])



name="case_ACTIVSg25k"
# name="case13659pegase"
# name="case_ACTIVSg10k"
# name="case9241pegase"
# name="case1354pegase"
# name="case_ACTIVSg500"
#  name="case300"
# name="case118"
# name="case57"
# name="case30"
name="case14"
# name="case9"