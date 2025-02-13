#%%
import gurobipy as gp
from gurobi_optimods import datasets, opf
import numpy as np
from gurobipy import GRB
import networkx as nx
from itertools import chain, combinations
import scipy.io
import re
import pandas as pd
#%% 

class network:

    def __init__(self):
        pass

    def load_case(self, case):
        self.case = case
        # dictionary generator: bus since one bus can have more than one generator
        self.generators = {i: case["gen"][i]['bus'] for i in range(len(case["gen"]))}
        generators = self.generators
        # list of buses
        self.buses = [case["bus"][i]['bus_i'] for i in range(len(case["bus"]))]
        buses = self.buses
        bus_to_index = dict(zip(self.buses, np.arange(len(buses)))) #maps bus_i to the index of the corresponding element in the list case["bus"]
        self.bus_to_index = bus_to_index 
        # list of branches
        self.branches = [(case["branch"][i]['fbus'], case["branch"][i]['tbus']) for i in range(len(case["branch"]))]
        branches = self.branches
        # list of branches with both orientations
        self.double_branches = branches + [(b, a) for (a, b) in branches]
        double_branches = self.double_branches
        # various limits of variables
        self.Pg_min = {i: case["gen"][i]['Pmin'] for i in generators.keys()}
        self.Pg_max = {i: case["gen"][i]['Pmax'] for i in generators.keys()}
        self.Qg_min = {i: case["gen"][i]['Qmin'] for i in generators.keys()}
        self.Qg_max = {i: case["gen"][i]['Qmax'] for i in generators.keys()}
        self.V_min = {i: case["bus"][bus_to_index[i]]['Vmin'] for i in buses}
        self.V_max = {i: case["bus"][bus_to_index[i]]['Vmax'] for i in buses}
        # demand
        self.Pd = {i: case["bus"][bus_to_index[i]]['Pd'] for i in buses}
        self.Qd = {i: case["bus"][bus_to_index[i]]['Qd'] for i in buses}

        # Loop constraints and extended flower inequalities
        G = nx.Graph()  # building the support graph of the network
        for (i, j) in self.branches:
            G.add_edge(i, j)
        loops_unord = nx.minimum_cycle_basis(G)  # evaluate minimum (in terms of length if unweighted graph) cycle basis of self
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
        self.loops = loops_ordered
        loops = self.loops
        nx.draw(G, with_labels = True)
        # loop_v = loops[0]
        self.loops_b = [[(loop_v[i-1], loop_v[i]) for i in np.arange(len(loop_v)) if (loop_v[i-1], loop_v[i])  in branches] + [(loop_v[i], loop_v[i-1]) for i in np.arange(len(loop_v)) if (loop_v[i-1], loop_v[i])  not in branches] for loop_v in loops]

    def init_base_model(self):
        generators = self.generators
        buses = self.buses
        bus_to_index = self.bus_to_index 
        branches = self.branches
        double_branches = self.double_branches
        Pg_min = self.Pg_min
        Pg_max = self.Pg_max
        Qg_min = self.Qg_min
        Qg_max = self.Qg_max
        V_min = self.V_min
        V_max = self.V_max
        Pd = self.Pd
        Qd = self.Qd

        # Writing the model
        m = gp.Model()
        self.m = m
        # Jabr variables
        Pg = m.addVars(generators.keys(), lb=Pg_min.values(), ub=Pg_max.values(),
                    vtype=GRB.CONTINUOUS, name=["Pg_{}".format(i) for i in generators])
        Qg = m.addVars(generators.keys(), lb=Qg_min.values(), ub=Qg_max.values(),
                    vtype=GRB.CONTINUOUS, name=["Qg_{}".format(i) for i in generators])
        v = m.addVars(buses, lb=[V_min[i]**2 for i in buses], ub=[V_max[i]**2 for i in buses],
                    vtype=GRB.CONTINUOUS, name=["v_{}".format(i) for i in buses])
        node_couples = list(set(branches))
        c = m.addVars(node_couples, lb=0, ub=[V_max[i]*V_max[j] for (i, j) in node_couples],
                    vtype=GRB.CONTINUOUS, name=["c_{}{}".format(i, j) for (i, j) in node_couples])
        s = m.addVars(node_couples, lb=[-V_max[i]*V_max[j] for (i, j) in node_couples], ub=[V_max[i]*V_max[j] for (i, j) in node_couples],
                    vtype=GRB.CONTINUOUS, name=["s_{}{}".format(i, j) for (i, j) in node_couples])
        P = m.addVars([i for i in range(len(double_branches))], lb=-GRB.INFINITY, ub=GRB.INFINITY,
                    vtype=GRB.CONTINUOUS, name=["P_{}{}{}".format(i, j, k) for k, (i, j) in enumerate(double_branches)])
        Q = m.addVars([i for i in range(len(double_branches))], lb=-GRB.INFINITY, ub=GRB.INFINITY,
                    vtype=GRB.CONTINUOUS, name=["Q_{}{}{}".format(i, j, k) for k, (i, j) in enumerate(double_branches)])

        self.Pg = Pg
        self.Qg = Qg
        self.v = v
        self.c = c
        self.s = s
        self.P = P
        self.Q = Q

    def init_bouquet_model(self, warm_start = False, MIPGap = None, TimeLimit = None):
        generators = self.generators
        buses = self.buses
        branches = self.branches
        loops = self.loops
        loops_b = self.loops_b
        double_branches = self.double_branches
        self.init_base_model()
        m = self.m
        Pg = self.Pg
        Qg = self.Qg
        v = self.v
        c = self.c
        s = self.s
        P = self.P
        Q = self.Q
        Pd = self.Pd
        Qd = self.Qd
        Pg_min = self.Pg_min
        Pg_max = self.Pg_max
        Qg_min = self.Qg_min
        Qg_max = self.Qg_max
        V_min = self.V_min
        V_max = self.V_max
        node_couples = list(set(branches))
        s_abs = m.addVars(node_couples, lb = 0, ub = 1,
                        vtype=GRB.CONTINUOUS, name=["s_abs_{}{}".format(i,j) for (i,j) in node_couples])
        u = m.addVars(node_couples, vtype=GRB.BINARY, name=["u_{}{}".format(i,j) for (i,j) in node_couples])
        z_C = m.addVars(np.arange(len(loops)), lb=0, ub=1, vtype=GRB.CONTINUOUS, name = ["z_cycle_{}".format(i) for i in np.arange(len(loops))])
        #for loop in loops:
        As_list = [[A for A in chain.from_iterable(combinations(loop, r) for r in range(0, len(loop)+1, 2))] for loop in loops_b]
        z_A_index = [(ii, this_A) for ii in range(len(As_list)) for this_A in As_list[ii]]
        z_A = m.addVars(z_A_index, lb=0, ub=1, vtype = GRB.CONTINUOUS, name = ["z_{}".format(A) for A in z_A_index])
        lambda_A = m.addVars(z_A_index, vtype=GRB.BINARY, name = ["lambda_{}".format(A) for A in z_A_index])
        m_A = m.addVars(z_A_index, vtype=GRB.INTEGER, lb = 0, ub = [len(A[1]) // 2 for A in z_A_index], name = ["m_{}".format(A) for A in z_A_index])  
        # constraints
        for (i, j) in node_couples:
            m.addConstr(c[(i,j)] ** 2 + s[(i,j)] ** 2 <= v[i] * v[j])  # (2)
            m.addConstr(s[(i,j)] == V_max[i] * V_max[j] * (2 * u[(i,j)] - 1) * s_abs[(i,j)]) # 

        def StdFormRelx(z, x, indices):
            for index in indices:
                m.addConstr(z <= x[index])
            m.addConstr(z + sum(1 - x[index] for index in indices) >= 1)
        


        #m.addLConstr(v[1] == 1)
        baseMVA = case["baseMVA"]
        for k, branch in enumerate(case['branch']):
            i = branch["fbus"]
            j = branch["tbus"]
            G, B = network.create_admittance_matrix(branch)
            # print("MIAO")
            m.addLConstr(P[k] / (baseMVA) == G[(i, i)] * v[i] + G[(i, j)] * c[(i,j)] + B[(i, j)] * s[(i,j)])   # (3)
            m.addLConstr(Q[k] / (baseMVA) == -B[(i, i)] * v[i] - B[(i, j)] * c[(i,j)] + G[(i, j)] * s[(i,j)])  # (4)
            m.addLConstr(P[k + len(branches)] / (baseMVA) == G[(j, j)] * v[j] + G[(j, i)] * c[(i,j)] - B[(j, i)] * s[(i,j)])   # (3)
            m.addLConstr(Q[k + len(branches)] / (baseMVA) == -B[(j, j)] * v[j] - B[(j, i)] * c[(i,j)] - G[(j, i)] * s[(i,j)])  # (4)

        for i in buses:
            gen_at_bus = [g for g in generators.keys() if generators[g] == i]
            linked_buses = [j for j in buses if (i, j) in double_branches]
            # m.addLConstr(sum(Pg[g] for g in gen_at_bus) - Pd[i] == sum(P[(i, j)] for j in linked_buses))  # (5) + (6) real
            # m.addLConstr(sum(Qg[g] for g in gen_at_bus) - Qd[i] == sum(Q[(i, j)] for j in linked_buses))  # (5) + (6) imaginary
            possible_ks = []
            for k in range(len(double_branches)):
                if double_branches[k][0] == i and double_branches[k][1] in linked_buses:
                    possible_ks.append(k)
            m.addLConstr(sum(Pg[g] for g in gen_at_bus) - Pd[i] == sum(P[k] for k in possible_ks))  # (5) + (6) real
            m.addLConstr(sum(Qg[g] for g in gen_at_bus) - Qd[i] == sum(Q[k] for k in possible_ks))  # (5) + (6) imaginary

        totalcost = 0
        for index, gen in enumerate(case["gen"]):
            gencost = case["gencost"][index]
            if gencost["costtype"] == 2: #polynomial cost
                costvector = gencost["costvector"]
                if len(costvector) == 3:
                    totalcost += costvector[0] * Pg[index] ** 2 + costvector[1] * Pg[index] + costvector[2]
        
                for i,A in z_A_index:
                    # TODO qui non so come chiamare s_abs perché non so che archi vengono pescati quando si fa il loop
                    m.addLConstr(lambda_A[(i,A)] + 2*m_A[(i,A)] == sum(u[h] for h in A))
                    # print(i, A, loops_b[i])
                    monomials = [s_abs[h] for h in A] + [c[k] / (V_max[k[0]] * V_max[k[1]]) for k in loops_b[i] if k not in A]
                    StdFormRelx(z_A[(i,A)], monomials , np.arange(len(monomials)))

        m.setObjective(totalcost)   
                
            
        if warm_start:
            print("Warmstaring with Jabr Relaxation")
            m.Params.LPwarmStart = 1
            m.optimize()
            print("finished warmstart")
  

        #U_C = np.prod(V_max[i]**2 for i in loop_v)
        #loop constraint
        for i, As in enumerate(As_list):
            m.addConstr(sum((-1)**(len(A) // 2) * (1 - 2*lambda_A[(i,A)]) * z_A[(i,A)] for A in As) == z_C[i])



        for i, loop in enumerate(loops):
            z = z_C[i]
            x = v
            for index in loop:
                m.addConstr(z <= x[index] / (V_max[index]**2)) #ARTHUR Should be z <= x[index]*(prod_index(V_max[index]))
            m.addConstr(z+ sum(1 - x[index] / (V_max[index]**2) for index in loop) >= 1) #ARTHUR z+sum(V_max[index]-x[index])>= prod_index(V_max[index]) ?

        
        if MIPGap is not None:
            print(f"Setting MIPGAP to {MIPGap}")
            m.Params.MIPGap = MIPGap    # 0.14%
        if TimeLimit is not None:
            print(f"Setting TimeLimit to {TimeLimit}")
            m.Params.TimeLimit = 3000  # 5 minutes
    
        self.s_abs = s_abs
        self.u = u
        self.z_C = z_C
        self.z_A = z_A
        self.lambda_A = lambda_A
        self.m_A = m_A

    def init_jabr_model(self):
        generators = self.generators
        buses = self.buses
        branches = self.branches
        loops = self.loops
        loops_b = self.loops_b
        double_branches = self.double_branches
        self.init_base_model()
        m = self.m
        Pg = self.Pg
        Qg = self.Qg
        v = self.v
        c = self.c
        s = self.s
        P = self.P
        Q = self.Q
        Pd = self.Pd
        Qd = self.Qd
        Pg_min = self.Pg_min
        Pg_max = self.Pg_max
        Qg_min = self.Qg_min
        Qg_max = self.Qg_max
        V_min = self.V_min
        V_max = self.V_max
        
        node_couples = list(set(branches))

        for (i, j) in node_couples:
            m.addConstr(c[(i,j)] ** 2 + s[(i,j)] ** 2 <= v[i] * v[j])  # (2)
            # m.addConstr(s[(i, j)] == V_max[i] * V_max[j] * (2 * u[(i,j)] - 1) * s_abs[(i,j)]) # 

        #m.addLConstr(v[1] == 1)
        baseMVA = case["baseMVA"]
        for k, branch in enumerate(case['branch']):
            i = branch["fbus"]
            j = branch["tbus"]
            G, B = network.create_admittance_matrix(branch)
            # print("MIAO")
            m.addLConstr(P[k] / (baseMVA) == G[(i, i)] * v[i] + G[(i, j)] * c[(i,j)] + B[(i, j)] * s[(i,j)])   # (3)
            m.addLConstr(Q[k] / (baseMVA) == -B[(i, i)] * v[i] - B[(i, j)] * c[(i,j)] + G[(i, j)] * s[(i,j)])  # (4)
            m.addLConstr(P[k + len(branches)] / (baseMVA) == G[(j, j)] * v[j] + G[(j, i)] * c[(i,j)] - B[(j, i)] * s[(i,j)])   # (3)
            m.addLConstr(Q[k + len(branches)] / (baseMVA) == -B[(j, j)] * v[j] - B[(j, i)] * c[(i,j)] - G[(j, i)] * s[(i,j)])  # (4)

        for i in buses:
            gen_at_bus = [g for g in generators.keys() if generators[g] == i]
            linked_buses = [j for j in buses if (i, j) in double_branches]
            # m.addLConstr(sum(Pg[g] for g in gen_at_bus) - Pd[i] == sum(P[(i, j)] for j in linked_buses))  # (5) + (6) real
            # m.addLConstr(sum(Qg[g] for g in gen_at_bus) - Qd[i] == sum(Q[(i, j)] for j in linked_buses))  # (5) + (6) imaginary
            possible_ks = []
            for k in range(len(double_branches)):
                if double_branches[k][0] == i: # and double_branches[k][1] in linked_buses
                    possible_ks.append(k)
            m.addLConstr(sum(Pg[g] for g in gen_at_bus) - Pd[i] == sum(P[k] for k in possible_ks))  # (5) + (6) real
            m.addLConstr(sum(Qg[g] for g in gen_at_bus) - Qd[i] == sum(Q[k] for k in possible_ks))  # (5) + (6) imaginary

        # TODO objective function
        totalcost = 0
        for index, gen in enumerate(case["gen"]):
            gencost = case["gencost"][index]
            if gencost["costtype"] == 2: #polynomial cost
                costvector = gencost["costvector"]
                if len(costvector) == 3:
                    totalcost += costvector[0] * Pg[index] ** 2 + costvector[1] * Pg[index] + costvector[2]
                    
        m.setObjective(totalcost)

    @staticmethod
    def create_admittance_matrix(branch):
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

#%%
def read_matpower_case(file_path):
    """
    Reads a MATPOWER case file in .m format and parses it into a dictionary of Pandas DataFrames.
    
    Args:
        file_path (str): Path to the MATPOWER case file.
        
    Returns:
        dict: A dictionary containing DataFrames for 'baseMVA', 'bus', 'branch', 'gen', etc.
    """
    data = {}
    column_headers = {}
    current_section = None
    
    # Regular expressions for identifying sections and their boundaries
    section_pattern = re.compile(r"mpc\.(\w+)\s*=\s*\[")
    end_pattern = re.compile(r"\];")
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            
            # Check for new section
            section_match = section_pattern.match(line)
            if section_match:
                current_section = section_match.group(1)
                print(
                    f"Found new section: {current_section}"    
                )
                data[current_section] = []
                continue
            
            # Detect and process column headers (lines starting with '%' in relevant sections)
            if line.startswith('% '):
                print('found colum')
                # Capture header if it looks like column names
                if not column_headers.get(current_section):
                    print(current_section)
                    columns = re.findall(r'\S+', line[1:])
                continue
            
            # Check for end of section
            if end_pattern.match(line):
                column_headers[current_section] = columns
                current_section = None
                continue
            
            # Collect data for the current section
            if current_section and line and not line.startswith('%'):
                line_data = [float(x) if x.replace('.', '', 1).replace('-', '', 1).isdigit() else x
                             for x in re.split(r'\s+', line.rstrip(';'))]
                data[current_section].append(line_data)
    
    # Convert each section to a Pandas DataFrame
    for key, rows in data.items():
        columns = column_headers.get(key, [f"col_{i}" for i in range(len(rows[0]))])
        if len(columns) > len(rows[0]):
            columns = colummns[:len(rows[0])]
        else:
            columns += [f"col_{i}" for i in range(len(columns),len(rows[0]))]
        data[key] = pd.DataFrame(rows, columns=columns)
        print("miao")
    
    return data, column_headers



#%%
def safe_float(x):
    try:
        return float(x)
    except ValueError:
        return np.nan  # Replace invalid entries with NaN

def numpy_float_array(string_array):
    """
    Convert a 2D array of strings to a 1D array of floats, with NaN values for invalid entries.
    """    
    return np.vectorize(safe_float)(string_array)


# %%

case = datasets.load_opf_example("case300")

# %%
nj = network()
nj.load_case(case)
nj.init_jabr_model()
nj.m.optimize()

#%%
n = network()
n.load_case(case)
n.init_bouquet_model(warm_start=True, MIPGap=0.0014, TimeLimit=3000)
n.m.optimize()

# %%
# %%
def distance_from_optimal(n, matpower_solution):
    #confront with solution

    data =  matpower_solution
    Pg = np.array([n.Pg[i].X for i in n.generators.keys()])
    Qg = np.array([n.Qg[i].X for i in n.generators.keys()])
    v = np.sqrt(np.array([n.v[i].X for i in n.v.keys()]))
    Pg0 = data['gen']['Pg'].values
    Qg0 = data['gen']['Qg'].values
    v0 = data['bus']['Vm'].values

    dist = np.sqrt(np.linalg.norm(Pg - Pg0)**2 + np.linalg.norm(Qg - Qg0)**2) #+ np.linalg.norm(v - v0)
    #print(dist)
    return dist
# %%
#data =  matpower_solution
Pg = np.array([nj.Pg[i].X for i in n.generators.keys()])
Qg = np.array([nj.Qg[i].X for i in n.generators.keys()])
v = np.sqrt(np.array([n.v[i].X for i in n.v.keys()]))
Pg0 = data['gen']['Pg'].values
Qg0 = data['gen']['Qg'].values
v0 = data['bus']['Vm'].values

distj = np.linalg.norm(Pg - Pg0) #+ np.linalg.norm(Qg - Qg0) #+ np.linalg.norm(v - v0)
# %%
np.save("118Pg_jabr00.npy", Pg)
np.save("118Qg_jabr00.npy", Qg)

# %%
