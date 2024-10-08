#%%
import gurobipy as gp
from gurobi_optimods import datasets, opf
import numpy as np
from gurobipy import GRB
import networkx as nx
from itertools import chain, combinations

#%% 

class network:

    def __init__(self):
        pass

    def load_case(self, case):
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
        c = m.addVars([i for i in range(len(branches))], lb=0, ub=[V_max[i]*V_max[j] for (i, j) in branches],
                    vtype=GRB.CONTINUOUS, name=["c_{}{}{}".format(i, j, k) for k, (i, j) in enumerate(branches)])
        s = m.addVars([i for i in range(len(branches))], lb=[-V_max[i]*V_max[j] for (i, j) in branches], ub=[V_max[i]*V_max[j] for (i, j) in branches],
                    vtype=GRB.CONTINUOUS, name=["s_{}{}{}".format(i, j, k) for k, (i, j) in enumerate(branches)])
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

    def init_bouquet_model(self):
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
 
        s_abs = m.addVars([i for i in range(len(branches))], lb = 0, ub = 1,
                        vtype=GRB.CONTINUOUS, name=["s_abs_{}{}".format(i,j) for (i,j) in branches])
        u = m.addVars([i for i in range(len(branches))], vtype=GRB.BINARY, name=["u_{}{}".format(i,j) for (i,j) in branches])
        z_C = m.addVars(np.arange(len(loops)), lb=0, ub=1, vtype=GRB.CONTINUOUS, name = ["z_cycle_{}".format(i) for i in np.arange(len(loops))])
        #for loop in loops:
        As_list = [[A for A in chain.from_iterable(combinations(loop, r) for r in range(0, len(loop)+1, 2))] for loop in loops_b]
        z_A_index = [(ii, this_A) for ii in range(len(As_list)) for this_A in As_list[ii]]
        z_A = m.addVars(z_A_index, lb=0, ub=1, vtype = GRB.CONTINUOUS, name = ["z_{}".format(A) for A in z_A_index])
        lambda_A = m.addVars(z_A_index, vtype=GRB.BINARY, name = ["lambda_{}".format(A) for A in z_A_index])
        m_A = m.addVars(z_A_index, vtype=GRB.INTEGER, lb = 0, ub = [len(A[1]) // 2 for A in z_A_index], name = ["m_{}".format(A) for A in z_A_index])  
        # constraints
        for k, (i, j) in enumerate(branches):
            m.addConstr(c[k] ** 2 + s[k] ** 2 <= v[i] * v[j])  # (2)
            m.addConstr(s[k] == V_max[i] * V_max[j] * (2 * u[k] - 1) * s_abs[k]) # 

        def StdFormRelx(z, x, indices):
            for index in indices:
                m.addConstr(z <= x[index])
            m.addConstr(z + sum(1 - x[index] for index in indices) >= 1)
        '''
        for i,A in z_A_index:
                # TODO qui non so come chiamare s_abs perché non so che archi vengono pescati quando si fa il loop
                m.addLConstr(lambda_A[(i,A)] + 2*m_A[(i,A)] == sum(u[h] for h in A))
                print(i, A, loops_b[i])
                monomials = [s_abs[h] for h in A] + [c[k] / (V_max[k[0]] * V_max[k[1]]) for k in loops_b[i] if k not in A]
                StdFormRelx(z_A[(i,A)], monomials , np.arange(len(monomials)))
                '''
                
            

        #U_C = np.prod(V_max[i]**2 for i in loop_v)
        #loop constraint
        for i, As in enumerate(As_list):
            m.addConstr(sum((-1)**(len(A) // 2) * (1 - 2*lambda_A[(i,A)]) * z_A[(i,A)] for A in As) == z_C[i])



        for i, loop in enumerate(loops):
            z = z_C[i]
            x = v
            for index in loop:
                m.addConstr(z <= x[index] / (V_max[index]**2))
            m.addConstr(z+ sum(1 - x[index] / (V_max[index]**2) for index in loop) >= 1)


        #m.addLConstr(v[1] == 1)
        baseMVA = case["baseMVA"]
        for k, branch in enumerate(case['branch']):
            i = branch["fbus"]
            j = branch["tbus"]
            G, B = network.create_admittance_matrix(branch)
            # print("MIAO")
            m.addLConstr(P[k] / (baseMVA) == G[(i, i)] * v[i] + G[(i, j)] * c[k] + B[(i, j)] * s[k])   # (3)
            m.addLConstr(Q[k] / (baseMVA) == -B[(i, i)] * v[i] - B[(i, j)] * c[k] + G[(i, j)] * s[k])  # (4)
            m.addLConstr(P[k + len(branches)] / (baseMVA) == G[(j, j)] * v[j] + G[(j, i)] * c[k] - B[(j, i)] * s[k])   # (3)
            m.addLConstr(Q[k + len(branches)] / (baseMVA) == -B[(j, j)] * v[j] - B[(j, i)] * c[k] - G[(j, i)] * s[k])  # (4)

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

        # TODO objective function
        totalcost = 0
        for index, gen in enumerate(case["gen"]):
            gencost = case["gencost"][index]
            if gencost["costtype"] == 2: #polynomial cost
                costvector = gencost["costvector"]
                if len(costvector) == 3:
                    totalcost += costvector[0] * Pg[index] ** 2 + costvector[1] * Pg[index] + costvector[2]
                    
        m.setObjective(totalcost)

        
        
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

        for k, (i, j) in enumerate(branches):
            m.addConstr(c[k] ** 2 + s[k] ** 2 <= v[i] * v[j])  # (2)
            # m.addConstr(s[(i, j)] == V_max[i] * V_max[j] * (2 * u[(i,j)] - 1) * s_abs[(i,j)]) # 

        #m.addLConstr(v[1] == 1)
        baseMVA = case["baseMVA"]
        for k, branch in enumerate(case['branch']):
            i = branch["fbus"]
            j = branch["tbus"]
            G, B = network.create_admittance_matrix(branch)
            # print("MIAO")
            m.addLConstr(P[k] / (baseMVA) == G[(i, i)] * v[i] + G[(i, j)] * c[k] + B[(i, j)] * s[k])   # (3)
            m.addLConstr(Q[k] / (baseMVA) == -B[(i, i)] * v[i] - B[(i, j)] * c[k] + G[(i, j)] * s[k])  # (4)
            m.addLConstr(P[k + len(branches)] / (baseMVA) == G[(j, j)] * v[j] + G[(j, i)] * c[k] - B[(j, i)] * s[k])   # (3)
            m.addLConstr(Q[k + len(branches)] / (baseMVA) == -B[(j, j)] * v[j] - B[(j, i)] * c[k] - G[(j, i)] * s[k])  # (4)

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


# %%

case = datasets.load_opf_example("case118")
n = network()
n.load_case(case)

# %%
n.init_jabr_model()
n.m.optimize()

#%%
# n.init_bouquet_model()
# n.m.optimize()