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
import ObjBound

def Solve(case,plot,theta_lim):
    duration=time.time()

    Bus_Gen={}
    for i in range(len(case["gen"])):
        bus = case["gen"][i]["bus"]
        if Bus_Gen.get(bus)==None:
            Bus_Gen[bus] = [i]
        else:
            if i in Bus_Gen[bus]:
                raise ValueError
            else:
                Bus_Gen[bus].append(i)

    N=len(case['branch'])
    print(f"Number of branches : {N}")

    Branch_Bus={}
    for i in range(len(case["branch"])):
        fbus_ = case["branch"][i]["fbus"]
        tbus_ = case["branch"][i]["tbus"]
        if fbus_ not in Branch_Bus:
            Branch_Bus[fbus_] = [i] # Branch i starts from fbus_
        else:
            if i not in Branch_Bus[fbus_]:
                Branch_Bus[fbus_].append(i) # Branch i starts from fbus_
            else:
                # print(f"{fbus_,tbus_} double")
                raise ValueError
        if tbus_ not in Branch_Bus:
            Branch_Bus[tbus_] = [i+N] # Branch i+N starts from tbus_
        else:
            if i+N not in Branch_Bus[tbus_]:
                Branch_Bus[tbus_].append(i+N) # Branch i+N starts from tbus_
            else:
                # print(f"{tbus_,fbus_} double")
                raise ValueError

    Bus_SH={} #Dict that makes a bus_id refers to its gs and bs
    for bus in case["bus"]:
        Bus_SH[bus["bus_i"]]=[bus["Gs"],bus["Bs"]]

    def create_admittance_matrix(branch):
        G={}
        B={}
        #for branch in case["branch"]:
        fbus = branch["fbus"]
        tbus = branch["tbus"]
        # if branch==case['branch'][1]:
        #     G[(fbus,fbus)] = 0
        #     G[(fbus,tbus)] =0
        #     G[(tbus,fbus)] = 0
        #     G[(tbus,tbus)] = 0

        #     B[(fbus,fbus)] = 0
        #     B[(fbus,tbus)] = 0
        #     B[(tbus,fbus)] = 0
        #     B[(tbus,tbus)] = 0
        #     return G,B
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
        Y11 = A/(t**2)
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
        if fbus==tbus:
            raise ValueError
        return G,B

    baseMVA=case["baseMVA"]

    # Create a new model
    model = gp.Model("Jabr-ACOPF")

    ### Variables definition & intervals

    Pg = []  # Active power generation
    Qg = []  # Reactive power generation

    for i in range(len(case['gen'])):
        Pg.append(model.addVar(lb=case['gen'][i]['Pmin'], ub=case['gen'][i]['Pmax'], name=f"Pg{i}"))
        Qg.append(model.addVar(lb=case['gen'][i]['Qmin'], ub=case['gen'][i]['Qmax'], name=f"Qg{i}"))

    print("Power Gen var generated")

    v = {}
    for i in range(len(case["bus"])):
        bus_i = int(case["bus"][i]["bus_i"])
        v[bus_i] = model.addVar(lb=case["bus"][i]["Vmin"], ub=case["bus"][i]["Vmax"], name=f"v {bus_i}")
        
    v2 = {}
    for i in range(len(case["bus"])):
        bus_i = int(case["bus"][i]["bus_i"])
        v2[bus_i] = model.addVar(lb=case["bus"][i]["Vmin"]**2, ub=case["bus"][i]["Vmax"]**2, name=f"v2 {bus_i}")

    theta={}
    for i in range(len(case["bus"])):
        bus_i = int(case["bus"][i]["bus_i"])
        theta[bus_i] = model.addVar(lb=-theta_lim, ub=theta_lim, name=f"theta {bus_i}")

    print("v2 var generated")

    VMAX={}
    for i in range(len(case["bus"])):
        VMAX[case["bus"][i]["bus_i"]]=case["bus"][i]["Vmax"]
    VMIN={}
    for i in range(len(case["bus"])):
        VMIN[case["bus"][i]["bus_i"]]=case["bus"][i]["Vmin"]

    P={}
    Q={}
    c={}
    s={}
    c_mult={}
    s_mult={}
    c_lim={}
    s_lim={}
    
    for i in range(len(case['branch'])):
        fbus=int(case['branch'][i]['fbus'])
        tbus=int(case['branch'][i]['tbus'])
        vk_max=VMAX[fbus]
        vm_max=VMAX[tbus]
        vk_min=VMIN[fbus]
        vm_min=VMIN[tbus]
        if vm_max*vk_max<=0:
            raise ValueError
        P[i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"P{(fbus,tbus,i)}") #P_km
        P[i+N]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"P{(tbus,fbus,i+N)}")
        Q[i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"Q{(fbus,tbus,i)}") #Q_km
        Q[i+N]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"Q{(tbus,fbus,i+N)}")
        c[i]=model.addVar(lb=math.cos(2*theta_lim),ub=1,name=f"c{(fbus,tbus,i)}")
        s[i]=model.addVar(lb=-math.sin(2*theta_lim),ub=math.sin(2*theta_lim),name=f"s{(fbus,tbus,i)}")
        c_lim[i]=[math.cos(2*theta_lim),1]
        s_lim[i]=[-math.sin(2*theta_lim),math.sin(2*theta_lim)]
        c_mult[i]=model.addVar(lb=vm_min*vk_min*math.cos(2*theta_lim),ub=vm_max*vk_max,name=f"c_mult{(fbus,tbus,i)}") #jabr c_km
        s_mult[i]=model.addVar(lb=-vm_max*vk_max*math.sin(2*theta_lim),ub=vm_max*vk_max*math.sin(2*theta_lim),name=f"s_mult{(fbus,tbus,i)}") #jabr s_km

    nb_step=20
    A={}
    for i in tqdm(range(len(case['branch'])),ncols=100,desc=f"Linear cos and sine"): ### Linear cos and sine
        fbus=int(case['branch'][i]['fbus'])
        tbus=int(case['branch'][i]['tbus'])
        gap=abs(math.sin(math.pi/10)-math.pi/10)
        model.addConstr(s[i]<=(theta[fbus]-theta[tbus])+gap)
        model.addConstr(s[i]>=(theta[fbus]-theta[tbus])-gap)
        for k in range(0,nb_step+1):
            x_bar=-2*theta_lim+k/nb_step*4*theta_lim
            a=-math.sin(x_bar)
            b=math.cos(x_bar)-x_bar*a
            model.addConstr(c[i]<=a*(theta[fbus]-theta[tbus])+b)
        # sum_A=0
        # nb_step=5
        # step=4*theta_lim/nb_step
        # for k in range(0,nb_step):
        #     A[(i,k)]=model.addVar(vtype=GRB.BINARY, name=f"A {i},{k}")
        #     sum_A+=A[(i,k)]
        #     x_bar=-2*theta_lim+k/nb_step*4*theta_lim
        #     a=(math.cos(x_bar+step)-math.cos(x_bar))/step
        #     b=math.cos(x_bar)-x_bar*a
        #     print(a,b)
        #     model.addConstr(c[i]>=a*(theta[fbus]-theta[tbus])+b-(1-A[(i,k)]),name=f"c_binary {i,k}")
        # model.addConstr(sum_A==1)




    v_lim={}
    nb_step=5
    for i in tqdm(range(len(case["bus"])),ncols=100,desc="Linear v->v²"):
        fbus=case["bus"][i]["bus_i"]
        vk_max=VMAX[fbus]
        vk_min=VMIN[fbus]
        v_lim[fbus]=[vk_min,vk_max]
        a=(vk_max**2-vk_min**2)/(vk_max-vk_min)
        b=-vk_min*a+vk_min**2
        model.addConstr(a*v[fbus]+b >= v2[fbus],name=f"v2 {fbus}")
        step=(v_lim[fbus][1]-v_lim[fbus][0])/nb_step
        for k in range(nb_step):
            x_bar=v_lim[fbus][0]+step*k
            a=2*x_bar
            b=-a*x_bar+x_bar**2
            model.addConstr(a*v[fbus]+b <= v2[fbus],name=f"v2 {fbus}")

    for i in tqdm(range(len(case['branch'])),ncols=100,desc="Linear Mult cv² / sv²"): ### Linear mult
        fbus=int(case['branch'][i]['fbus'])
        tbus=int(case['branch'][i]['tbus'])
        vars=[c[i],v[fbus],v[tbus]]
        cuboid=[c_lim[i],v_lim[fbus],v_lim[tbus]]
        HP=hyperplanes.generate_all_neighbor_hp(cuboid)
        for hp in HP : #COMPUTE HYPERPLANES INEQUALITIES
            hyperplane=hp[0][:-1]
            constant=hp[0][-1]
            constraint=np.dot(hyperplane,vars)+constant
            if hp[1]>0:
                model.addConstr(constraint >= c_mult[i],name=f"c_mult {fbus} {tbus}")
                # print(f"HP, {loop_number,k,id_comb}")
            elif hp[1]<0:
                model.addConstr(constraint <= c_mult[i],name=f"c_mult {fbus} {tbus}")
                # print(f"HP, {loop_number,k,id_comb}")
            elif hp[1]==0:
                model.addConstr(constraint == c_mult[i],name=f"c_mult {fbus} {tbus}")
        vars=[s[i],v[fbus],v[tbus]]
        cuboid=[s_lim[i],v_lim[fbus],v_lim[tbus]]
        HP=hyperplanes.generate_all_neighbor_hp(cuboid)
        for hp in HP : #COMPUTE HYPERPLANES INEQUALITIES
            hyperplane=hp[0][:-1]
            constant=hp[0][-1]
            constraint=np.dot(hyperplane,vars)+constant
            if hp[1]>0:
                model.addConstr(constraint >= s_mult[i],name=f"s_mult {fbus} {tbus}")
                # print(f"HP, {loop_number,k,id_comb}")
            elif hp[1]<0:
                model.addConstr(constraint <= s_mult[i],name=f"s_mult {fbus} {tbus}")
                # print(f"HP, {loop_number,k,id_comb}")
            elif hp[1]==0:
                model.addConstr(constraint == s_mult[i],name=f"s_mult {fbus} {tbus}")



    print("Branch var generated")

    ### Constraints 
    # (1b) & (1c)

    for bus in tqdm(range(len(case['bus'])), desc="Const. (1b) & (1c) ...", ncols=100):
        bus_number = int(case['bus'][bus]['bus_i'])
        Pd=case['bus'][bus]['Pd']  # Active load
        Qd=case['bus'][bus]['Qd']  # Reactive load
        GenBus=Bus_Gen.get(bus_number)
        Pg_bus=0
        Qg_bus=0
        if GenBus!=None:
            Pg_bus=sum(Pg[i] for i in GenBus) # Active power gen on bus k
            Qg_bus=sum(Qg[i] for i in GenBus) # Reactive power gen on bus k
        P_bus=0
        Q_bus=0
        BranchBus=Branch_Bus.get(bus_number)
        if BranchBus!=None:
            P_bus=sum(P[i] for i in BranchBus) 
            Q_bus=sum(Q[i] for i in BranchBus) 
        # for i,branch in enumerate(case["branch"]): Conventional way
        #     if branch["fbus"]==bus_number:
        #         P_bus+=P[i]
        #         Q_bus+=Q[i]
        #     if branch["tbus"]==bus_number:
        #         P_bus+=P[i+N]
        #         Q_bus+=Q[i+N]
        model.addLConstr(P_bus == Pg_bus - Pd , name=f"ActivePowerBalance_{bus_number}")
        model.addLConstr(Q_bus == Qg_bus - Qd , name=f"ReactivePowerBalance_{bus_number}")

    print("Const. (1b) & (1c) generated")

    # (6) Jabr const.

    for i in tqdm(range(len(case['branch'])), desc="Const. (6) Jabr ..." ,ncols=100):
        
        branch_=case["branch"][i]
        G,B=create_admittance_matrix(branch_)
        k=branch_['fbus']
        m=branch_['tbus']

        model.addLConstr(P[i]/baseMVA==G[(k,k)]*v2[k]+G[(k,m)]*c_mult[i]+B[(k,m)]*s_mult[i], name=f"6a_{k},{m},{i}")
        model.addLConstr(P[i+N]/baseMVA==G[(m,m)]*v2[m]+G[(m,k)]*c_mult[i]-B[(m,k)]*s_mult[i], name=f"6b_{m},{k},{i}")
        model.addLConstr(Q[i]/baseMVA==-B[(k,k)]*v2[k]-B[(k,m)]*c_mult[i]+G[(k,m)]*s_mult[i], name=f"6c_{k},{m},{i}")
        model.addLConstr(Q[i+N]/baseMVA==-B[(m,m)]*v2[m]-B[(m,k)]*c_mult[i]-G[(m,k)]*s_mult[i], name=f"6d_{m},{k},{i}")

    # for i in range(len(case["branch"])): #SOCP Jabr
    #     fbus=case["branch"][i]["fbus"]
    #     tbus=case["branch"][i]["tbus"]
    #     model.addConstr(c[i] ** 2 + s[i] ** 2 <= v2[fbus] * v2[tbus]) # (2)

    print("Const. (6) Jabr generated")

    alpha_tot=0
    beta_tot=0
    gamma_tot=0
    Ki_tot=0
    for i in range(len(case['branch'])):
        branch_=case["branch"][i]
        rateA=case["branch"][i]["rateA"]
        r = branch_["r"] # real part of impendence
        x = branch_["x"] # complex part of impendence
        z = r + x*1j # impendence
        y = 1/z
        g=y.real
        b=y.imag

        fbus=case['branch'][i]["fbus"]
        tbus=case['branch'][i]["tbus"]

        ratio = case["branch"][i]["ratio"]
        bshunt = case["branch"][i]["b"]
        angle = case["branch"][i]["angle"]

        if ratio==0:
            ratio=1

        # mp_cbusf=v2[fbus]
        # mp_cbust=v2[tbus]
        # mp_c=c[i]
        # mp_s=s[i]

        # i2_var  = 0
        # i2_var += (g*g + b*b)/(ratio*ratio) * ( (mp_cbusf/(ratio*ratio)) + mp_cbust - (2/ratio) * ( mp_c * math.cos(angle) + mp_s * math.sin(angle) ) )
        # i2_var += b*bshunt/(ratio**3) * ( (mp_cbusf/ratio) - (mp_c * math.cos(angle) + mp_s * math.sin(angle) ))
        # i2_var += g*bshunt/(ratio**3) * ( mp_s * math.cos(angle) - mp_c * math.sin(angle) )
        # i2_var += (bshunt*bshunt*mp_cbusf/(4*(ratio**4)) )

        tau=ratio
        sigma=angle

        gsh=0
        bsh=bshunt

        Alpha=(1/(tau**4))*((g**2+b**2)+(g*gsh+b*bsh)+((gsh**2+bsh**2)/4))
        alpha_tot+=Alpha
        Beta=(g**2+b**2)/(tau**2)
        beta_tot+=Beta/Alpha
        Gamma=(1/(tau**3))*(math.cos(sigma)*(-2*(g**2+b**2)-(g*gsh+b*bsh))+math.sin(sigma)*(b*gsh-g*bsh))
        gamma_tot+=Gamma/Alpha
        Ki=(1/(tau**3))*(math.sin(sigma)*(-2*(g**2+b**2)-(g*gsh+b*bsh))-math.cos(sigma)*(b*gsh-g*bsh))
        Ki_tot+=Ki/Alpha

        i2like=Alpha*v2[fbus]+Beta*v2[tbus]+Gamma*c_mult[i]+Ki*s_mult[i]

        rho=100
        
        i2={}
        if Alpha<rho:
            i2[i]=model.addVar(lb=0,ub=GRB.INFINITY,name=f"i2{(fbus,tbus,i)}") #i2km
            # model.addConstr(i2_var==i2[i],name=f"goodi2 {fbus},{tbus},{i}")
            model.addConstr(i2like==i2[i],name=f"goodi2 {fbus},{tbus},{i}") #TRUE

        else:
            # i2_var=i2_var/Alpha
            # model.addConstr(i2_var<=rateA/Alpha,name=f"badi2-1 {fbus},{tbus},{i}")
            # model.addConstr(i2_var>=0,name=f"badi2-2 {fbus},{tbus},{i}")
            i2like=i2like/Alpha
            model.addConstr(i2like<=rateA/Alpha,name=f"badi2-1 {fbus},{tbus},{i}")
            model.addConstr(i2like>=0,name=f"badi2-2 {fbus},{tbus},{i}") #TRUE

    # print("ALPHA MEAN "+str(alpha_tot/len(case["branch"])))
    # print("BETA/ALPHA MEAN "+str(beta_tot/len(case["branch"])))
    # print("GAMMA/ALPHA MEAN "+str(gamma_tot/len(case["branch"])))
    # print("KI/ALPHA MEAN "+str(Ki_tot/len(case["branch"])))

    print("Const. i2 generated")

    # for i in range(len(case['branch'])): #SCOP limit
    #     rateA=case['branch'][i]["rateA"]
    #     if rateA!=0:
    #         model.addConstr(P[i] ** 2 + Q[i] ** 2 <= rateA**2/(baseMVA**2),name=f"limit {i},{rateA}")
    #         model.addConstr(P[i+N] ** 2 + Q[i+N] ** 2 <= rateA**2/(baseMVA**2),name=f"limit {i+N},{rateA}")

    print("Const. limit generated")

    # model.setParam('OutputFlag', 0)
    # cs_bounds=ObjBound.GetBounds(model,s,c,v2)
    # model.setParam('OutputFlag', 1)

    # cs_bounds=[[],[],[],[]] #Activate to delete variable bounding
    # for i in range(len(case['branch'])):
    #     cs_bounds[0].append(s[i].LB)
    #     cs_bounds[1].append(s[i].UB)
    #     cs_bounds[2].append(c[i].LB)
    #     cs_bounds[3].append(c[i].UB)

    # for i in range(len(case["branch"])):
    #     fbus=case["branch"][i]['fbus']
    #     tbus=case["branch"][i]['tbus']
    #     model.addConstr(cs_bounds[0][i]<=s[i],name=f"LB s{fbus,tbus,i}")
    #     model.addConstr(cs_bounds[1][i]>=s[i],name=f"LB s{fbus,tbus,i}")
    #     model.addConstr(cs_bounds[2][i]<=c[i],name=f"LB c{fbus,tbus,i}")
    #     model.addConstr(cs_bounds[3][i]>=c[i],name=f"LB c{fbus,tbus,i}")

    def score(s_lim,c_lim,v_lim=[]):
        score_=0
        for j in s_lim:
            i=s_lim[j]
            score_+=i[1]-i[0]
        for j in c_lim:
            i=c_lim[j]
            score_+=i[1]-i[0]
        for j in v_lim:
            i=v_lim[j]
            score_+=i[1]-i[0]
        return score_
    
    eps_bound=1
    score_=score(s_lim,c_lim,v_lim)
    print(f"Score = {score_:.4f}")
    Bound_constr=[]
    Bounds_save=[]
    while eps_bound>0.01:
    # for _ in range(0,1000):


        model.setParam('OutputFlag', 0)
        # new_bounds=ObjBound.GetBoundsTheta(model,case,s,c)
        new_bounds=ObjBound.GetBoundsTheta(model,case,s,c,v)
        Bounds_save.append(new_bounds)
        model.setParam('OutputFlag', 1)
        mem_s_lim=s_lim.copy()
        mem_c_lim=c_lim.copy()
        mem_v_lim=v_lim.copy()
        s_lim=new_bounds[0]
        c_lim=new_bounds[1]
        if len(new_bounds)>=3:
            v_lim=new_bounds[2]
        mem_score=score_
        score_=score(s_lim,c_lim,v_lim)
        eps_bound=(mem_score-score_)/mem_score
        print(f"Score = {score_:.4f}")
        print(f"Eps = {eps_bound:.4f}")

        ###DELETE ALL HP
        # for constr in Bound_constr:
        #     model.remove(constr)

        nb_step=5
        for i in range(len(case["bus"])):
            fbus=case["bus"][i]["bus_i"]
            a=(vk_max+vk_min)
            b=-vk_min*a+vk_min**2
            Bound_constr.append(model.addConstr(a*v[fbus]+b >= v2[fbus],name=f"v2 {fbus}"))
            step=(v_lim[fbus][1]-v_lim[fbus][0])/nb_step
            for k in range(nb_step):
                x_bar=v_lim[fbus][0]+step*k
                a=2*x_bar
                b=-a*x_bar+x_bar**2
                Bound_constr.append(model.addConstr(a*v[fbus]+b <= v2[fbus],name=f"v2 {fbus}"))

        for i in tqdm(range(len(case['branch'])),ncols=100,desc="Linear Mult cv² / sv²"): ### Linear mult
            fbus=int(case['branch'][i]['fbus'])
            tbus=int(case['branch'][i]['tbus'])
            vars=[c[i],v[fbus],v[tbus]]
            cuboid=[c_lim[i],v_lim[fbus],v_lim[tbus]]
            HP=hyperplanes.generate_all_neighbor_hp(cuboid)
            for hp in HP : #COMPUTE HYPERPLANES INEQUALITIES
                hyperplane=hp[0][:-1]
                constant=hp[0][-1]
                constraint=np.dot(hyperplane,vars)+constant
                if hp[1]>0:
                    Bound_constr.append(model.addConstr(constraint >= c_mult[i],name=f"c_mult {fbus} {tbus}"))
                    # print(f"HP, {loop_number,k,id_comb}")
                elif hp[1]<0:
                    Bound_constr.append(model.addConstr(constraint <= c_mult[i],name=f"c_mult {fbus} {tbus}"))
                    # print(f"HP, {loop_number,k,id_comb}")
                elif hp[1]==0:
                    Bound_constr.append(model.addConstr(constraint == c_mult[i],name=f"c_mult {fbus} {tbus}"))
            vars=[s[i],v[fbus],v[tbus]]
            cuboid=[s_lim[i],v_lim[fbus],v_lim[tbus]]
            HP=hyperplanes.generate_all_neighbor_hp(cuboid)
            for hp in HP : #COMPUTE HYPERPLANES INEQUALITIES
                hyperplane=hp[0][:-1]
                constant=hp[0][-1]
                constraint=np.dot(hyperplane,vars)+constant
                if hp[1]>0:
                    Bound_constr.append(model.addConstr(constraint >= s_mult[i],name=f"s_mult {fbus} {tbus}"))
                    # print(f"HP, {loop_number,k,id_comb}")
                elif hp[1]<0:
                    Bound_constr.append(model.addConstr(constraint <= s_mult[i],name=f"s_mult {fbus} {tbus}"))
                    # print(f"HP, {loop_number,k,id_comb}")
                elif hp[1]==0:
                    Bound_constr.append(model.addConstr(constraint == s_mult[i],name=f"s_mult {fbus} {tbus}"))


    objective = gp.QuadExpr()

    for i in range(len(case['gen'])):
        ca,cb,cc = case['gencost'][i]['costvector']
        objective += ca*Pg[i]*Pg[i]+cb*Pg[i]+cc

    model.setObjective(objective, GRB.MINIMIZE)

    model.optimize()

    if model.status == GRB.INFEASIBLE:
        print("Model is infeasible.")
        # You can try to find infeasibilities by calling the following:
        model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
        model.write("infeasible_model.ilp")

    print(f'Temps : {time.time()-duration}')

    Cuts=[]

    def JabrCuts(eps,pT,Const_to_add): 
        violation_amount=[]
        for i in range(len(case['branch'])):
            fbus=case['branch'][i]["fbus"]
            tbus=case['branch'][i]["tbus"]
            amount=c_mult[i].X*c_mult[i].X+s_mult[i].X*s_mult[i].X-v2[fbus].X*v2[tbus].X
            if amount>eps:
                violation_amount.append([amount,i,c_mult[i].X,s_mult[i].X,v2[fbus].X,v2[tbus].X])
        violation_amount.sort(reverse=True) # Order the violation amounts
        violation_amount=violation_amount[:round(len(violation_amount)*pT)] # Keep the top p_T percent
    # for i in range(len(violation_amount)):
        for i in range(len(violation_amount)):
            branch_i=violation_amount[i][1]
            fbus=case['branch'][branch_i]["fbus"]
            tbus=case['branch'][branch_i]["tbus"]
            c_=violation_amount[i][2]
            s_=violation_amount[i][3]
            v2k=violation_amount[i][4]
            v2m=violation_amount[i][5]
            if v2k+v2m>eps:
                n0=((2*c_)**2+(2*s_)**2+(v2k-v2m)**2)**(1/2)
                Const_to_add.append([4*c_*c_mult[branch_i]+4*s_*s_mult[branch_i]+(v2k-v2m-n0)*v2[fbus]+(v2m-v2k-n0)*v2[tbus],f"CutJabr_{fbus},{tbus},{branch_i}"])
        return len(violation_amount)

    def f_eps_par(new_cut, Cuts):
        dict_new_cut={}
        for term in range(new_cut.size()):
            dict_new_cut[new_cut.getVar(term)]=new_cut.getCoeff(term)
        new_cut_norm=0
        for i in dict_new_cut:
            new_cut_norm+=dict_new_cut[i]**2
        new_cut_norm=new_cut_norm**(1/2)

        for cuts in Cuts:
            dot_product=0
            cut=cuts[0]
            row = model.getRow(cut)
            for k in range(row.size()):
                if row.getVar(k) in dict_new_cut:
                    dot_product+=dict_new_cut[row.getVar(k)]*row.getCoeff(k)
            if dot_product!=0:
                cut_norm=0
                for k in range(row.size()):
                    cut_norm+=row.getCoeff(k)**2
                cut_norm=cut_norm**(1/2)

                cos_theta = dot_product / (new_cut_norm * cut_norm)

                if cos_theta > 1 - eps_par:
                    return False 
        return True

    T=0
    Tlim=1000
    pjabr=0.55
    pi2=0.15
    plim=1
    Tage=5
    Tftol=5
    eps=10E-5
    eps_par=10E-8
    eps_ftol=10E-8
    eps_slack=10E-8
    r=0


    nb_const=[[0]]
    nb_const_tot=[0]
    nb_const_deleted=[0]
    obj=model.getObjective()
    val=[obj.getValue()]
    temps=[time.time()-duration]

    initval=val[0]
    ftol=True
    init=True #I am still on the root ?
    while(T<Tlim and (r<Tftol or init)):

        nb_const.append([])
        T+=1
        print(f"{T}--------------------CUTS-------------------")
        
        to_remove=[]
        for i,cut in enumerate(Cuts): #Cut old ones that aren't slack
            constraint = cut[0]
            slack=-model.getRow(constraint).getValue()
            if slack > eps_slack and (T - cut[1]) >= Tage and not init: # Removing condition
                to_remove.append([constraint,cut])
        remove_compt=0
        for rr in range(len(to_remove)):
            remove=to_remove[rr]
            model.remove(remove[0])  # Remove from the model
            Cuts.remove(remove[1])  # Remove from Cuts
            remove_compt+=1
        nb_const_deleted.append(remove_compt)
        
        Const_to_add=[]
        nb_const[-1].append(JabrCuts(eps,pjabr,Const_to_add))
        # nb_const[-1].append(i2Cuts(eps,pi2,Const_to_add))
        # nb_const[-1].append(LimitCuts(eps,plim,Const_to_add))
        print(nb_const[-1])
        nb_const_tot.append(sum(nb_const[-1]))
        if(nb_const_tot[-1]==0):
            nb_const_tot.remove(nb_const_tot[-1])
            nb_const_deleted.remove(nb_const_deleted[-1])
            break
        cuts_to_add=[True for i in range(len(Const_to_add))]
        # for i,new_cuts in enumerate(Const_to_add): # Test de parallélisme, enelever la boucle.
        #     if not f_eps_par(new_cuts[0],Cuts):
        #         cuts_to_add[i]=False
        #         nb_const_tot[-1]-=1
        for i in range(len(cuts_to_add)):
            if cuts_to_add[i]:
                Cuts.append([model.addLConstr(Const_to_add[i][0] <= 0, name=Const_to_add[i][1]),T])

        model.update()  # Synchronize the model after removal
        
        model.optimize()

        temps.append(time.time()-duration)
        if model.status == GRB.INFEASIBLE:
            print("Model is infeasible.")
            model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
            model.write("infeasible_model.ilp")
            break
        obj=model.getObjective()
        val.append(obj.getValue())
        if(val[-1]-val[-2]<val[-2]*eps_ftol):
            r+=1
        else:
            r=0
        if val[-1]-initval>10**-5:
            init=False
        
        
    print(val)
    print(nb_const_tot)
    print(nb_const)

    fig, ax1 = plt.subplots()

    ax1.plot(temps,val, color='blue', label='Lower Bound')
    ax1.scatter(temps, val, color='blue', marker='|')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Lower Bound', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')

    ax2 = ax1.twinx()
    ax2.plot(temps,nb_const_tot, color='red', label='#Added Constraints')
    ax2.scatter(temps, nb_const_tot, color='red', marker='|')
    ax2.plot(temps,nb_const_deleted, color='green', label='#Deleted Constraints')
    ax2.scatter(temps, nb_const_deleted, color='green', marker='|')
    ax2.set_ylabel('#Constraints', color='black')
    ax2.tick_params(axis='y', labelcolor='black')

    ax2.legend(loc='lower right')

    last_val = val[-1]
    ax1.text(0.95, 0.90, f'{last_val:.5e}', transform=ax1.transAxes, ha='right', va='top', color='blue')

    plt.title(f'Instance : {case["casename"]}')

    print("Jabr = violation")

    for i in range(len(case["branch"])):
        branch=case["branch"][i]
        fbus=branch["fbus"]
        tbus=branch["tbus"]
        violation=abs(c_mult[i].X**2+s_mult[i].X**2-v2[fbus].X*v2[tbus].X)
        if violation<=1e-04:
            print(f"{i} {violation:.4e}")
        else:
            print(f"{i} {violation:.4e} WARNING")


    # for bus in range(len(case['bus'])):
    #     bus_number = int(case['bus'][bus]['bus_i'])
    #     Pd=case['bus'][bus]['Pd']  # Active load
    #     Qd=case['bus'][bus]['Qd']  # Reactive load
    #     GenBus=Bus_Gen.get(bus_number)
    #     Pg_bus=0
    #     Qg_bus=0
    #     if GenBus!=None:
    #         Pg_bus=sum(Pg[i].X for i in GenBus) # Active power gen on bus k
    #         Qg_bus=sum(Qg[i].X for i in GenBus) # Reactive power gen on bus k
    #     P_bus=0
    #     Q_bus=0
    #     BranchBus=Branch_Bus.get(bus_number)
    #     if BranchBus!=None:
    #         P_bus=sum(P[i].X for i in BranchBus) 
    #         Q_bus=sum(Q[i].X for i in BranchBus) 
        # print(f"Violation P : {abs(P_bus-(Pg_bus - Pd)):.3e}")
        # print(f"Violation P : {abs(Q_bus-(Qg_bus - Qd)):.3e}")

    # for i in range(len(case["branch"])):
    #     fbus=case["branch"][i]["fbus"]
    #     tbus=case["branch"][i]["tbus"]
    #     print(f"{fbus},{tbus} = {s[i].X}")

 
    ssinmin=5
    ssinmax=-5
    ccosmin=5
    ccosmax=-5
    for i in range(len(case["branch"])):
        fbus=case["branch"][i]["fbus"]
        tbus=case["branch"][i]["tbus"]
        qtte=s[i].X-math.sin(theta[fbus].X-theta[tbus].X)
        ssinmin=min(ssinmin,qtte)
        ssinmax=max(ssinmax,qtte)
        qtte=c[i].X-math.cos(theta[fbus].X-theta[tbus].X)
        ccosmin=min(ccosmin,qtte)
        ccosmax=max(ccosmax,qtte)
        # print(f"s-sin {fbus} : {s_mult[fbus].X-math.sin(theta[fbus].X)}")
        # print(f"c-cos {fbus} : {c_mult[fbus].X-math.cos(theta[fbus].X)}")
    print(f"s-sin : [{ssinmin:.4f},{ssinmax:.4f}]")
    print(f"c-cos : [{ccosmin:.4f},{ccosmax:.4f}]")

    vvmax=v2[case["bus"][0]["bus_i"]].X-v[case["bus"][0]["bus_i"]].X
    vvmin=v2[case["bus"][0]["bus_i"]].X-v[case["bus"][0]["bus_i"]].X
    for i in range(len(case["bus"])):
        fbus=case["bus"][i]["bus_i"]
        qtte=v2[fbus].X-v[fbus].X**2
        vvmax=max(qtte,vvmax)
        vvmin=min(qtte,vvmin)
        # print(f"v²-v {fbus} : {v2[fbus].X-v[fbus].X}")
    print(f"v2-v² : [{vvmin:.4f},{vvmax:.4f}]")

    ssv2min=c_mult[i].X-c[i].X*v[case["branch"][0]["fbus"]].X*v[case["branch"][0]["tbus"]].X
    ssv2max=c_mult[i].X-c[i].X*v[case["branch"][0]["fbus"]].X*v[case["branch"][0]["tbus"]].X
    ccv2min=s_mult[i].X-s[i].X*v[case["branch"][0]["fbus"]].X*v[case["branch"][0]["tbus"]].X
    ccv2max=s_mult[i].X-s[i].X*v[case["branch"][0]["fbus"]].X*v[case["branch"][0]["tbus"]].X
    for i in range(len(case["branch"])):
        fbus=case["branch"][i]["fbus"]
        tbus=case["branch"][i]["tbus"]
        qtte=c_mult[i].X-c[i].X*v[fbus].X*v[tbus].X
        ccv2min=min(ccv2min,qtte)
        ccv2max=max(ccv2max,qtte)
        qtte=s_mult[i].X-s[i].X*v[fbus].X*v[tbus].X
        ssv2min=min(ssv2min,qtte)
        ssv2max=max(ssv2max,qtte)
        # print(f"c_mult-cv² {fbus} {tbus} : {c_mult[i].X-c[i].X*v[fbus].X*v[tbus].X}")
        # print(f"s_mult-sv² {fbus} {tbus} : {s_mult[i].X-s[i].X*v[fbus].X*v[tbus].X}")
    print(f"s_mult-sv² : [{ssv2min:.4f},{ssv2max:.4f}]")
    print(f"c_mult-cv² : [{ccv2min:.4f},{ccv2max:.4f}]")
    ### Computing violation quantity ###

    tot_val_viol=0
    for i in range(len(case['branch'])):
        
        branch_=case["branch"][i]
        fbus=branch_["fbus"]
        tbus=branch_["tbus"]
        G,B=create_admittance_matrix(branch_)
        k=branch_['fbus']
        m=branch_['tbus']
        aureol=c_mult[i].X**2+s_mult[i].X**2
        new_c=c_mult[i].X/math.sqrt(aureol)
        new_s=s_mult[i].X/math.sqrt(aureol)

        val_viol=0
        val_viol+=abs(P[i].X/baseMVA-(G[(k,k)]*v2[k].X+G[(k,m)]*new_c+B[(k,m)]*new_s))
        val_viol+=abs(P[i+N].X/baseMVA-(G[(m,m)]*v2[m].X+G[(m,k)]*new_c-B[(m,k)]*new_s))
        val_viol+=abs(Q[i].X/baseMVA-(-B[(k,k)]*v2[k].X-B[(k,m)]*new_c+G[(k,m)]*new_s))
        val_viol+=abs(Q[i+N].X/baseMVA-(-B[(m,m)]*v2[m].X-B[(m,k)]*new_c-G[(m,k)]*new_s))
        tot_val_viol+=val_viol
        # print(f"Branch {fbus},{tbus} = {val_viol:.3e}")
    print(f"TOT VIOLATION ={tot_val_viol}")

    if len(case["bus"])<=500:
        G = nx.Graph()  # building the support graph of the network
        for i in range(len(case["branch"])):
            fbus=case["branch"][i]["fbus"]
            tbus=case["branch"][i]["tbus"]
            G.add_edge(fbus,tbus)
        loops = nx.minimum_cycle_basis(G) 

        fbustbus_ID={}
        for i in range(len(case["branch"])):
            fbus=case["branch"][i]["fbus"]
            tbus=case["branch"][i]["tbus"]
            fbustbus_ID[(fbus,tbus)]=i
            fbustbus_ID[(tbus,fbus)]=i+len(case["branch"])

        sum_tot=0
        for loop in loops:
            sum_=0
            for i in range(len(loop)):
                index=fbustbus_ID[(loop[i],loop[(i+1)%len(loop)])]
                if index>=len(case["branch"]):
                    index=index-len(case["branch"])
                aureol=c_mult[index].X**2+s_mult[index].X**2
                new_c=c_mult[index].X/math.sqrt(aureol)
                new_s=s_mult[index].X/math.sqrt(aureol)
                diffTheta=math.asin(new_s)
                sum_+=diffTheta
            sum_tot+=sum_
        print(f"Contraintes de loop : {sum_tot} grad")

    for i in range(len(case["bus"])):
        fbus=case["bus"][i]["bus_i"]
        val_=theta[fbus].X
        val_2=v[fbus].X
        print(f"{fbus} : Theta = {val_:.4f}, Tension = {val_2:.4f}")

    for i in range(len(Bounds_save[0])):
        for j in range(len(Bounds_save)):
            print(f"{Bounds_save[j][i] }")

    if plot:
        plt.show()
    return [val[-1],tot_val_viol]

if __name__=="__main__": 
    local_path="C:/Users/ytsha/AppData/Local/Programs/Python/Python310/Lib/site-packages/gurobi_optimods/data/opf/"
    
    # name="case9"
    name="case14"
    name="case30"
    # name="case57"
    # name="case118"
    # name="caseNY"
    # name="case_ACTIVSg500"
    # name="case1354pegase"
    # name="case9241pegase"
    # name="case13659pegase"
    path=local_path+name+".m"
    case = conversion.conversion(path)
    theta_lim=math.pi/10
    print(Solve(case,True,theta_lim))