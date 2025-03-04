import re

def conversion(path):
    # Lire le fichier .m
    with open(path, 'r') as file:
        content = file.read()

    # Fonction pour extraire les données sous forme de dictionnaire
    def extract_matrix_to_dict(text, headers, costtype=False):
        matrix_pattern = re.compile(r'\[([^\]]+)\]', re.DOTALL)
        match = matrix_pattern.search(text)
        if match:
            matrix_text = match.group(1).strip()
            rows = matrix_text.split(';')
            data = []
            for row in rows:
                values = list(map(float, row.split()))
                if costtype:
                    if(len(values)>1):
                        # For gencost, format it in the specific dictionary
                        cost_dict = {
                            'costtype': values[0],  # First value is costtype
                            'startup': values[1],   # Second value is startup
                            'shutdown': values[2],  # Third value is shutdown
                            'n': values[3],         # Fourth value is n
                            'costvector': values[4:]  # Remaining values form costvector
                        }
                        data.append(cost_dict)
                else:
                    data.append(dict(zip(headers, values)))
            return data
        return None

    # Dictionnaire pour stocker toutes les données
    mpc = {}

    # Extraire le nom du cas (casename)
    casename_match = re.match(r'function\s+mpc\s*=\s*(\w+)', content)
    if casename_match:
        mpc['casename'] = casename_match.group(1)

    # Extraire la donnée des buses (bus data)
    bus_data_text = re.search(r'mpc\.bus\s*=\s*\[.*?\];', content, re.DOTALL)
    if bus_data_text:
        bus_headers = [
            'bus_i', 'type', 'Pd', 'Qd', 'Gs', 'Bs', 'area', 'Vm', 'Va', 'baseKV', 'zone', 'Vmax', 'Vmin'
        ]
        mpc['bus'] = extract_matrix_to_dict(bus_data_text.group(0), bus_headers)
    if(mpc['bus'][-1]=={}):
        mpc['bus']=mpc['bus'][:-1]


    # Extraire les données des générateurs (generator data)
    gen_data_text = re.search(r'mpc\.gen\s*=\s*\[.*?\];', content, re.DOTALL)
    if gen_data_text:
        gen_headers = [
            'bus', 'Pg', 'Qg', 'Qmax', 'Qmin', 'Vg', 'mBase', 'status', 'Pmax', 'Pmin', 'Pc1', 'Pc2', 'Qc1min', 
            'Qc1max', 'Qc2min', 'Qc2max', 'ramp_agc', 'ramp_10', 'ramp_30', 'ramp_q', 'apf'
        ]
        mpc['gen'] = extract_matrix_to_dict(gen_data_text.group(0), gen_headers)
    if(mpc['gen'][-1]=={}):
        mpc['gen']=mpc['gen'][:-1]

    # Extraire les données des branches (branch data)
    branch_data_text = re.search(r'mpc\.branch\s*=\s*\[.*?\];', content, re.DOTALL)
    if branch_data_text:
        branch_headers = [
            'fbus', 'tbus', 'r', 'x', 'b', 'rateA', 'rateB', 'rateC', 'ratio', 'angle', 'status', 
            'angmin', 'angmax'
        ]
        mpc['branch'] = extract_matrix_to_dict(branch_data_text.group(0), branch_headers)
    if(mpc['branch'][-1]=={}):
        mpc['branch']=mpc['branch'][:-1]
    for branch in range(len(mpc['branch'])):
        mpc['branch'][branch]["fbus"]=int(mpc['branch'][branch]["fbus"])
        mpc['branch'][branch]["tbus"]=int(mpc['branch'][branch]["tbus"])

    # Extraire la donnée de base MVA
    baseMVA_text = re.search(r'mpc\.baseMVA\s*=\s*(\d+);', content)
    if baseMVA_text:
        mpc['baseMVA'] = float(baseMVA_text.group(1))

    # Extraire la donnée de cost (gencost)
    gencost_text = re.search(r'mpc\.gencost\s*=\s*\[.*?\];', content, re.DOTALL)
    if gencost_text:
        # Extract gencost and format it
        mpc['gencost'] = extract_matrix_to_dict(gencost_text.group(0), None, costtype=True)
    if(mpc['gencost'][-1]=={}):
        mpc['gencost']=mpc['gencost'][:-1]

    # Afficher les données extraites
    return mpc

# Exemple d'utilisation
path="C:/Users/ytsha/AppData/Local/Programs/Python/Python310/Lib/site-packages/gurobi_optimods/data/opf/case_ACTIVSg25k.m"
case = conversion(path)

path="C:/Users/ytsha/AppData/Local/Programs/Python/Python310/Lib/site-packages/gurobi_optimods/data/opf/case_ACTIVSg25k.m"


import itertools
import numpy as np
from tqdm import tqdm

def calculateHyperplane(dim):
    coeffs=[]
    for i in range(2**dim):
        coeffs.append([])
        n=i
        while n != 0 or len(coeffs[-1]) != dim:
            if n % 2 == 0:
                coeffs[-1].append(-1)
            else:
                coeffs[-1].append(1)  # Utilisation de 1 au lieu de n%2
            n = n // 2
        coeffs[-1].append(1)

    print("Coefficients:", coeffs)
    
    # Génération des combinaisons de triplets
    triplets = list(itertools.combinations(coeffs, dim + 1))

    for i in range(len(triplets)):
        triplets[i] = np.array(triplets[i])
    valid=0
    tot=0
    compteur=0
    flat=0
    nb_neighbor=0
    # Résolution du système linéaire pour chaque triplet
    for i in tqdm(range(len(triplets)), desc='Nb N+1-plets', ncols=100):
        compteur+=1
        # On calcule b avec la multiplication de chaque élément du triplet
        b = np.array([mult(j) for j in triplets[i]])
        # print("")
        # print(f"Résolution pour le triplet {i}:")
        # print(f"Triplet: {triplets[i]}")
        # print(f"Valeur b: {b}")

        # Vérification de la forme de la matrice avant la résolution
        if triplets[i].shape[0] == triplets[i].shape[1]:
            try:
                solution = np.linalg.solve(triplets[i], b)
                # print("Solution:", solution)
                tot+=1
                if verification(solution,coeffs)==True:
                    # print("ALERTE")
                    # print(f"n+1-plet: {triplets[i]}")
                    # print(f"Valeur b: {b}")
                    # print("Solution:", solution)
                    valid+=1
                    b0=b[0]
                    b_flat=True
                    for t in range(len(b)):
                        if b[t]!=b0:
                            b_flat=False
                    if b_flat:
                        flat+=1
                    else:
                        print("ALERTE")
                        print(f"n+1-plet: {triplets[i]}")
                        print(f"Valeur b: {b}")
                        print("Solution:", solution)
                        if is_neighbor(triplets[i]):
                            nb_neighbor+=1
            except np.linalg.LinAlgError as e:
                print(f"Erreur lors de la résolution pour le triplet {triplets[i]} : {e}")
                continue
        # else:
            # print(f"Le triplet {i} n'est pas une matrice carrée. Impossible de résoudre.")
        # print(len(triplets)*tot/compteur)
    print("Nb sol = "+str(tot)+"/"+str(compteur))
    print("Nb valides = "+str(valid)+"/"+str(tot))
    print("Nb flat among valid = "+str(flat)+"/"+str(valid))
    print("Nb neighbors among valid non flat = "+str(nb_neighbor))

def is_neighbor(coeffs):
    coeffs=np.array(coeffs)
    distances=[]
    for init in coeffs:
        distances.append(0)
        for i in coeffs:
            distances[-1]=max(distances[-1],sum(abs(init-i)))
    if 2 in distances:
        return True
    return False

def voisin(x):
    L=[x]
    for k in range(len(x)-1):
        a=x.copy()
        a[k]=x[k]*(-1)
        L.append(a)
    return L

def NeighborHyperplane(dim):
    coeffs=[]
    voisins=[]
    for i in range(2**dim):
        coeffs.append([])
        n=i
        while n != 0 or len(coeffs[-1]) != dim:
            if n % 2 == 0:
                coeffs[-1].append(-1)
            else:
                coeffs[-1].append(1)  # Utilisation de 1 au lieu de n%2
            n = n // 2
        coeffs[-1].append(1)
        voisins.append(voisin(coeffs[-1]))
    # print("Coefficients:", coeffs)

    tot=0
    valid=0
    nb_flat=0
    compteur=0
    nb_neighbors=0
    # Résolution du système linéaire pour chaque triplet
    for i in range(len(voisins)):
        compteur+=1
        # On calcule b avec la multiplication de chaque élément du triplet
        b = np.array([mult(j) for j in voisins[i]])
        # print("")
        # print(f"Résolution pour le triplet {i}:")

        try:
            solution = np.linalg.solve(voisins[i], b)
            tot+=1
            if verification(solution,coeffs)==False:
                print(f"n+1-plet: {voisins[i]}")
                print(f"Valeur b: {b}")
                print("Solution:", solution)
            else:
                # print(f"n+1-plet: {voisins[i]}")
                # print(f"Valeur b: {b}")
                print("Solution:", solution)
                valid+=1
                b0=b[0]
                b_flat=True
                if is_neighbor(voisins[i]):
                    nb_neighbors+=1
                for t in range(len(b)):
                    if b[t]!=b0:
                        b_flat=False
                if b_flat:
                    print(b)
                    nb_flat+=1
        except np.linalg.LinAlgError as e:
            # print(f"Erreur lors de la résolution pour le triplet {triplets[i]} : {e}")
            continue
        # else:
            # print(f"Le triplet {i} n'est pas une matrice carrée. Impossible de résoudre.")
        # print(len(triplets)*tot/compteur)
    print("Nb sol = "+str(tot)+"/"+str(compteur))
    print("Nb valides = "+str(valid)+"/"+str(compteur))
    print("Nb neighbors = "+str(nb_neighbors)+"/"+str(compteur))
    print("Nb flat = "+str(nb_flat)+"/"+str(compteur))

def mult(x):
    tot = 1
    for k in range(len(x)-1):
        tot = tot * x[k]
    return tot

def verification(equa,coeffs):
    dim=len(equa)-1
    lmin=0
    lmax=0
    for k in coeffs:
        val=np.dot(equa,k)-mult(k)
        lmin=min(val,lmin)
        lmax=max(val,lmax)
    if lmin*lmax<0:
        return False
    # print(max(abs(lmin),abs(lmax)))
    return True

# Exemple d'appel
# calculateHyperplane(3)
# NeighborHyperplane(7)