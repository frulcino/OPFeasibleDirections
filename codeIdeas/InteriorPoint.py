import numpy as np
from scipy.optimize import minimize

# Fonction objectif
def objective(x):
    return -(x[0] + x[1])  # On maximise x1 + x2, donc on minimise -(x1 + x2)

# Contrainte de norme
def constraint(x):
    return 4 - (x[0]**3 + x[1]**2)  # x1^2 + x2^2 <= 9

# Point de départ
x0 = np.array([1.0, 1.0])

# Contrainte sous forme de dictionnaire
con = {'type': 'ineq', 'fun': constraint}

# Appel à l'optimiseur pour maximiser x1 + x2
result = minimize(objective, x0, constraints=con)

# Affichage des résultats
print("Solution optimale : ", result.x)
print("Valeur de la fonction objectif : ", -result.fun)  # On maximise, donc on prend l'opposé
