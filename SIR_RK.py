import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from ipywidgets import interact, IntSlider
from IPython.display import display
import copy

# Initialisation des parametres
buildings = {
    1:  {"S": 100, "I": 1,  "R": 0, "neighbors": [2, 3]},
    2:  {"S": 120, "I": 0,  "R": 0, "neighbors": [1, 4]},
    3:  {"S": 80,  "I": 0,  "R": 0, "neighbors": [1, 5]},
    4:  {"S": 150, "I": 0,  "R": 0, "neighbors": [2, 6]},
    5:  {"S": 90,  "I": 0,  "R": 0, "neighbors": [3, 7]},
    6:  {"S": 110, "I": 0,  "R": 0, "neighbors": [4, 8]},
    7:  {"S": 95,  "I": 0,  "R": 0, "neighbors": [5, 9]},
    8:  {"S": 130, "I": 0,  "R": 0, "neighbors": [6, 10]},
    9:  {"S": 85,  "I": 0,  "R": 0, "neighbors": [7, 11]},
    10: {"S": 100, "I": 0,  "R": 0, "neighbors": [8, 12]},
    11: {"S": 75,  "I": 0,  "R": 0, "neighbors": [9, 13]},
    12: {"S": 140, "I": 0,  "R": 0, "neighbors": [10, 14]},
    13: {"S": 105, "I": 0,  "R": 0, "neighbors": [11, 15]},
    14: {"S": 115, "I": 0,  "R": 0, "neighbors": [12, 16]},
    15: {"S": 125, "I": 0,  "R": 0, "neighbors": [13, 17]},
    16: {"S": 135, "I": 0,  "R": 0, "neighbors": [14, 18]},
    17: {"S": 145, "I": 0,  "R": 0, "neighbors": [15, 19]},
    18: {"S": 155, "I": 0,  "R": 0, "neighbors": [16, 20]},
    19: {"S": 165, "I": 0,  "R": 0, "neighbors": [17]},
    20: {"S": 175, "I": 0,  "R": 0, "neighbors": [18]}
}

beta = 0.3      # taux d'infection
gamma = 0.1     # taux de guérison
DS = 0.01       # diffusion S
DI = 0.01       # diffusion I
DR = 0.01       # diffusion R
dx = 1.0        # distance spatiale
temps_o = 0
temps_total = 100.0  # temps total de simulation
iteration = 500 #Nonbre d'iteration
dt = (temps_total - temps_o)/iteration        # pas de temps


#Discrétisation spatiale (EDP)
def laplacian(var, i, buildings):
    neighbors = buildings[i]["neighbors"]
    return sum((var[j] - var[i]) / dx**2 for j in neighbors)


# Fonction de la methode de Runge Kutta a une seule ittération
def rk4_step(buildings, dt):
    def get_states():
        S = {i: buildings[i]["S"] for i in buildings}
        I = {i: buildings[i]["I"] for i in buildings}
        R = {i: buildings[i]["R"] for i in buildings}
        return S, I, R

    S0, I0, R0 = get_states()

    def deriv(S, I, R):
        dS, dI, dR = {}, {}, {}
        for i in buildings:
            N = S[i] + I[i] + R[i]
            dS[i] = -beta * S[i] * I[i] / N + DS * laplacian(S, i, buildings)
            dI[i] = beta * S[i] * I[i] / N - gamma * I[i] + DI * laplacian(I, i, buildings)
            dR[i] = gamma * I[i] + DR * laplacian(R, i, buildings)
        return dS, dI, dR

    # k1
    k1S, k1I, k1R = deriv(S0, I0, R0)

    # k2
    S2 = {i: S0[i] + k1S[i]/2 for i in buildings}
    I2 = {i: I0[i] + k1I[i]/2 for i in buildings}
    R2 = {i: R0[i] + k1R[i]/2 for i in buildings}
    k2S, k2I, k2R = deriv(S2, I2, R2)

    # k3
    S3 = {i: S0[i] + k2S[i]/2 for i in buildings}
    I3 = {i: I0[i] + k2I[i]/2 for i in buildings}
    R3 = {i: R0[i] + k2R[i]/2 for i in buildings}
    k3S, k3I, k3R = deriv(S3, I3, R3)

    # k4
    S4 = {i: S0[i] + k3S[i] for i in buildings}
    I4 = {i: I0[i] + k3I[i] for i in buildings}
    R4 = {i: R0[i] + k3R[i] for i in buildings}
    k4S, k4I, k4R = deriv(S4, I4, R4)

    # Mise à jour
    for i in buildings:
        buildings[i]["S"] += (dt/6) * (k1S[i] + 2*k2S[i] + 2*k3S[i] + k4S[i])
        buildings[i]["I"] += (dt/6) * (k1I[i] + 2*k2I[i] + 2*k3I[i] + k4I[i])
        buildings[i]["R"] += (dt/6) * (k1R[i] + 2*k2R[i] + 2*k3R[i] + k4R[i])
    
    return buildings
# Calcul sur plusieurs itterations
temps = []
dict_list = []
for t in range(iteration):
    temps.append(temps_o + t*dt)
    dict_list.append(rk4_step(buildings, dt))
    buildings = copy.deepcopy(dict_list[t])
    
import numpy as np

def construire_matrices(dict_list):
    # Extraire et trier les identifiants des bâtiments
    cles = sorted(dict_list[0].keys())
    
    S_matrix = []
    I_matrix = []
    R_matrix = []

    for i in cles:
        ligne_S = [d[i]["S"] for d in dict_list]
        ligne_I = [d[i]["I"] for d in dict_list]
        ligne_R = [d[i]["R"] for d in dict_list]
        
        S_matrix.append(ligne_S)
        I_matrix.append(ligne_I)
        R_matrix.append(ligne_R)
        
    return np.array(S_matrix), np.array(I_matrix), np.array(R_matrix)

S_matrix, I_matrix, R_matrix = construire_matrices(dict_list)


# Representation de la solution avec ipwidgets
# Thème esthétique
# plt.style.use('seaborn-darkgrid')
matplotlib.rcParams.update({'font.size': 12})


x_values = np.arange(S_matrix.shape[1])

# --- Fonction d'affichage ---
def plot_sir(batiment_index):
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(temps, S_matrix[batiment_index], label="Sains (S)", color="#1f77b4", linewidth=2.5)
    ax.plot(temps, I_matrix[batiment_index], label="Infectés (I)", color="#ff7f0e", linewidth=2.5)
    ax.plot(temps, R_matrix[batiment_index], label="Rétablis (R)", color="#2ca02c", linewidth=2.5)

    ax.set_title(f"Évolution SIR – Bâtiment {batiment_index + 1}", fontsize=15, weight='bold')
    ax.set_xlabel("Temps", fontsize=13)
    ax.set_ylabel("Population", fontsize=13)
    ax.legend(loc="best", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.set_ylim(0, max(S_matrix.max(), I_matrix.max(), R_matrix.max()) + 10)

    plt.tight_layout()
    plt.show()

# --- Interface avec slider ---
interact(plot_sir, batiment_index=IntSlider(min=0, max=S_matrix.shape[0]-1, step=1, value=0, description="Bâtiment"))
