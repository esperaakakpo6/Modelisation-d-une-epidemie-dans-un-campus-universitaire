import matplotlib.pyplot as plt
import ipywidgets as widgets
from ipywidgets import interact
import copy

# Initialisation des batiments
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
    20: {"S": 175, "I": 0,  "R": 0, "neighbors": [18]},
}

# Parametres du modele
beta = 0.3     # Taux de transmission
gamma = 0.1    # Taux de guerison

D_S = 0.01     # Diffusion des S
D_I = 0.01     # Diffusion des I
D_R = 0.01     # Diffusion des R

dt = 0.1       # Pas de temps
dx = 1.0       # Pas spatial (non utilise ici mais pour reference)
T = 50         # Duree totale
steps = int(T / dt)

# Implementation avec Euler explicite
# On suppose que "history" est une liste de snapshots de l'etat des batiments a chaque iteration
# Chaque element de history est un dict : {1: {"S": ..., "I": ..., "R": ...}, 2: {...}, ..., 20: {...}}

history = []

for t in range(steps):
    new_buildings = copy.deepcopy(buildings)
    
    for i in buildings:
        S = buildings[i]["S"]
        I = buildings[i]["I"]
        R = buildings[i]["R"]
        N = S + I + R

        # Terme de diffusion
        diffusion_S = sum((buildings[j]["S"] - S)/dx**2 for j in buildings[i]["neighbors"])
        diffusion_I = sum((buildings[j]["I"] - I)/dx**2  for j in buildings[i]["neighbors"])
        diffusion_R = sum((buildings[j]["R"] - R)/dx**2  for j in buildings[i]["neighbors"])

        # Mises a jour selon Euler explicite
        dS = -beta * S * I / N + D_S * diffusion_S
        dI = beta * S * I / N - gamma * I + D_I * diffusion_I
        dR = gamma * I + D_R * diffusion_R

        new_buildings[i]["S"] += dt * dS
        new_buildings[i]["I"] += dt * dI
        new_buildings[i]["R"] += dt * dR

    buildings = new_buildings
    history.append(copy.deepcopy(buildings))

# Visualisation
def plot_SIR(building_id):
    S_vals = [snapshot[building_id]["S"] for snapshot in history]
    I_vals = [snapshot[building_id]["I"] for snapshot in history]
    R_vals = [snapshot[building_id]["R"] for snapshot in history]

    plt.figure(figsize=(10, 5))
    plt.plot(S_vals, label="S (Susceptibles)", color='blue')
    plt.plot(I_vals, label="I (Infectés)", color='red')
    plt.plot(R_vals, label="R (Rétablis)", color='green')
    plt.title(f"Évolution SIR - Bâtiment {building_id}")
    plt.xlabel("Temps (itérations)")
    plt.ylabel("Nombre d’individus")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

#  Widget interactif : selection du batiment (entre 1 et 20)
interact(plot_SIR, building_id=widgets.IntSlider(min=1, max=20, step=1, value=1, description="Batiment ID"))