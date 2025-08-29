import matplotlib.pyplot as plt
import ipywidgets as widgets
from ipywidgets import interact
import copy

# Paramètres
beta = 0.3
gamma = 0.1
D_S = 0.01
D_I = 0.01
D_R = 0.01
dt = 0.1
dx = 1.0
T = 50
steps = int(T / dt)
dx2 = 1.0

# Initialisation
initial_buildings = {
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

def simulate_euler_implicite(buildings_init, steps):
    # Copie de l'état initial des bâtiments pour ne pas modifier l'entrée originale
    buildings = copy.deepcopy(buildings_init)
    
    # Historique de l'évolution SIR à chaque pas de temps
    history = []

    # Boucle principale de temps (de 0 à T avec un pas de dt)
    for _ in range(steps):
        # Initialisation des nouvelles valeurs de S, I, R pour tous les bâtiments
        S_new = {i: buildings[i]["S"] for i in buildings}
        I_new = {i: buildings[i]["I"] for i in buildings}
        R_new = {i: buildings[i]["R"] for i in buildings}

        # Méthode d'Euler implicite : on applique une résolution itérative par point fixe
        for _ in range(20):  # Nombre d’itérations de point fixe (convergence approximative)
            # Sauvegarde des valeurs précédentes (à l’itération k-1)
            S_prev = S_new.copy()
            I_prev = I_new.copy()
            R_prev = R_new.copy()

            # Mise à jour de chaque bâtiment individuellement
            for i in buildings:
                neighbors = buildings[i]["neighbors"]  # Voisins spatiaux
                Ni = S_prev[i] + I_prev[i] + R_prev[i]  # Population totale du bâtiment i

                # Terme de diffusion : somme des différences avec les voisins
                diff_S = sum((S_prev[j] - S_prev[i])/dx**2  for j in neighbors)
                diff_I = sum((I_prev[j] - I_prev[i])/dx**2  for j in neighbors)
                diff_R = sum((R_prev[j] - R_prev[i])/dx**2  for j in neighbors)

                # Mise à jour des variables S, I, R selon le schéma implicite
                # S^{n+1}_i = S^n_i + dt × [ -β S^{n+1} I^{n+1} / N + D × (diffusion) ]
                S_new[i] = buildings[i]["S"] + dt * (
                    -beta * S_prev[i] * I_prev[i] / Ni + D_S * diff_S / dx2
                )
                I_new[i] = buildings[i]["I"] + dt * (
                    beta * S_prev[i] * I_prev[i] / Ni - gamma * I_prev[i] + D_I * diff_I / dx2
                )
                R_new[i] = buildings[i]["R"] + dt * (
                    gamma * I_prev[i] + D_R * diff_R / dx2
                )

        # Mise à jour finale des valeurs des bâtiments avec les valeurs implicites convergées
        for i in buildings:
            buildings[i]["S"] = S_new[i]
            buildings[i]["I"] = I_new[i]
            buildings[i]["R"] = R_new[i]

        # Enregistrement de l’état du système à cet instant dans l’historique
        history.append(copy.deepcopy(buildings))

    # Retourne l'évolution complète sur tout l'intervalle de temps
    return history


# Simulation
history_implicit = simulate_euler_implicite(initial_buildings, steps)

# Affichage interactif
def plot_SIR_implicit(building_id):
    S_vals = [snapshot[building_id]["S"] for snapshot in history_implicit]
    I_vals = [snapshot[building_id]["I"] for snapshot in history_implicit]
    R_vals = [snapshot[building_id]["R"] for snapshot in history_implicit]

    plt.figure(figsize=(10, 5))
    plt.plot(S_vals, label="S (Susceptibles)", color='blue')
    plt.plot(I_vals, label="I (Infectés)", color='red')
    plt.plot(R_vals, label="R (Rétablis)", color='green')
    plt.title(f"Évolution SIR (Euler implicite) - Bâtiment {building_id}")
    plt.xlabel("Temps (itérations)")
    plt.ylabel("Nombre d’individus")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

interact(plot_SIR_implicit, building_id=widgets.IntSlider(min=1, max=20, step=1, value=1, description="Bâtiment"));
