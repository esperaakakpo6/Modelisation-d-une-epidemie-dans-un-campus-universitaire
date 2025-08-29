import matplotlib.pyplot as plt
import ipywidgets as widgets
import matplotlib
from ipywidgets import interact, IntSlider, FloatSlider
import numpy as np
import copy

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


        # Debut euler implicite pour plusieur combinason de dt et dx
# Paramètres
beta = 0.3
gamma = 0.1
D_S = 0.05
D_I = 0.05
D_R = 0.05
T = 50
To = 0
list_dx = [1/2, 2/2, 3/2, 4/2, 5/2]
list_dt = [0.1*i for i in range(1, 6)]
# D = max(D_I, D_R, D_S)
# list_dx = np.linspace(0.01, 0.1, 101)[:5]
# list_dt = []
# for i in range(5):
#     dx_i = list_dx[i]
#     list_dt.append(round(min(2/(gamma+beta), (dx_i**2)/(2*D)), 5))

for dt in list_dx:
    steps = int(T / dt)
    for dx in list_dx:
        dx2 = dx**2   
        def simulate_euler_implicite(buildings_init, steps):
            buildings = copy.deepcopy(buildings_init)
            history = []
            for _ in range(steps):
                S_new = {i: buildings[i]["S"] for i in buildings}
                I_new = {i: buildings[i]["I"] for i in buildings}
                R_new = {i: buildings[i]["R"] for i in buildings}

                for _ in range(20):  # Itération de point fixe
                    S_prev = S_new.copy()
                    I_prev = I_new.copy()
                    R_prev = R_new.copy()

                    for i in buildings:
                        neighbors = buildings[i]["neighbors"]
                        Ni = S_prev[i] + I_prev[i] + R_prev[i]

                        diff_S = sum(S_prev[j] - S_prev[i] for j in neighbors)
                        diff_I = sum(I_prev[j] - I_prev[i] for j in neighbors)
                        diff_R = sum(R_prev[j] - R_prev[i] for j in neighbors)

                        S_new[i] = buildings[i]["S"] + dt * (-beta * S_prev[i] * I_prev[i] / Ni + D_S * diff_S / dx2)
                        I_new[i] = buildings[i]["I"] + dt * (beta * S_prev[i] * I_prev[i] / Ni - gamma * I_prev[i] + D_I * diff_I / dx2)
                        R_new[i] = buildings[i]["R"] + dt * (gamma * I_prev[i] + D_R * diff_R / dx2)

                for i in buildings:
                    buildings[i]["S"] = S_new[i]
                    buildings[i]["I"] = I_new[i]
                    buildings[i]["R"] = R_new[i]

                history.append(copy.deepcopy(buildings))
            return history

        # Simulation
        history_implicit = simulate_euler_implicite(initial_buildings, steps)

        S_matrix = np.zeros((len(initial_buildings), steps))
        I_matrix = np.zeros((len(initial_buildings), steps))
        R_matrix = np.zeros((len(initial_buildings), steps))
        for i in range(1, 21):
            S_matrix[i-1, :] = [snapshot[i]["S"] for snapshot in history_implicit]
            I_matrix[i-1, :] = [snapshot[i]["I"] for snapshot in history_implicit]
            R_matrix[i-1, :] = [snapshot[i]["R"] for snapshot in history_implicit]
            
# Global storage for all results
all_S_matrices = {}
all_I_matrices = {}
all_R_matrices = {}

for dt in list_dt:
    steps = int(T / dt)
    for dx in list_dx:
        dx2 = dx**2

        def simulate_euler_implicite(buildings_init, steps):
            buildings = copy.deepcopy(buildings_init)
            history = []
            for _ in range(steps):
                S_new = {i: buildings[i]["S"] for i in buildings}
                I_new = {i: buildings[i]["I"] for i in buildings}
                R_new = {i: buildings[i]["R"] for i in buildings}

                for _ in range(20):  # Fixed-point iteration
                    S_prev = S_new.copy()
                    I_prev = I_new.copy()
                    R_prev = R_new.copy()

                    for i in buildings:
                        neighbors = buildings[i]["neighbors"]
                        Ni = S_prev[i] + I_prev[i] + R_prev[i]

                        diff_S = sum(S_prev[j] - S_prev[i] for j in neighbors)
                        diff_I = sum(I_prev[j] - I_prev[i] for j in neighbors)
                        diff_R = sum(R_prev[j] - R_prev[i] for j in neighbors)

                        S_new[i] = buildings[i]["S"] + dt * (-beta * S_prev[i] * I_prev[i] / Ni + D_S * diff_S / dx2)
                        I_new[i] = buildings[i]["I"] + dt * (beta * S_prev[i] * I_prev[i] / Ni - gamma * I_prev[i] + D_I * diff_I / dx2)
                        R_new[i] = buildings[i]["R"] + dt * (gamma * I_prev[i] + D_R * diff_R / dx2)

                for i in buildings:
                    buildings[i]["S"] = S_new[i]
                    buildings[i]["I"] = I_new[i]
                    buildings[i]["R"] = R_new[i]

                history.append(copy.deepcopy(buildings))
            return history

        # Simulation
        history_implicit = simulate_euler_implicite(copy.deepcopy(initial_buildings), steps)

        S_matrix = np.zeros((len(initial_buildings), steps))
        I_matrix = np.zeros((len(initial_buildings), steps))
        R_matrix = np.zeros((len(initial_buildings), steps))

        for i in range(1, len(initial_buildings) + 1): # Iterate through building IDs
            S_matrix[i-1, :] = [snapshot[i]["S"] for snapshot in history_implicit]
            I_matrix[i-1, :] = [snapshot[i]["I"] for snapshot in history_implicit]
            R_matrix[i-1, :] = [snapshot[i]["R"] for snapshot in history_implicit]

        # Store the matrices using a unique key (e.g., a tuple of dt and dx)
        key = (dt, dx)
        all_S_matrices[key] = S_matrix
        all_I_matrices[key] = I_matrix
        all_R_matrices[key] = R_matrix

# Example of how to access a specific matrix:
# If you want the S_matrix for dt=0.1 and dx=1.0:
s_matrix_for_0_1_1_0 = all_S_matrices[(0.1, 1.0)]
print(f"Shape of S_matrix for dt=0.1, dx=1.0: {s_matrix_for_0_1_1_0.shape}")

    # Debut runge kutta
DS = D_S
DI = D_I
DR = D_R

dx = 1.0        # distance spatiale
dt = 0.01        # pas de temps
temps_o = To
temps_total = T  # temps total de simulation
iteration = (temps_total - temps_o)/dt        # pas de temps

def laplacian(var, i, buildings):
    neighbors = buildings[i]["neighbors"]
    return sum((var[j] - var[i]) / dx**2 for j in neighbors)

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

temps = []
dict_list = []
buildings = copy.deepcopy(initial_buildings)
for t in range(int(iteration)):
    temps.append(temps_o + t*dt)
    dict_list.append(rk4_step(buildings, dt))
    buildings = copy.deepcopy(dict_list[t])
    


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

def return_diff_squared(i, dt, dt_var, dx_var, matrix_rk, big_matrix):
    error_i = 0
    for j in range(big_matrix[(dt_var, dx_var)].shape[1]): 
        pos = int((j*dt_var)/dt)
        error_i += (big_matrix[(dt_var, dx_var)][i, j] - matrix_rk[i, pos])**2
    error_i = (dt_var*error_i)**(1/2)
    return error_i

all_error_matrix = np.zeros(len(list_dx))
all_S_error = {}
all_I_error = {}
all_R_error = {}

for matrix_rk, big_matrix, lettre in zip([S_matrix, I_matrix, R_matrix], [all_S_matrices, all_I_matrices, all_R_matrices], ["S", "I", "R"]):
    for j in range(len(list_dx)):
        dx_var = list_dx[j]
        error_matrix = np.zeros((len(initial_buildings), len(list_dt)))
        for i in range(len(initial_buildings)):
            for k in range(len(list_dt)):
               dt_var = list_dt[k]
               error_matrix[i,k] = return_diff_squared(i, dt, dt_var, dx_var, matrix_rk, big_matrix)
                
        key = (dx_var)
        if lettre == "S":
             all_S_error[key] = error_matrix
        if lettre == "I":
             all_I_error[key] = error_matrix
        if lettre == "R":
             all_R_error[key] = error_matrix
        
        
# X-axis values
x_values = list_dt

# Create sliders
key_slider = widgets.FloatSlider(
    value=0.5,
    min=0.5,
    max=2.5,
    step=0.5,
    description='Delta x:',
    continuous_update=False
)

row_slider = widgets.IntSlider(
    value=0,
    min=0,
    max=19,
    step=1,
    description='Batiment:',
    continuous_update=False
)

# Output widget for the plot
output = widgets.Output()

# Function to update the plot
def update_plot(change):
    with output:
        clear_output(wait=True)
        key = key_slider.value
        row = row_slider.value
        S_values = all_S_error[key][row]
        I_values = all_I_error[key][row]
        R_values = all_R_error[key][row]
        
        plt.figure(figsize=(15, 5))
        plt.plot(x_values, S_values, label="Sains (S)", color="#0d00ff", linewidth=1)
        plt.plot(x_values, I_values, label="Infectés (I)", color="#ff0000", linewidth=1)
        plt.plot(x_values, R_values, label="Rétablis (R)", color="#04ff04", linewidth=1)
        plt.title(f'Batiment {row+1} - Delta {key+1}')
        plt.xlabel('Delta t')
        plt.ylabel('Erreur')
        plt.legend()
        plt.grid(True)
        plt.show()

# Connect sliders to the update function
key_slider.observe(update_plot, names='value')
row_slider.observe(update_plot, names='value')

# Display sliders and initial plot
display(key_slider, row_slider, output)

# Initial plot
update_plot(None)