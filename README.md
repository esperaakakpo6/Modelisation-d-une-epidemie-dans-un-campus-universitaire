# Modélisation d'une Épidémie dans un Campus Universitaire

## Description
Ce projet modélise la propagation d'une épidémie de grippe dans un campus composé de 20 bâtiments interconnectés, en utilisant une version spatialisée du modèle SIR (Susceptibles - Infectés - Rétablis). Il intègre des termes de diffusion spatiale pour simuler les déplacements entre bâtiments.

Les méthodes numériques implémentées incluent :
- Euler explicite
- Euler implicite
- Runge-Kutta d'ordre 4

Le projet analyse la stabilité numérique (condition CFL), la convergence (erreur L2), et propose des visualisations interactives avec ipywidgets. Un rapport détaillé au format PDF est fourni.

## Structure des Fichiers
- **erreur.py** : Calcul des erreurs de convergence (L2) entre méthodes, avec visualisations interactives des erreurs pour différents pas de temps/spatiaux.
- **explicite_expose.py** : Implémentation de la méthode Euler explicite, avec visualisation interactive SIR par bâtiment.
- **implicite_expose.py** : Implémentation de la méthode Euler implicite, avec visualisation interactive SIR par bâtiment.
- **SIR_RK.py** : Implémentation de la méthode Runge-Kutta d'ordre 4, avec visualisation interactive SIR par bâtiment.
- **Modélisation_d_une_épidémie_dans_un_campus_universitaire.pdf** : Rapport complet (38 pages) couvrant la modélisation, les méthodes numériques, l'analyse et des applications (ex. : mesures sanitaires localisées).

## Prérequis
- Python 3.x
- Bibliothèques : 
