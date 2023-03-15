# INF580 - DGP

Dans le cadre du projet INF580 - DGP, nous avons décidé de stocker les codes et les algorithmes pertinents dans un dossier nommé "utils". Ce dossier comprend quatre fichiers de code:

1. affichage.py - pour l'affichage de la trajectoire
2. annexe.py - pour la gestion des données (génération, rotation, projection, etc.)
3. optimization.py - pour résoudre le problème KEDM et pour la reconstruction (Application de Procrustes)
4. parametres.py - pour la gestion des paramètres tels que le nombre de points, le mode à utiliser, etc.

Les codes AMPL pour la méthode de base sont stockés dans un dossier nommé "ampl". Ce dossier comprend cinq fichiers de code:
1. Deux fichiers .mod, l'un pour la méthode sans vitesse et l'autre avec vitesse.
2. Deux fichiers .run, l'un pour la méthode sans vitesse et l'autre avec vitesse.
3. Un fichier .dat servant de support pour stocker les données générées avec annexe.generate_trajectory().

Le fichier principal INF580.ipynb est utilisé pour résoudre le problème d'optimisation avec les deux méthodes en sélectionnant les paramètres souhaités.
