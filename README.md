# INF580 - DGP

Dans le cadre du projet INF580 - DGP, nous avons décidé de stocker les codes et les algorithmes pertinents dans un dossier nommé "utils". Ce dossier comprend quatre fichiers de code:

\begin{enumerate}
\item affichage.py - pour l'affichage de la trajectoire
\item annexe.py - pour la gestion des données (génération, rotation, projection, etc.)
\item optimization.py - pour résoudre le problème KEDM et pour la reconstruction (Application de Procrustes)
\item parametres.py - pour la gestion des paramètres tels que le nombre de points, le mode à utiliser, etc.
\end{enumerate}

Les codes AMPL pour la méthode de base sont stockés dans un dossier nommé "ampl". Ce dossier comprend cinq fichiers de code:
\begin{enumerate}
\item Deux fichiers .mod, l'un pour la méthode sans vitesse et l'autre avec vitesse.
\item Deux fichiers .run, l'un pour la méthode sans vitesse et l'autre avec vitesse.
\item Un fichier .dat servant de support pour stocker les données générées avec annexe.generate_trajectory().
\end {enumerate}

Le fichier principal INF580.ipynb est utilisé pour résoudre le problème d'optimisation avec les deux méthodes en sélectionnant les paramètres souhaités.