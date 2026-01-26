THERMONEUTRONICS — SIMULATION NEUTRONIQUE ET THERMIQUE

CONTEXTE PHYSIQUE
- Un réacteur nucléaire génère de la chaleur via des réactions de fission induites par des neutrons. La modélisation de ce processus implique la résolution couplée d’équations neutroniques et thermiques. Ce projet implémente un modèle simplifié en 2D cartésien pour un assemblage nucléaire typique, en résolvant les équations de diffusion neutronique et thermique couplées. On se place dans un assemblage nucléaire typique avec des barres de combustible, de l'eau qui modère le flux (et transporte la chaleur dans un réacteur) et des barres de contrôle qui absorbent totalement les neutrons. Dans ce milieu, les neutrons se diffusent, sont absorbés, et provoquent des fissions qui génèrent de la chaleur, chaque réaction de fission libére aussi des neutrons qui alimentent une réaction en chaîne, le but du réacteur étant de maintenir cette réaction de manière contrôlée en maximisant la production de chaleur. Un réacteur nucléaire réel est bien plus complexe, il comporte un grand nombre d'assemblages qui ont des géométries 3D et des configurations spécifiques de matériaux. De plus, de nombreux phénomènes physiques sont omis ici pour simplifier le modèle.

Les hypothèses principales utilsiées pour simplifier le modèle sont:
- Équation neutronique (1 groupe): diffusion avec pertes par absorption et production par fission. Le flux neutronique phi évolue sous l’effet de D (coefficient de diffusion), Sigma_a (section d’absorption), et nu*Sigma_f (source de fission).
- Équation thermique: équation de la chaleur avec conduction (lambda), capacité (rho*Cp) et source volumique Q proportionnelle à E_f * nu * Sigma_f * phi. Les grandeurs thermiques sont mises à l’échelle en unités SI via length_unit.
- Couplage: unidirectionnel dans la version actuelle (phi -> Q thermique -> T). La rétroaction de la température sur la neutronique n’est pas activée.
- Maillage et matériaux: assemblage 2D cartésien; chaque maille porte un matériau (UO2, H2O, B4C) avec voisins en grille (haut, bas, gauche, droite). 
- L'assemblage est construit avec des canaux carrés pouvant contenir des tubes de combustible ou des barres de contrôle (tubes "guide"). Dans le code, les paramètres N (nombre de canaux), pitch (pas entre canaux), r_fuel (rayon des barres de combustible), r_guide (rayon des barres de contrôle), guide_positions (positions des barres de contrôle) définissent la géométrie de l'assemblage. Le facteur d'insertion f_ins détermine si des barres de contrôle sont présentes (1) ou non (0). Dans un "pitch", un retrouve un crayon de combustible entouré d'eau, (ou potentiellement une barre de contrôle sisi la position correpsond à guide_positions). Le paramètre cells_per_pitch contrôle la résolution spatiale en définissant le nombre de mailles par pitch (pour un côté du pitch).


FONCTIONNEMENT GLOBAL
- Objectif: Simuler un assemblage nucléaire 2D couplé neutronique–thermique avec schéma en temps implicite et sorties VTK.
- Modules principaux:
  - src/material.* : Définit les matériaux standards (UO2, H2O, B4C) et écrit un VTK des matériaux du maillage.
  - src/mesh.* : Génère un maillage cartésien d’assemblage(s) avec positions de guides et voisinages.
  - src/solver.* : Implémente NeutronicSolver (diffusion + absorption + fission) et ThermalSolver (conduction + source volumique).
  - src/scheme.* : Gère les champs (phi, T), l’initialisation, l’export VTK, et le pas de temps (Euler implicite).
  - src/config.* : Lecture du fichier de configuration et affichage d’un résumé.
  - src/test_modes.* : Test avec des solutions analytiques (modes propres) pour valider le schéma neutronique.

ORDRE DU PROGRAMME PRINCIPAL
1. Lecture du fichier de configuration (ou déclenchement du mode test eigenmodes).
2. Construction du maillage et génération VTK des matériaux.
3. Initialisation du schéma numérique (flux neutronique phi, température T).
4. Boucle temporelle: avance du schéma (neutronique puis thermique) à chaque pas dt jusqu’à tf.
5. Sauvegarde périodique des résultats au format VTK dans results/.
6. Affichage de métriques finales (moyennes de phi et T).

COMPILATION
- Prérequis: Eigen3 installé et accessible par CMake. (sudo apt install libeigen3-dev avec Linux)
- Exécutable généré: build/simulation
- Note: Adapter le chemin vers Eigen si nécessaire dans CMakeLists.txt.
- Commandes :
  mkdir build
  cd build
  cmake ..
  make

EXÉCUTION
- Lancer avec le fichier par défaut:
  ./build/simulation
- Lancer avec un fichier spécifique:
  ./build/simulation config/config.txt
- Mode test eigenmodes (Neumann, UO2 uniforme):
  ./build/simulation test_eigenmode
  Affiche la métrique L2 relative entre solution numérique et solution exacte continue, pour plusieurs modes (n, m).

SORTIES ET VISUALISATION
- Résultats VTK dans results/:
  solution_N{N}_c{N*cells_per_pitch}_fins{f_ins}_t{step}.vtk
- Maillage matériaux VTK: mesh_materials.vtk (généré à l’initialisation).
- Visualisation:
  - Ouvrir les .vtk dans ParaView ou VisIt.
  - Choisir la donnée à afficher (phi ou T).
  
NOTES DE MODÈLE ET NUMÉRIQUE
- Neutronique: diffusion harmonisée, absorption Sigma_a, source fission nu*Sigma_f; opérateur discret assemblé dans NeutronicSolver::buildMatrix().
- Thermique: conduction (lambda), capacité volumique (rho*Cp), source volumique Q proportionnelle à E_f * nu * Sigma_f * phi, avec mises à l’échelle SI via length_unit.
- Conditions aux bords:
  - NEUMANN: flux normal nul (majoritairement utilisé ici)
- Schéma temporel: Euler implicite dans ImplicitEulerScheme pour les deux solveurs.

VALIDATION ET TEST EIGENMODES
- Le mode test_eigenmode initialise des modes cosinus cos(n*pi*x/L) * cos(m*pi*y/L) et calcule l’erreur L2 relative après intégration jusqu’à tf.
- Les erreurs augmentent avec les modes élevés à cause de la dispersion spatiale et de la discrétisation temporelle.

LIMITATIONS DU PROJET
- Modèle neutronique 1 groupe: pas de multi-groupe (neutrons rapides/thermiques séparés).
- Propriétés constantes: sections efficaces et propriétés thermiques constantes dans le temps et uniformes par maille (pas de dépendance en T, pas de combustion/burnup).
- Discrétisation spatiale simple: stencil de diffusion sur grille cartésienne; précision réduite aux bords (Neumann/Dirichlet simples).
- Schéma temporel: Euler implicite (ordre 1); pas de schémas d’ordre supérieur activés pour le neutronique.
- Couplage unidirectionnel: pas de rétroaction thermique sur le neutronique (pas de réactivité dépendante de T).
- Géométrie 2D: pas de 3D, pas de géométries complexes, pas de maillages non-cartésiens.
- Validation limitée par la difficulté d’obtenir des solutions analytiques exactes pour des configurations réalistes.
