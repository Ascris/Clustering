# **Test et développement d'un nouvel outil de classification à base de k-mers**
### Stage de fin de M1 Informatique au sein de l'équipe bioinformatique de l'IRHS d'Angers
##### Durée du stage (10 semaines) : _du 04 avril 2016 au 17 juin 2016_

---
Packages utilisées :
+  cluster (version 2.0.3)
+  seqinr (version 3.1-3)
+  ade4 (version 1.7-4)
+  datasets (version 3.2.4)
+  graphics (version 3.2.4)
+  grDevices (version 3.2.4)
+  grid (version 3.2.4)
+  methods (version 3.2.4)
+  stats (version 3.2.4)
+  utils (version 3.2.4)
+  igraph (version 1.0.1)
+  VennDiagram (version 1.6.17)

---

Pour toute première utilisation du programme, décommentez les premières lignes de scripts/main.r pour charger le fichier souhaité (fichier 403 ou fichier 1751 par exemple) Cette opération prend quelques minutes.
Une fois les données stockées dans la session R, re-commentez ces lignes pour pouvoir utiliser les fonctions sans avoir re-charger les données à chaque exécution.


#### **Protocole expérimental :**

Pour parvenir à catégoriser les X protéines en N groupes, nous avons dù passer par différentes étapes :
- Lecture de la matrice de distances (fichier .raw, matrice diagonale inférieure) et stockage des informations dans une matrice de taille X
- Application de PAM avec un nombre de clusters variant de min à max (pour une meilleure précision, on a appliqué chaque n-PAM k-fois)
- Création d'une matrice dite 'de robustesse' de taille X dans laquelle la case [i,j] correspond au nombre de fois où les protéines i et j ont été classées dans le même cluster par PAM
- Visualisation graphique des groupes sous forme de cliques (dont chaque arc lie 2 protéines par une valeur 'robustesse' variable)
- Formation des groupes faite par 'recherche de réseau d'amis' (ensemble des protéines s'étant retrouvées au moins 'robustesse' fois dans le même groupe)
- Stockage des réseaux d'amis dans un fichier txt avec ajout du groupe des 'singletons', protéines non liées à quelconque groupe

#### Étapes importantes :
- [x] lecture fichier .raw
- [x] création matrice de robustesse
- [x] visualisation graphique sous forme de clique
- [x] recherche des réseaux d'amis
- [x] écriture/stockage des résultats dans les fichiers

Détails/difficultés des étapes :
=================================
1. **Lecture du fichier raw :**
	Le fichier d'entrée étant une matrice diagonale inférieure (problème de lecture au début, la dernière protéine (taille n) a été placée en première position dans le .raw)
	- extraction des noms des protéines en première colonne
	- chaque ligne contient de 1 à n 'mots' mais chaque protéine partage une valeur avec les n-1 autres

2. **Création de la matrice de robustesse :**
    
	PAM (Partitional Around Medoids) est un algorithme permettant, à partir d'une matrice de distances de taille X et d'un paramètre k, de répartir les X élements de la matrice en k-groupes distincts. L'idée a alors été d'appliquer PAM à notre jeu de données en faisant varier k et ainsi d'observer si certaines protéines se retrouvaient systématiquement dans le même groupe. Si tel est le cas, on peut en déduire que les deux protéines en question ont une forme forte.

	- On différenciera la 'matrice de robustesse globale' des 'matrices de robustesse locales'
		- on a appliqué X fois l'algorithme PAM sur chaque k fixé pour éviter les erreurs dues à l'aléatoire. On a obtenu une 'matrice de robustesse locale' pour chaque k testé.
		- la matrice de robustesse globale 'matG' somme les matrices de robustesse locales 'matL'. On obtient donc pour chaque case [i,j] de 'matG' le nombre de fois que les protéines i et j ont été classées ensemble par les k-applications de PAM.

3. **Visualisation graphique**

	Pour pouvoir simplement observer les différents groupes révélés par les itérations de PAM, deux fichiers ont été créés : l'un contenant les sommets (les protéines présentes dans la matrice de robustesse) et l'autre contenant les arcs (lien de valeur 'robustesse' (variable) entre deux protéines A et B). On crée ensuite un graphe à partir de ces deux fichiers auquel on ajoute deux paramètres rendant les protéines affichées plus parlantes : la couleur pour son milieu (T, M, P ou H) ou son type (TOP ou RG) (selon le graphique) et la taille pour son règne (A, B ou E).

	Dans le but de tester la stabilité des groupes formés (que ce soit par PAM, hclust et cutree ou par la recherche de réseaux d'amis), plusieurs méthodes ont été testées.On s'intéresse ici uniquement à la méthode avec le cutree et on appellera 'coupe' l'ensemble des groupes trouvés qui découlent de cette méthode. Les coupes auxquelles on s'intéresse sont créées par la méthode "ward.D2" demandée en paramètre par hclust.

---
	Le paramètre que l'on va faire varier est la matrice de données. On va alors comparer les groupes obtenus directement grâce à la matrice de distance (obtenue après calculs par la méthode VLD), ceux obtenus avec la matrice de robustesse (obtenue après k-itérations de PAM) et ceux obtenus avec une sous-matrice de robustesse à laquelle on lui a retiré un certain pourcentage d'individus (10% ici).

Pour tester les groupes, on a alors :
- une matrice de distances
-  une matrice de robustesse
- une sous-matrice de robustesse 

---
Ces 3 types de coupes sont regroupées au sein d'un même diagramme de Venn pour pouvoir être comparées.

4. **Recherche des réseaux d'amis**

	*"Les amis de mes amis sont mes amis"*

	+ On utilisera le mot 'ami' pour désigner qu'il existe un lien de robustesse fort entre deux protéines A et B (i.e. matG[A,B] >= 'robustesse'). Deux protéines amies se retrouvent systématiquement dans le même groupe lorsque l'on applique PAM avec n'importe quel paramètre k.
	+ On appelera réseau d'amies de X, avec X une protéine quelconque, l'ensemble des protéines qui sont amies avec X ou une amie de X.
	
	- fonction récursive nécessaire ici, coût en mémoire relativement important - à surveiller
	- on utilise une fonction qui récupère la liste des amies de X et qui l'ajoute à une liste totale d'amies de X
	- pour chaque amie Y de cette liste, on ajoute à la liste totale d'amies de X la liste d'amies de Y
	- on réitère jusqu'à ne plus avoir de nouvelle amie à ajouter aux amies de X
	- on applique cet algorithme à toutes les protéines, on obtient alors une liste d'amies pour chacune d'elles
	- on retire les doublons
	
5. **Écriture/stockage des résultats dans les fichiers**

	L'organisation de mon espace de travail s'est faite de cette manière : un dossier par type de fichier traité. Les répertoires 'scripts' et 'fonctions' contiennent les algorithmes utilisés ainsi que les fonctions appelées par ces derniers. Le répertoire 'data' contient les données : fichiers .raw contenant les matrices de distance ainsi que les fichiers sommets et arcs créés pour tracer les graphiques. Enfin, le répertoire 'resultats' contient tous les fichiers graphiques (heatmaps, cliques, diagrammes de Venn), les fichiers liés au temps de calcul et les matrices de robustesse locales qui ont été préservées.
	
	- gérer les noms de fichiers en fonction de ceux déjà présents dans le répertoire
	- parser les identifiants des protéines pour en extraire les caractéristiques (nom, règne, type et milieu)
