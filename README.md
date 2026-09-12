# Introduction au machine learning pour la génomique de l'antibiorésistance

Support du module **DU Bioinformatique, UM6SS (2026)** consacré à
l'apprentissage supervisé et non supervisé appliqué à la résistance aux
antimicrobiens chez *Escherichia coli*.

Le projet consiste à prédire la résistance à la ciprofloxacine à partir de
caractéristiques d'assemblage, de métadonnées de collecte et des résultats AST
pour sept autres antibiotiques. Il combine :

- un phénotypage non supervisé avec `phynotype` ;
- une comparaison de modèles supervisés avec `funcml` ;
- une imputation adaptée aux données mixtes avec `missknn` ;
- une lecture biologique et méthodologique des résultats.

> **Document de travail obligatoire**
>
> Toutes les questions du projet se trouvent dans
> [`project/projet-enonce.pdf`](project/projet-enonce.pdf). Ce PDF est
> l'énoncé de référence : l'étudiant doit le suivre du début à la fin et
> répondre aux questions dans l'ordre. Le présent README explique comment
> s'organiser et exécuter le travail ; il ne remplace pas l'énoncé.

## Contenu du dépôt

```text
.
├── README.md
├── project/
│   ├── project.Rproj
│   ├── projet-enonce.qmd
│   ├── projet-enonce.html
│   ├── projet-enonce.pdf
│   ├── ecoli_amr_isolate_v2.csv
│   ├── build_clean_data.R
│   └── raw/
│       ├── ecoli_amr_raw.csv
│       └── ecoli_metadata_raw.csv
├── slides/
│   ├── cours-ml-genomique-amr.qmd
│   └── cours-ml-genomique-amr.pdf
└── refs/
    ├── references.bib
    └── articles et documentation du cours
```

Le corrigé reste réservé à l'enseignant pendant la réalisation du projet et
n'est pas accessible dans ce dépôt. Il sera partagé avec les étudiants après
la remise des travaux, afin qu'ils puissent comparer leur démarche, comprendre
leurs erreurs et reprendre l'analyse pas à pas. Jusqu'à cette date, seul
l'énoncé fait foi.

## Prérequis

Avant de commencer, installer :

- une version récente de [R](https://cran.r-project.org/) ;
- au choix, [RStudio Desktop](https://posit.co/download/rstudio-desktop/),
  [Visual Studio Code](https://code.visualstudio.com/) ou
  [Positron](https://positron.posit.co/) ;
- [Quarto](https://quarto.org/docs/get-started/) ;
- facultativement, une distribution LaTeX pour produire le PDF. L'export HTML
  ne nécessite pas LaTeX.

Avec Visual Studio Code, installer les extensions **R** et **Quarto** depuis le
catalogue des extensions. RStudio intègre directement la console R et reconnaît
le fichier `project.Rproj`. Positron intègre nativement R et Quarto.

Une connexion Internet est nécessaire lors de la première installation des
packages.

## Installation des packages R

Ouvrir R, puis exécuter une seule fois :

```r
install.packages("pak")

pak::pak(c(
  "ielbadisy/funcml",
  "ielbadisy/missknn",
  "ielbadisy/phynotype"
))
```

`fastgbm` et `densemlp` sont des dépendances de `funcml`. Il n'est pas
nécessaire de les installer ou de les charger séparément.

Pour vérifier l'installation :

```r
library(funcml)
library(missknn)
library(phynotype)
```

## Démarrage recommandé

1. Cloner le dépôt :

   ```bash
   git clone https://github.com/ielbadisy/intro-ml4genomics-amr-2026.git
   cd intro-ml4genomics-amr-2026
   ```

2. Choisir un environnement de travail :

   - avec **RStudio**, ouvrir `project/project.Rproj` ;
   - avec **Visual Studio Code**, ouvrir le dossier racine du dépôt, puis
     utiliser un terminal R avec `project/` comme répertoire courant ;
   - avec **Positron**, ouvrir directement le dossier `project/`, puis démarrer
     la console R intégrée.

3. Ouvrir [`project/projet-enonce.pdf`](project/projet-enonce.pdf) et lire
   entièrement l'énoncé avant de commencer. La
   [version HTML](project/projet-enonce.html) et le
   [source Quarto](project/projet-enonce.qmd) contiennent les mêmes consignes.
4. Lire `project/build_clean_data.R` afin de comprendre la provenance et
   l'échantillonnage des données.
5. Créer dans `project/` un nouveau document, par exemple
   `projet-nom-prenom.qmd`.
6. Travailler avec `project/` comme répertoire courant. Le fichier de données
   se charge alors avec :

   ```r
   raw <- basetable::btread("ecoli_amr_isolate_v2.csv")
   ```

Les trois environnements conviennent. Les commandes R, le document Quarto et
les résultats attendus sont identiques dans RStudio, Visual Studio Code et
Positron.

Ne pas modifier les fichiers bruts. Toute transformation doit être exprimée
dans le document reproductible.

## Méthode de travail pas à pas

Le projet n'est pas une liste d'analyses indépendantes. Chaque étape produit
les objets et les arguments nécessaires à l'étape suivante. Pour chaque
question du PDF :

1. recopier son numéro et son intitulé comme section du rapport personnel ;
2. identifier les données d'entrée et l'objet R attendu en sortie ;
3. écrire et exécuter le code correspondant ;
4. afficher le tableau, la métrique ou le graphique demandé ;
5. ajouter immédiatement une interprétation en quelques phrases ;
6. rendre régulièrement le document pour détecter les erreurs de
   reproductibilité avant de poursuivre.

Ne pas passer à la partie supervisée tant que le prétraitement et le
phénotypage non supervisé ne fonctionnent pas. Le rapport final doit conserver
la même progression que le PDF : **prétraitement → non supervisé → supervisé
→ lecture croisée → conclusion**.

## Travail demandé

Le document final doit répondre, dans l'ordre, à toutes les questions
numérotées de [`project/projet-enonce.pdf`](project/projet-enonce.pdf). Le
résumé ci-dessous sert de feuille de route ; en cas de doute, le texte complet
du PDF fait foi.

### Partie 1 : Prétraitement

1. Construire la cible binaire `cipro_R` à partir de `ciprofloxacin`, en
   conservant uniquement `Resistant` et `Susceptible`.
2. Distinguer les deux mécanismes de données manquantes :
   - AST absent : coder explicitement `Not_tested`, sans imputation ;
   - métadonnée absente : imputer avec `missknn` et justifier le choix.
3. Regrouper les modalités rares des variables de pays et d'hôte.
4. Calculer la prévalence de la résistance et expliquer ses conséquences sur
   les métriques.

### Partie 2 : Apprentissage non supervisé

1. Préparer les variables mixtes sans inclure `cipro_R`.
2. Comparer et valider plusieurs nombres de groupes ; justifier le choix de
   `k` avec les indices internes et la stabilité bootstrap.
3. Examiner les profils, tailles, projections et au moins une méthode
   d'interprétabilité.
4. Croiser ensuite seulement les phénotypes obtenus avec `cipro_R`.

La cible ne doit jamais intervenir dans la construction des groupes. Son usage
est autorisé uniquement pour l'évaluation a posteriori de leur association à
la résistance.

### Partie 3 : Apprentissage supervisé

1. Mettre de côté environ 25 % des observations comme jeu de test.
2. Comparer exactement ces sept apprenants dans un seul appel à
   `funcml::compare()` :

   ```r
   c("glm", "rpart", "ranger", "fastgbm", "densemlp", "kknn", "e1071_svm")
   ```

   Rapporter `auc`, `brier`, `accuracy` et `f1`, avec leurs incertitudes, puis
   expliquer l'apport et les limites de l'exactitude et du F1-score.
3. Choisir le modèle sur l'AUC validée croisée, régler ses hyperparamètres,
   puis l'ajuster sur tout le jeu d'entraînement.
4. Comparer la validation standard à une validation groupée par pays avec
   `group_cv()`.
5. Évaluer une seule fois le modèle final sur le jeu de test : ROC/AUC, Brier,
   calibration, importance globale et explication locale.
6. Inclure le phénotype non supervisé comme variable candidate et juger sa
   contribution.

### Partie 4 : Lecture croisée et hypothèses

Répondre explicitement à la question suivante : le phénotype non supervisé
apporte-t-il une information au modèle supervisé ? Appuyer la réponse sur
l'importance par permutation ou les valeurs SHAP.

Discuter ensuite, en deux à quatre phrases chacune :

- **H1, co-résistance** : association avec les autres antibiotiques et
  plausibilité biologique ;
- **H2, signal d'assemblage** : rôle de `genome_length` et `gc_content`,
  direct ou indirect ;
- **H3, confusion géographique/hôte** : effet de collecte, diffusion clonale
  ou combinaison des deux.

## Garde-fous méthodologiques

Le rapport doit rendre visibles les cinq contrôles suivants :

1. **Décalage de distribution** : comparer la composition des ensembles
   d'entraînement et de test.
2. **Dépendance entre exemples** : utiliser une validation groupée adaptée.
3. **Confusion** : discuter la structure de population, le pays et l'hôte.
4. **Fuite de données** : ajuster l'imputation, le regroupement ou toute
   sélection sur l'entraînement uniquement et, lorsque nécessaire, dans les
   plis de validation croisée.
5. **Déséquilibre des classes** : interpréter conjointement discrimination,
   calibration et performances par classe.

## Rendu du rapport

Depuis le dossier `project/`, rendre le document avec :

```bash
quarto render projet-nom-prenom.qmd --to html
```

Pour produire un PDF lorsque LaTeX est installé :

```bash
quarto render projet-nom-prenom.qmd --to pdf
```

Le rendu doit être relancé depuis une session R propre avant la remise. Une
analyse qui fonctionne uniquement grâce à des objets préexistants dans
l'environnement R n'est pas reproductible.

## Livrables

Remettre :

- le fichier source Quarto ou R Markdown (`.qmd` ou `.Rmd`) ;
- l'export HTML ou PDF produit à partir de ce fichier ;
- une conclusion d'environ une demi-page répondant à H1–H3 et à la lecture
  croisée ;
- tous les fichiers supplémentaires indispensables à l'exécution, s'il y en
  a.

Convention de nommage conseillée :

```text
projet-nom-prenom.qmd
projet-nom-prenom.html
```

Ne pas remettre l'historique R (`.Rhistory`), le dossier `.Rproj.user/`, le
cache Quarto ou les sorties intermédiaires.

## Checklist avant remise

- [ ] Le document repart d'une session R vide et charge lui-même ses packages.
- [ ] La graine aléatoire est fixée avec `set.seed()`.
- [ ] La cible est absente du clustering non supervisé.
- [ ] Le jeu de test n'intervient dans aucun choix de modèle ou hyperparamètre.
- [ ] Les transformations apprises sur les données sont confinées aux plis
      d'entraînement appropriés.
- [ ] Les sept apprenants demandés sont présents dans `compare()`.
- [ ] AUC, Brier, accuracy et F1 sont rapportés pour la comparaison.
- [ ] La validation groupée par pays est comparée à la validation standard.
- [ ] ROC/AUC, Brier et calibration sont rapportés sur le test final.
- [ ] Une explication globale et une explication locale sont incluses.
- [ ] La contribution du phénotype est discutée explicitement.
- [ ] H1, H2 et H3 reçoivent chacune une réponse fondée sur les résultats.
- [ ] Les limites sont rappelées : pas de preuve causale et pas de validation
      externe.
- [ ] Le document complet est rendu sans erreur depuis une session propre.

## Dépannage

### Le fichier CSV est introuvable

Vérifier le répertoire courant :

```r
getwd()
file.exists("ecoli_amr_isolate_v2.csv")
```

Ouvrir `project/project.Rproj` ou placer le document dans `project/`. Éviter
les chemins absolus propres à un ordinateur.

### Un package est introuvable

Relancer l'installation avec `pak::pak()` puis redémarrer la session R. Les
packages du cours sont des versions de développement installées depuis GitHub.

### Le rendu PDF échoue

Commencer par rendre en HTML. Pour le PDF, installer une distribution LaTeX ou
TinyTeX :

```bash
quarto install tinytex
```

### Un modèle produit un avertissement

Un avertissement n'est pas toujours une erreur. Lire son contenu : un niveau
catégoriel absent d'un pli, une variable constante ou une séparation parfaite
peuvent signaler une propriété importante des données. Ne pas masquer les
avertissements sans en comprendre la cause.

### Le calcul SHAP est lent

SHAP repose sur des simulations et peut prendre plusieurs minutes. Commencer
par peu d'observations et une petite valeur de `nsim` pendant le développement,
puis augmenter ces valeurs pour le rendu final.

## Utilisation des références

Les articles et la bibliographie se trouvent dans `refs/`. Toute affirmation
biologique ou méthodologique importante doit être rattachée à une source. Le
rapport doit distinguer clairement :

- une association prédictive ;
- une hypothèse biologique plausible ;
- une conclusion causale, qui n'est pas établie par ce projet observationnel.

## Aide

Avant de demander de l'aide, conserver le message d'erreur complet et noter :

- la commande exécutée ;
- le répertoire courant (`getwd()`) ;
- la version de R (`R.version.string`) ;
- le résultat de `sessionInfo()`.

Ces informations permettent de reproduire et diagnostiquer rapidement le
problème.
