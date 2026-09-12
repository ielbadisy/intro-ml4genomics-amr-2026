set.seed(2026)

# Correspondance avec les questions de projet-enonce.qmd
#
# P1.1–P1.4 : Prétraitement
# P2.1–P2.4 : Apprentissage non supervisé (phynotype)
# P3.1–P3.6 : Apprentissage supervisé (funcml)
# P4         : Question de lecture croisée
# H1–H3      : Hypothèses de travail
# G1–G5      : Garde-fous méthodologiques
#
# Chaque bloc reprend le code de projet-correction.qmd ; les commentaires
# donnent la réponse/interprétation correspondante sans prose Quarto.

# ---- Données et contrôles préliminaires (contexte) ----

required_packages <- c("basetable", "missknn", "phynotype", "funcml")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Packages manquants : ", paste(missing_packages, collapse = ", "),
    ". Installez-les avec pak avant d'exécuter la correction."
  )
}

library(basetable)

# Fonctionne depuis la racine du dépôt comme depuis le dossier project/.
data_candidates <- c(
  "ecoli_amr_isolate_v2.csv",
  file.path("project", "ecoli_amr_isolate_v2.csv")
)
data_path <- data_candidates[file.exists(data_candidates)][1]
if (is.na(data_path)) {
  stop("Fichier ecoli_amr_isolate_v2.csv introuvable (racine ou project/).")
}

raw <- btread(data_path)
dim(raw)

preview(raw)
count(raw, by = "ciprofloxacin")

# ============================================================================
# Partie 1 : Prétraitement
# ============================================================================

# ---- P1.1 : Construire la cible binaire cipro_R ----

library(missknn)

abx_other <- c("ampicillin", "cefotaxime", "meropenem", "sxt",
               "ceftazidime", "gentamicin", "pip_taz")

iso <- raw |>
  uniquerows(cols = "genome_id", .keep_all = TRUE) |>
  subset(ciprofloxacin %in% c("Resistant", "Susceptible")) |>
  transform(
    cipro_R = factor(ifelse(ciprofloxacin == "Resistant", "R", "S"),
                     levels = c("S", "R"))
  )

# ---- P1.2 : Distinguer AST non testé et métadonnée manquante ----

# isolation_country : métadonnée imputable (missknn), contrairement aux AST
country_imp <- iso[, c("isolation_country", "genome_length", "gc_content",
                        "contigs", "collection_year")]
country_imp$isolation_country <- factor(country_imp$isolation_country)
imp <- missknn(country_imp, k = 5, seed = 2026)
iso$isolation_country <- complete(imp)$isolation_country

# co-résistance : NA = "non testé", catégorie structurelle, non imputée
for (a in abx_other) iso[[a]] <- factor(nato(iso[[a]], "Not_tested"))

# Réponse : l'AST absent signifie "non testé" et reste une catégorie
# explicite. Le pays manquant est une métadonnée incomplète, donc imputable
# par voisins proches avec missknn sans fabriquer de résultat clinique.

# ---- P1.3 : Regrouper les modalités rares (pays et hôte) ----

# pays rares regroupés ("Other" = rare, pas manquant) pour éviter des niveaux absents d'un pli de CV

frequent_countries <- names(which(table(iso$isolation_country) >= 25))
iso$country_grp <- factor(ifelse(
  iso$isolation_country %in% frequent_countries,
  as.character(iso$isolation_country), "Other"
))

iso$host_grp <- factor(ifelse(
  is.na(iso$host_name), "Unknown",
  ifelse(grepl("Homo sapiens", iso$host_name), "Human", "Animal/other")
))

iso <- pick(iso, c("genome_id", "cipro_R",
                   "genome_length", "gc_content", "contigs",
                   abx_other, "host_grp", "country_grp"))
iso <- iso[complete.cases(iso), ]

dim(iso)
count(iso, by = "cipro_R")

# Réponse P1.1/P1.3 : 683 isolats sans statut intermédiaire sont conservés.
# Les pays avec moins de 25 isolats deviennent "Other" ; l'hôte est ramené
# à Human / Animal-other / Unknown afin de stabiliser les niveaux dans les plis.

# ---- P1.4 : Prévalence et choix des métriques ----

prev <- mean(iso$cipro_R == "R")
round(prev, 3)

# Réponse : 144/683 = 21,1 % des isolats sont résistants. La cible est
# déséquilibrée : accuracy seule est trompeuse ; compléter par AUC, Brier,
# calibration, F1 et sensibilité de la classe R.

# ============================================================================
# Partie 2 : Apprentissage non supervisé (phynotype)
# ============================================================================

# ---- Installation et chargement de phynotype ----

library(phynotype)

# ---- P2.1 : Données mixtes, sans la cible cipro_R ----

cluster_vars <- c("genome_length", "gc_content", "contigs",
                   abx_other, "host_grp", "country_grp")

iso_mixed <- prepare_mixed_data(iso[, cluster_vars], center = TRUE, scale = TRUE)
dim(iso_mixed)

# ---- P2.2 : Clustering, validation et choix de k ----

fit_phen <- cluster(iso_mixed, method = "kmeans", k = 3, seed = 2026)
fit_phen
sizes(fit_phen)

val <- validate(fit_phen)
val$metrics_table

grid_val <- validate(iso_mixed, method = "kmeans", k = 2:6)
grid_val$metrics_table

# Stabilité et accord entre solutions candidates.
mfit <- metacluster(
  iso_mixed,
  methods = c("kmeans", "pam", "hclust"),
  k = 2:5,
  seed = 2026
)
mfit$stability_summary
mfit$selection_summary

plot_coassoc(mfit)
plot_consensus(mfit)
plot_dendrogram(mfit)

# Réponse : k = 3 est retenu. Sa silhouette est modeste (0,155), mais son
# ARI bootstrap est très élevé (0,976/0,981). k = 4 ou 6 gagne à peine en
# silhouette, tandis que la stabilité chute au-delà de k = 3.

# ---- P2.3 : Profils, tailles, biplot et interprétabilité ----

expl <- explore(fit_phen, data = iso[, cluster_vars])
names(expl)

plot_clusters(fit_phen, data = iso_mixed)

plot_silhouette(fit_phen, data = iso_mixed)

plot_feature_profiles(expl)
expl$separation_table

plot_cluster_sizes(fit_phen)

plot_biplot(fit_phen)
plot_biplot(fit_phen, variant = "cos2")

imp_clust <- feature_importance(fit_phen, data = iso_mixed, metric = "silhouette",
                                 n_repeats = 10, seed = 2026)
imp_clust$summary
plot(imp_clust)

cp <- ceteris_paribus(
  fit_phen, iso_mixed[1:3, ],
  features = c("genome_length", "gc_content"),
  grid_size = 20
)
plot(cp)

lx <- lime_explain(fit_phen, iso_mixed[1:3, ], n_permutations = 100,
                    n_features = 6, seed = 2026)
lx$explanations[, c("observation", "feature", "estimate", "direction")]
plot(lx)

# Réponse P2.3 : la géométrie est surtout portée par la collecte :
# host_grp_Unknown et country_grp_Norway dominent l'importance (~0,030),
# suivis de host_grp_Human (~0,027), England (~0,018) et UK (~0,017).
# Les co-résistances sont intermédiaires ; longueur et GC sont faibles.

# ---- P2.4 : Association a posteriori entre phénotype et cipro_R ----

new_isolates <- iso_mixed[1:5, ]
pred_phen <- predict(fit_phen, new_isolates)
pred_phen$clusters
head(pred_phen$distances)

iso$phenotype <- factor(clusters(fit_phen))
count(iso, by = c("phenotype", "cipro_R"))

# Réponse : le taux R varie de 11,4 % (phénotype 1) à 38,0 % (phénotype 2)
# et 27,0 % (phénotype 3). Le phénotype 2 concentre donc plus de trois fois
# plus de résistance que le phénotype 1, sans que cipro_R ait servi au cluster.

# ============================================================================
# Partie 3 : Apprentissage supervisé (funcml)
# ============================================================================

# ---- Installation et chargement de funcml ----

library(funcml)

# ---- P3.1 : Mettre de côté environ 25 % des données ----

set.seed(2026)
idx_test <- sample(nrow(iso), floor(0.25 * nrow(iso)))
train <- iso[-idx_test, ]
test  <- iso[ idx_test, ]
c(train = nrow(train), test = nrow(test))

# ---- P3.2 : Comparer exactement sept apprenants par CV ----

f <- cipro_R ~ genome_length + gc_content + contigs +
     ampicillin + gentamicin + cefotaxime + meropenem + sxt +
     ceftazidime + pip_taz + host_grp + country_grp +
     phenotype # P3.6 : le phénotype est bien une variable candidate

candidate_models <- c("glm", "rpart", "ranger", "kknn", "e1071_svm",
                       "fastgbm", "densemlp")

cmp <- compare(
  data = train, formula = f,
  models = candidate_models,
  resampling = cv(5, seed = 2026),
  metrics = c("auc", "brier", "accuracy", "f1"),
  seed = 2026
)
cmp$results
plot(cmp)

# Réponse : ranger domine (AUC 0,889 [0,877 ; 0,900], Brier 0,207,
# accuracy 0,864, F1 0,782). Le F1 macro révèle mieux une faiblesse sur la
# classe minoritaire que l'accuracy, mais ne remplace pas la sensibilité.

# ---- P3.3 : Sélection sur l'AUC, réglage et ajustement final ----

auc_tbl <- subset(cmp$results, metric == "auc")
best_model <- auc_tbl$model[which.max(auc_tbl$mean)]
best_model

tuning_grid <- switch(best_model,
  ranger    = expand.grid(num.trees = c(300, 500, 800),
                           mtry = c(2, 4, 6),
                           min.node.size = c(1, 5, 10)),
  rpart     = expand.grid(cp = c(0.001, 0.01, 0.05),
                           minsplit = c(10, 20, 40)),
  kknn      = expand.grid(k = c(5, 7, 11, 15), distance = c(1, 2)),
  e1071_svm = expand.grid(cost = c(0.1, 1, 10), gamma = c(0.01, 0.1, 1)),
  fastgbm   = expand.grid(ntrees = c(100, 200, 400),
                           learning_rate = c(0.05, 0.1, 0.2),
                           max_depth = c(3, 5, 7)),
  densemlp  = expand.grid(dropout = c(0, 0.2), lr = c(1e-3, 5e-4)),
  glm       = NULL
)

if (!is.null(tuning_grid)) {
  tuned <- tune(
    data = train, formula = f, model = best_model,
    grid = tuning_grid,
    resampling = cv(5, seed = 2026),
    metric = "auc", type = "prob"
  )
  tuned
  tuned$best
  plot(tuned)
  fit_best <- tuned$fit_best
} else {
  fit_best <- fit(f, data = train, model = best_model, seed = 2026)
}

fit_best

# Réponse : ranger est retenu sur l'AUC moyenne et son intervalle ; kknn et
# glm restent les concurrents proches. La grille règle ensuite uniquement
# le meilleur apprenant avant l'ajustement sur tout le jeu d'entraînement.

# ---- P3.4 / G2 : Validation croisée groupée par pays ----

gcv <- group_cv(v = 5, group = train$country_grp, seed = 2026)
ev_grouped <- evaluate(
  data = train, formula = f, model = best_model, spec = fit_best$spec,
  resampling = gcv, metrics = c("auc", "brier", "accuracy", "f1"), seed = 2026
)
ev_grouped$summary
plot(ev_grouped)

# Réponse : l'AUC groupée par pays tombe à 0,714 [0,567 ; 0,860], contre
# 0,889 en CV standard ; accuracy = 0,691, F1 = 0,668 et Brier = 0,412.
# Cette forte baisse indique un raccourci lié à la collecte (H3).

# ---- P3.5 : Évaluation et interprétation sur le jeu de test ----

p_mat <- predict(fit_best, newdata = test, type = "prob")
head(p_mat)

p_R <- p_mat[, "R"]
cls <- predict(fit_best, newdata = test)
table(predicted = cls, observed = test$cipro_R)

c(accuracy = accuracy(test$cipro_R, cls),
  f1       = f1(test$cipro_R, cls))

# Réponse : 135 vrais S, 13 vrais R, 1 faux R et 21 résistances manquées.
# Accuracy = 0,871 et F1 macro = 0,778, mais la sensibilité R n'est que
# 13/34 = 38,2 %, contre une spécificité de 135/136 = 99,3 %.

# Importance globale par permutation.

vip <- interpret(fit_best, data = train, method = "permute", seed = 2026)
vip$result$scores
plot(vip)

# Réponse : gc_content arrive en tête, puis country_grp, genome_length,
# gentamicin, sxt, host_grp et contigs. phenotype est avant-dernier : la
# contribution du regroupement non supervisé est réelle mais secondaire.

plot(interpret(fit_best, data = train, method = "pdp", features = "gc_content"))
plot(interpret(fit_best, data = train, method = "pdp", features = "genome_length"))

# SHAP local sur une observation, puis SHAP global.
shap_local <- interpret(
  fit_best, data = train, method = "shap",
  newdata = test[1, , drop = FALSE], nsim = 30, seed = 2026
)
print(shap_local)

plot(shap_local, kind = "waterfall")

shap_global <- interpret(
  fit_best, data = train, method = "shap",
  newdata = test, nsim = 10, seed = 2026
)
plot(shap_global, kind = "bar")

plot(shap_global, kind = "beeswarm")
plot(shap_global, kind = "dependence", v = "gc_content")

# Réponse SHAP globale : country_grp (~0,042) et gc_content (~0,039)
# dominent ; phenotype contribue ~0,014, devant pip_taz et meropenem seulement.

# Calibration.

calib <- interpret(fit_best, data = test, method = "calibration")
summary(calib)
plot(calib)

# Analyse de courbe de décision.
dca_res <- interpret(fit_best, data = test, method = "dca")
plot(dca_res)

# Réponse : le modèle garde un bénéfice net positif sur une large plage
# de seuils et dépasse nettement "traiter tout le monde" au-delà de ~0,20.

# Courbe ROC et AUC.

roc_obj <- roc_curve(test$cipro_R, p_R)
cat("AUC test :", round(roc_obj$auc, 3), "\n")
plot(roc_obj)

# Réponse : AUC test = 0,814, cohérente avec la CV standard (0,889), mais
# plus proche de la validation groupée (0,714). Calibration imparfaite :
# ECE = 0,104 et MCE = 0,251 ; les probabilités demanderaient recalibration.

# ---- P3.6 : Vérifier l'inclusion du phénotype ----

"phenotype" %in% attr(terms(f), "term.labels")

# ============================================================================
# Partie 4 : Question de lecture croisée
# ============================================================================

# ---- P4 : Le phénotype contribue-t-il au modèle supervisé ? ----

# Réponse : phenotype apporte un signal réel mais modeste. Il est
# avant-dernier en permutation et contribue environ 0,014 en SHAP moyen,
# derrière les variables de co-résistance, d'assemblage et de collecte.

# ============================================================================
# Hypothèses de travail
# ============================================================================

# ---- H1 : Co-résistance ----
# Réponse : H1 est soutenue par le rang intermédiaire/élevé de sxt,
# gentamicin, ceftazidime, ampicillin et cefotaxime. L'effet ajusté du
# céfotaxime est estimé à +0,314 [0,299 ; 0,330], compatible avec une
# co-sélection sur éléments mobiles, sans constituer une preuve causale.

d_cef <- iso |>
  subset(cefotaxime %in% c("Resistant", "Susceptible")) |>
  transform(cefotaxime_AST = factor(as.character(cefotaxime),
                                    levels = c("Susceptible", "Resistant")))

ate <- estimate(
  data = d_cef,
  formula = cipro_R ~ genome_length + gc_content + contigs + cefotaxime_AST + host_grp,
  model = best_model,
  treatment = "cefotaxime_AST",
  estimand = "ATE",
  treatment_level = "Resistant",
  control_level = "Susceptible",
  seed = 2026
)
ate
plot(ate)

# ---- H2 : Signal d'assemblage ----
# Réponse : H2 est soutenue avec réserve. gc_content arrive en tête et
# genome_length au milieu du classement, mais leurs PDP sont modestes et
# sans direction mécanistique claire : signal probablement indirect.

# ---- H3 / G3 : Confusion géographique et liée à l'hôte ----
# Réponse : H3 est fortement soutenue par le clustering dominé par pays/hôte,
# le rang élevé de country_grp et la chute d'AUC de 0,889 à 0,714 en CV
# groupée. L'effet peut mêler biais de collecte et diffusion clonale réelle.

# ============================================================================
# Garde-fous méthodologiques à expliciter dans le rapport
# ============================================================================

# G1 : Décalage de distribution : comparer la composition train/test.
# G2 : Dépendance : voir la CV groupée par pays en P3.4.
# G3 : Confusion : voir H3 et les importances pays/hôte.
# G4 : Fuite : ajuster imputation et sélection dans les plis d'entraînement.
# G5 : Déséquilibre : voir P1.4, P3.2 et les métriques de P3.5.

# ---- Conclusion du corrigé ----
# H1, H2 et surtout H3 sont soutenues. Le phénotype contribue réellement,
# mais modestement. Ces associations ne prouvent pas un mécanisme causal
# et leur généralisation exige une validation externe.

if ("package:basetable" %in% search()) {
  detach("package:basetable", unload = TRUE, force = TRUE)
}

message("Analyse terminée avec succès.")
