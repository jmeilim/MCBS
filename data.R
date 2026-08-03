# 1. Installation de PharmacoGx (si ce n'est pas déjà fait)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("PharmacoGx", quietly = TRUE))
  BiocManager::install("PharmacoGx")

library(PharmacoGx)


gdsc <- downloadPSet("GDSC1000")


X_se <- molecularProfiles(gdsc, "rna")
X_brut <- t(assays(X_se)$exprs) 


mut_se <- molecularProfiles(gdsc, "mutation")
mut_brut <- t(assays(mut_se)$exprs)
A_brut <- mut_brut[, "BRAF"] 


sens_matrix <- summarizeSensitivityProfiles(gdsc, sensitivity.measure = "aac_recomputed")
Y_brut <- sens_matrix["Vemurafenib", ]


cellules_communes <- intersect(rownames(X_brut), names(Y_brut))
cellules_communes <- intersect(cellules_communes, names(A_brut))

X_aligne <- X_brut[cellules_communes, ]
A_aligne <- A_brut[cellules_communes]
Y_aligne <- Y_brut[cellules_communes]


valides <- !is.na(A_aligne) & !is.na(Y_aligne)

X <- X_aligne[valides, ]
A <- A_aligne[valides]
Y <- Y_aligne[valides]


cat( nrow(X),  ncol(X) )
print(table(A))
