library(readxl)

path_glio <- "/Users/jmeilim/Desktop/MCBS/Radiomics_paper/Gliosarcoma/GBM_GBS.xlsx"


Gdat <- as.data.frame(read_excel(path_glio, skip = 3))
cat("Dimensions brutes:", nrow(Gdat), "x", ncol(Gdat), "\n")
cat("6 premières colonnes:", paste(names(Gdat)[1:6], collapse=", "), "\n\n")


clin_cols  <- intersect(c("name","Sex","Age","Localization","Necrosis","Edema"), names(Gdat))
feat_names <- setdiff(names(Gdat), clin_cols)


X_raw <- as.matrix(sapply(Gdat[, feat_names], function(z) as.numeric(as.character(z))))
rownames(X_raw) <- Gdat[["name"]]


ok_col <- apply(X_raw, 2, function(z) all(is.finite(z)) && sd(z) > 0)
X_raw  <- X_raw[, ok_col, drop = FALSE]
cat("Matrice X finale:", nrow(X_raw), "patients x", ncol(X_raw), "features\n")


sds <- apply(X_raw, 2, sd)
cat("Ratio sd_max/sd_median:", round(max(sds)/median(sds), 1), "\n")