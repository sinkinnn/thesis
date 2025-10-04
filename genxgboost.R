install.packages("xgboost")
install.packages("caret")
library("data.table")
library(xgboost)
library(caret)
library(data.table)

# 1. Simuloidaan genotyyppi- ja fenotyyppidata
set.seed(123)
n <- 200 # yksilöitä
p <- 1000 # SNP:iä

geno <- matrix(rbinom(n * p, 2, 0.3), nrow = n, ncol = p)
colnames(geno) <- paste0("SNP", 1:p)
pheno <- rnorm(n, mean = geno[,1]*0.5 + geno[,2]*0.3, sd = 1)  # fenotyyppi riippuu kahdesta SNP:stä

# 2. Jaetaan data
train_idx <- createDataPartition(pheno, p = 0.8, list = FALSE)
X_train <- geno[train_idx, ]
X_test <- geno[-train_idx, ]
y_train <- pheno[train_idx]
y_test <- pheno[-train_idx]

# 3. XGBoost-malli
dtrain <- xgb.DMatrix(data = X_train, label = y_train)
dtest <- xgb.DMatrix(data = X_test, label = y_test)

params <- list(objective = "reg:squarederror", max_depth = 5, eta = 0.05)
model <- xgb.train(params = params, data = dtrain, nrounds = 200, verbose = 0)

preds <- predict(model, dtest)
real_mse <- mean((preds - y_test)^2)

# 4. Permutaatio: satunnaistetaan fenotyyppi ja toistetaan
n_perm <- 100
perm_mse <- numeric(n_perm)

for (i in 1:n_perm) {
  y_perm <- sample(y_train)
  dtrain_perm <- xgb.DMatrix(data = X_train, label = y_perm)
  model_perm <- xgb.train(params = params, data = dtrain_perm, nrounds = 200, verbose = 0)
  preds_perm <- predict(model_perm, dtest)
  perm_mse[i] <- mean((preds_perm - y_test)^2)
}

# 5. P-arvo: kuinka moni permutaatio oli yhtä hyvä tai parempi
pval <- mean(perm_mse <= real_mse)
cat("Todellinen MSE:", round(real_mse, 3), "\n")
cat("Permutaatioihin perustuva p-arvo:", round(pval, 4), "\n")