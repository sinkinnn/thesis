# Pakettien asentaminen
install.packages("BGLR")
install.packages("mice")
install.packages("ggplot2")
install.packages("olsrr")
install.packages("poolr")
install.packages("fitdistrplus")
install.packages("extRemes")
install.packages("evd")
install.packages("qqman")
install.packages("DescTools")  # Jos paketti ei ole jo asennettu
library(DescTools)
library(qqman)
library(mice)
library(evd)
library(extRemes)
library(fitdistrplus)
library(poolr)
library(ggplot2)
library(BGLR)  
library(dplyr)
library(olsrr)
install.packages("GGally")
library(GGally)

data(package="BGLR")

# 3. Lataa erityisesti 'mice' datasetti
data(mice)

# 4. Tarkastele 'mice.pheno' dataa
head(mice.pheno)

BMI <- mice.pheno$Obesity.BMI
str(BMI)


install.packages("bigPSD", repos="http://R-Forge.R-project.org")
install.packages("BiocManager")
BiocManager::install("snpStats")

ls()
?ls


dim(mice.X)
dim(mice.pheno)
head(mice.pheno)
??mice.pheno

print(bmi)

summary(data(mice))
########### Vastemuuttujan arvojen selvittäminen (obesity.BMI)

# Lasketaan BMI käyttäen obesityEndnormalBW ja bodylength muuttujia
BMI_calculated <- mice.pheno$Obesity.EndNormalBW / (mice.pheno$Obesity.BodyLength ^ 2)
head(BMI_calculated)
data(mice)
mice.pheno$Obesity.BMI <- BMI_calculated

#Laske p-arvojen jakauma
#Kun saat selville todelliset positiiviset, niin siinä sulla on 10% DE geenit ilman korjausmenetelmiä
LDL <- mice.pheno

# esimerkillisiä BMI-arvoja hiiridatasta
bmi_values <- c(0.3762641, 0.4699584, 0.4298125, 0.4691828, 0.3435240, 0.5079082)

bmi <- mice.pheno$Obesity.BMI
mice_df <- as.data.frame(mice)
print(mice_df)
head(mice_df)
head(BMI_calculated)

summary(bmi)
sd(bmi)


# Keskiarvo
mean_bmi <- mean(bmi, na.rm = TRUE)
cat("Keskiarvo (BMI):", mean_bmi, "\n")

# Keskihajonta
sd_bmi <- sd(BMI_calculated, na.rm = TRUE)
cat("Keskihajonta (BMI):", sd_bmi, "\n")

# Varianssi
var_bmi <- var(BMI_calculated, na.rm = TRUE)
cat("Varianssi (BMI):", var_bmi, "\n")

nrow(mice.X)


# Simuloidaan BMI-arvot normaalijakaumasta käyttäen laskettua keskiarvoa ja hajontaa
set.seed(123)  # Varmistaa, että tulokset ovat toistettavissa
simulated_bmi <- rnorm(1000, mean = mean_bmi, sd = sd_bmi)  # Luodaan 1000 simuloitua BMI-arvoa

# Voimme tarkastella simuloitua dataa
summary(simulated_bmi)

# Histogrammi simuloidusta datasta
hist(simulated_bmi, main = "Simuloidut BMI-arvot", xlab = "BMI", col = "skyblue", border = "black")
hist(BMI_calculated)

################

# Luo dataframe tuloksia varten
manhattan_data <- data.frame(
  SNP = 1:ncol(mice.X),
  p_value = p_values,
  neg_log10_p = -log10(p_values)
)
print(manhattan_data)


# Suurin p-arvo
min_value <- min(manhattan_data$p_value)

# Rivi, jossa suurin p-arvo esiintyy
min_index <- which.min(manhattan_data$p_value)

# SNP-tunniste kyseiseltä riviltä
snp_id <- manhattan_data$SNP[min_index]

print(min_value)
print(snp_id)





# Manhattan plot -log10(p-arvo) jokaiselle SNP:lle, graduun 
p <- plot(-log10(p_values), main="Manhattan kuvio ylipainoisten hiirien BMI", xlab="SNP Markkerin tunniste", ylab="-log10(p-arvo)")
abline(h = -log10(0.0001), col="pink") #rotta
abline(h = -log10(0.05), col = "red")  # Merkitsevä taso
abline(h = -log10((0.05/10346)), col = "blue") #Bonferroni korjattu taso, ihan ok, usein tiukka
abline(h = -log10((0.05/2816)), col = "green") #PCA, simpleM Bonferroni korjattu taso
ggsave("kuvio.pdf", plot = p, width = 6, height = 4)

############# Permutointi

set.seed(123)  # Toistettavuus
n <- length(p_values)
group_labels <- sample(c("A", "B"), size = n, replace = TRUE)

group_A <- p_values[group_labels == "A"]
group_B <- p_values[group_labels == "B"]

# Tarkastellaan dimensioita
cat("Group A:", length(group_A), "\n")
cat("Group B:", length(group_B), "\n")

median_val <- median(p_values)
group_low <- p_values[p_values <= median_val]
group_high <- p_values[p_values > median_val]

cat("Low group:", length(group_low), "\n")
cat("High group:", length(group_high), "\n")

# Alkuperäinen ero ryhmien keskiarvoissa
observed_diff <- mean(group_low) - mean(group_high)

# Permutointien määrä
n_perm <- 10000
perm_diffs <- numeric(n_perm)

# Yhdistetään kaikki p-arvot
all_values <- c(group_low, group_high)
n_low <- length(group_low)

# Permutointisilmukka
set.seed(42)
for (i in 1:n_perm) {
  permuted <- sample(all_values)
  perm_low <- permuted[1:n_low]
  perm_high <- permuted[(n_low + 1):length(permuted)]
  perm_diffs[i] <- mean(perm_low) - mean(perm_high)
}

# Laske permutointip-arvo
p_perm <- mean(abs(perm_diffs) >= abs(observed_diff))
cat("Permutointip-arvo:", p_perm, "\n")

hist(perm_diffs, breaks = 50, main = "Permutointijakauma", xlab = "Erotus keskiarvoissa")
abline(v = observed_diff, col = "red", lwd = 2)



################################### TESTI123
set.seed(123)
n_perm <- 10000
n_snps <- ncol(mice.X)
perm_pvals <- numeric(n_snps)

for (snp in 1:n_snps) {
  geno <- mice.X[, snp]
  
  # Poistetaan puuttuvat arvot
  complete_idx <- complete.cases(geno, BMI)
  geno_clean <- geno[complete_idx]
  pheno_clean <- BMI[complete_idx]
  
  # Alkuperäinen ero: keskiarvo BMI genotyyppiryhmien välillä
  median_val <- median(geno_clean)
  group_low <- pheno_clean[geno_clean <= median_val]
  group_high <- pheno_clean[geno_clean > median_val]
  observed_diff <- mean(group_low) - mean(group_high)
  
  # Permutointi
  all_values <- c(group_low, group_high)
  n_low <- length(group_low)
  perm_diffs <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    permuted <- sample(all_values)
    perm_low <- permuted[1:n_low]
    perm_high <- permuted[(n_low + 1):length(permuted)]
    perm_diffs[i] <- mean(perm_low) - mean(perm_high)
  }
  
  # Empiirinen p-arvo
  perm_pvals[snp] <- mean(abs(perm_diffs) >= abs(observed_diff))
}

# Muutetaan p-arvot log-muotoon
log_pvals <- -log10(perm_pvals)

# Jos sinulla ei ole kromosomi- ja positio-tietoja, käytetään indeksiä
snp_positions <- 1:n_snps
chromosome <- rep(1, n_snps)  # Oletetaan kaikki SNP:t yhdellä kromosomilla

# Laitetaan data data.frame-muotoon
manhattan_data1 <- data.frame(
  SNP = paste0("SNP", snp_positions),
  CHR = chromosome,
  BP = snp_positions,
  P = perm_pvals
)

# Asennetaan tarvittava paketti (jos ei ole)
if (!requireNamespace("qqman", quietly = TRUE)) {
  install.packages("qqman")
}
library(qqman)

# Piirretään Manhattan-plot
manhattan(manhattan_data1, 
          col = c("skyblue", "navy"), 
          main = "Manhattan Plot - BMI vs SNPs",
          ylim = c(0, max(-log10(perm_pvals)) + 1),
          suggestiveline = -log10(1e-4), 
          genomewideline = -log10(5e-8))

# Korvataan nollat pienellä arvolla
perm_pvals[perm_pvals == 0] <- 1e-10

# Lasketaan log-p-arvot
log_pvals <- -log10(perm_pvals)

# Päivitetään ylim
ylim_val <- c(0, max(log_pvals) + 1)
ylim_val <- c(0, max(log_pvals) + 1)
# Piirretään uudelleen
manhattan_data1$P <- perm_pvals  # päivitetään p-arvot

manhattan(manhattan_data1, 
          col = c("skyblue", "navy"), 
          main = "Manhattan Plot - BMI vs SNPs",
          ylim = ylim_val,
          suggestiveline = -log10(1e-4), 
          genomewideline = -log10(5e-8))
head(manhattan_data1)

any(!is.finite(log_pvals))
all(rownames(mice.map) == colnames(mice.X))
# Korvataan nollat ja puuttuvat p-arvot pienellä arvolla
perm_pvals[!is.finite(perm_pvals) | perm_pvals == 0] <- 1e-10

# Lasketaan log-p-arvot uudelleen
log_pvals <- -log10(perm_pvals)

# Varmistetaan että kaikki ovat nyt kelvollisia
any(!is.finite(log_pvals))  # pitäisi nyt antaa FALSE

# Järjestetään perm_pvals vastaamaan mice.map:n järjestystä
snp_order <- match(rownames(mice.map), colnames(mice.X))

# Tarkistetaan että ei ole NA-arvoja
if (any(is.na(snp_order))) {
  stop("Jotkin SNP-nimet eivät täsmää mice.map ja mice.X välillä.")
}

# Järjestetään p-arvot
perm_pvals_ordered <- perm_pvals[snp_order]

head(rownames(mice.map))       # Kartan SNP-nimet
head(colnames(mice.X))         # Genotyyppimatriisin SNP-nimet









################################################### ilman ryhmittelyä
set.seed(123)
n_perm <- 10000
n_snps <- ncol(mice.X)
perm_pvals <- numeric(n_snps)

for (snp in 1:n_snps) {
  geno <- mice.X[, snp]
  
  # Poistetaan puuttuvat arvot
  complete_idx <- complete.cases(geno, BMI)
  geno_clean <- geno[complete_idx]
  pheno_clean <- BMI[complete_idx]
  
  # Alkuperäinen korrelaatio
  observed_corr <- cor(geno_clean, pheno_clean)
  
  # Permutaatio: satunnaistetaan fenotyyppi
  perm_corrs <- numeric(n_perm)
  for (i in 1:n_perm) {
    perm_pheno <- sample(pheno_clean)
    perm_corrs[i] <- cor(geno_clean, perm_pheno)
  }
  
  # Empiirinen p-arvo
  perm_pvals[snp] <- mean(abs(perm_corrs) >= abs(observed_corr))
}

print(perm_pvals)
plot(-log10(perm_pvals), main = "Manhattan kuvio")

min(p_values)

length(perm_pvals)

######################################################


threshold <- quantile(p_values, 0.01)
significant_snps <- which(p_values <= threshold)
mean(threshold)
plot(-log10(p_values), 
     main = "Manhattan kuvio – 0,1 % häntä",
     xlab = "SNP",
     ylab = "-log10(p-arvo)",
     col = "skyblue", pch = 20)

points(significant_snps, -log10(p_values[significant_snps]), 
       col = "red", pch = 19, cex = 1.2)

abline(h = -log10(threshold), col = "darkred", lty = 2)


######################################################

# Asetetaan siemen toistettavuuden varmistamiseksi
set.seed(123)

# Oletetaan, että p_values on valmiiksi olemassa oleva vektori, esim:
# p_values <- runif(1000) # (vain esimerkki)

# Määritellään kierrosten määrä
n_iter <- 10000

# Tallennetaan jokaisen permutaation pienin p-arvo
min_p_values <- numeric(n_iter)

for (i in seq_len(n_iter)) {
  # Permutoidaan fenotyyppivektori (tässä simuloidaan sekoittamalla p-arvot)
  permuted <- sample(p_values)
  
  # Tallennetaan pienin p-arvo
  min_p_values[i] <- min(permuted)
}

# Lasketaan 5 % kvantiili pienimmistä p-arvoista
threshold <- quantile(min_p_values, 0.05)

threshold

min(p_values)




########## Merkitsevät SNP:t
significance_threshold <- 0.05 / 1740  # valitaan Bonferroni-korjattu taso, simpleM 

########## Askeltava muuttujien valinta
# Valitse alkuperäiset merkittävät muuttujat
above_threshold_snps <- which(p_values <= significance_threshold)
selected_predictors <- colnames(mice.X)[above_threshold_snps]
subdata <- as.data.frame(mice.X[, above_threshold_snps])
colnames(subdata) <- selected_predictors

# Aloitetaan askeltava muuttujien valinta ja tarkistus
final_selected_predictors <- c()
remaining_predictors <- selected_predictors

# Askeltava muuttujien valinta, jossa etsitään aina paras merkittävä muuttuja
while (length(remaining_predictors) > 0) {
  # Laske p-arvot jäljellä oleville muuttujille
  p_values_model <- sapply(remaining_predictors, function(pred) {
    current_predictors <- c(final_selected_predictors, pred)
    model <- lm(ob_bmi ~ ., data = subdata[, current_predictors, drop = FALSE])
    p_value <- summary(model)$coefficients[nrow(summary(model)$coefficients), 4]
    return(p_value)
  })
  
  # Valitse paras merkittävä muuttuja (pienin p-arvo)
  best_predictor <- remaining_predictors[which.min(p_values_model)]
  
  # Lisää paras muuttuja listalle riippumatta p-arvosta
  final_selected_predictors <- c(final_selected_predictors, best_predictor)
  
  # Päivitä jäljellä olevat muuttujat, poista valittu muuttuja
  remaining_predictors <- setdiff(remaining_predictors, best_predictor)
  
  # Tulosta kierroksen valitut muuttujat
  cat("Kierroksen jälkeen valitut muuttujat:", final_selected_predictors, "\n")
  
  # Tarkista, jos enää ei ole merkittäviä muuttujia jäljellä
  above_threshold_snps <- which(p_values_model <= significance_threshold)
  if (length(above_threshold_snps) == 0) {
    break
  }
}

# Lopulliset tulostukset
cat("Lopullisesti valitut muuttujat:", final_selected_predictors, "\n")
cat("Valittujen muuttujien määrä:", length(final_selected_predictors), "\n")

# Lasketaan suhteellinen osuus kaikista muuttujista
total_predictors <- ncol(mice.X)
relative_proportion_above <- length(final_selected_predictors) / total_predictors
cat("Suhteellinen osuus lopullisista muuttujista:", relative_proportion_above * 100, "%\n")

# Piirretään alkuperäinen Manhattan-kuvio
manhattan_data <- data.frame(
  SNP = 1:ncol(mice.X),
  p_value = p_values,
  neg_log10_p = -log10(p_values)
)

# Valitse alkuperäiset merkittävät muuttujat (vihreän viivan yläpuolella)
significant_snps <- which(manhattan_data$neg_log10_p >= -log10(significance_threshold))

# Piirrä alkuperäinen Manhattan-kuvio
plot(
  manhattan_data$SNP, manhattan_data$neg_log10_p,
  pch = 19, cex = 0.5,
  xlab = "SNP Marker", ylab = "-log10(p-value)",
  main = "Askeltava muuttujien valinta"
)

# Lisää kynnykset (Bonferroni, SimpleM ja vihreä viiva)
abline(h = -log10(0.05 / ncol(mice.X)), col = "red", lty = 2)    # Bonferroni
abline(h = -log10(significance_threshold), col = "green")         # SimpleM vihreä viiva
# Tarkistetaan valittujen muuttujien indeksit
final_selected_indices <- match(final_selected_predictors, colnames(mice.X))
print(final_selected_indices)
# Tulostetaan tarkempaa tietoa
cat("Lopullisesti valitut muuttujat:", final_selected_predictors, "\n")
cat("Niiden indeksit:", sort(final_selected_indices), "\n")

# Tarkista, ovatko indeksit NA-arvoja (eli puuttuvatko ne alkuperäisestä data-aineistosta)
missing_indices <- which(is.na(final_selected_indices))
if (length(missing_indices) > 0) {
  cat("Seuraavat valitut muuttujat eivät löytyneet alkuperäisestä aineistosta:\n")
  cat(final_selected_predictors[missing_indices], "\n")
}

# Piirretään siniset pisteet vain niille, jotka löytyvät
valid_indices <- final_selected_indices[!is.na(final_selected_indices)]
points(
  valid_indices, manhattan_data$neg_log10_p[valid_indices],
  col = "blue", pch = 19, cex = 1
)

# Tarkistus: jos pisteitä puuttuu, päivitä ja varmista tiedot.
cat("Validien indeksien määrä:", length(valid_indices), "\n")
cat("-log10(p-arvot) valituille muuttujille:", manhattan_data$neg_log10_p[valid_indices], "\n")


######################