# Pakettien asentaminen
install.packages("BGLR")
install.packages("mice")
install.packages("readr")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("dagitty")
install.packages("ggdag")
install.packages("dplyr")
install.packages("pbapply")
install.packages("pheatmap")

library(dplyr)
library(pbapply)

library(mice)
library(BGLR)  
library(ggplot2)
library(tidyverse)

library(dagitty)
library(ggdag)
library(pheatmap)


library(readr)
library(dplyr)
pheno <- read_tsv("http://mtweb.cs.ucl.ac.uk/mus/www/GSCAN/HS_PHENOTYPES/Obesity.txt")
resids <- read_tsv("http://mtweb.cs.ucl.ac.uk/mus/www/GSCAN/HS_PHENOTYPES/Obesity.txt.residuals")
insulin <- read_tsv("http://mtweb.cs.ucl.ac.uk/mus/www/GSCAN/HS_PHENOTYPES/Insulin.txt")
geno <- read_tsv("http://mtweb.cs.ucl.ac.uk/mus/www/GSCAN/HS_GENOTYPES/Obesity_genotypes.txt")

phenobio <- read_tsv("http://mtweb.cs.ucl.ac.uk/mus/www/GSCAN/HS_PHENOTYPES/Biochemistry.txt")

tail(pheno)

?md
# Tarkastellaan puuttuvia arvoja
md.pattern(pheno)

unique(pheno$Obesity.BMI)

# Imputointi (oletus = menetelmä "pmm", predictive mean matching)
imputed_data <- mice(pheno, m = 5, method = 'pmm', seed = 123)

# Luodaan täydellinen data
pheno_imputed <- complete(imputed_data)

install.packages("missForest")
library(missForest)

# Imputointi
imputed <- missForest(pheno)
pheno_imputed <- imputed$ximp

glimpse(pheno_imputed)

# Esikatsele data
glimpse(pheno)
head(pheno)

colnames(pheno)

BW <- pheno$EndNormalBW
BL <- pheno$Obesity.BodyLength

BMI <- ((BW/BL)^2)
pheno$Obesity.BMI <- BMI
head(pheno$Obesity.BMI)

plot(BMI)

view(BW)
table(insulin$EndNormalBW)
ggplot(pheno, aes(x = GENDER)) +
  geom_bar(fill = "skyblue") +
  labs(title = "Sukupuolijakauma tutkimuksessa",
       x = "Sukupuoli", y = "Henkilömäärä")


pheno %>%
  group_by(GENDER) %>%
  summarise(
    n = n(),
    mean_weight = mean(EndNormalBW, na.rm = TRUE),
    sd_weight = sd(EndNormalBW, na.rm = TRUE),
    median_weight = median(EndNormalBW, na.rm = TRUE)
  )

pheno %>%
  group_by(GENDER) %>%
  summarise(
    n = n(),
    mean_weight = mean(BMI, na.rm = TRUE),
    sd_weight = sd(BMI, na.rm = TRUE),
    median_weight = median(BMI, na.rm = TRUE)
  )


ggplot(pheno, aes(x = GENDER, y = BMI, fill = GENDER)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("M" = "#1f77b4", "F" = "#e377c2")) +
  labs(
    title = "BMI sukupuolen mukaan",
    x = "Sukupuoli",
    y = "BMI"
  ) +
  theme_minimal()


# Yhdistetään
combined <- tibble(
  Obesity.BMI = pheno$Obesity.BMI,
  resid.Obesity.BMI = resids$resid.Obesity.BMI
)

head(combined)
# Poistetaan rivit, joissa on NA-arvoja
combined_clean <- combined %>% drop_na()

dim(combined_clean)

# Katsotaan yhteenveto
summary(model)


# Tulostetaan
head(combined)


ggplot(merged, aes(x = Obesity.BMI, y = resid.Obesity.BMI)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Residuaalit vs. Obesity.BMI",
       x = "Obesity.BMI",
       y = "Residuaali (resid.Obesity.BMI)") +
  theme_minimal()


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


pheno <- pheno %>%
  mutate(
    Date.StudyDay = as.Date(Date.StudyDay)  # muutetaan päivämääräksi jos ei jo ole
  )

pheno_summary <- pheno %>%
  group_by(GENDER, Date.StudyDay) %>%
  summarise(
    mean_BMI = mean(BMI, na.rm = TRUE),
    sd_BMI   = sd(BMI, na.rm = TRUE),
    n        = n(),
    se_BMI   = sd_BMI / sqrt(n)
  )

ggplot(pheno_summary, aes(x = Date.StudyDay, y = mean_BMI, color = GENDER, group = GENDER)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_BMI - se_BMI, ymax = mean_BMI + se_BMI), width = 0.2) +
  scale_color_manual(values = c("M" = "#1f77b4", "F" = "#e377c2")) +
  labs(
    title = "Keskimääräinen BMI tutkimuspäivittäin sukupuolen mukaan",
    x = "Tutkimuspäivä",
    y = "Keskimääräinen BMI ± SE",
    color = "Sukupuoli"
  ) +
  theme_minimal(base_size = 14)



dag_triangle <- dagitty('
dag {
  GENDER -> BMI
  GENDER -> GENO
  GENO -> BMI
}
')

# Asetetaan layout niin, että muodostuu kolmio
coords <- data.frame(
  name = c("GENO", "BMI", "GENDER"),
  x = c(0, 1, 0.5),
  y = c(0, 0, 1)
)

# Lisätään koordinaatit DAG-objektiin
coordinates(dag_triangle) <- coords

# Piirretään DAG kolmion muodossa
ggdag(dag_triangle) +
  theme_dag() +
  ggtitle("Kolmion muotoinen DAG: GENDER, GENO ja BMI") +
  theme(plot.title = element_text(hjust = 0.5))





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


# Manhattan plot -log10(p-arvo) jokaiselle SNP:lle, graduun 
p <- plot(-log10(p_values), main="Manhattan kuvio ylipainoisten hiirien BMI", xlab="SNP Markkerin tunniste", ylab="-log10(p-arvo)")
abline(h = -log10(0.0001), col="pink") #rotta
abline(h = -log10(0.05), col = "red")  # Merkitsevä taso
abline(h = -log10((0.05/10346)), col = "blue") #Bonferroni korjattu taso, ihan ok, usein tiukka
abline(h = -log10((0.05/2816)), col = "green") #PCA, simpleM Bonferroni korjattu taso
ggsave("kuvio.pdf", plot = p, width = 6, height = 4)


Gender <- pheno$GENDER



manhattan_data$Gender <- pheno$GENDER
manhattan_data <- rbind(
  data.frame(SNP = 1:ncol(mice.X), p_value = p_values_males, Gender = "Male"),
  data.frame(SNP = 1:ncol(mice.X), p_value = p_values_females, Gender = "Female")
)
manhattan_data$neg_log10_p <- -log10(manhattan_data$p_value)


# Suurin p-arvo
min_value <- min(manhattan_data$p_value)

# Rivi, jossa suurin p-arvo esiintyy
min_index <- which.min(manhattan_data$p_value)

# SNP-tunniste kyseiseltä riviltä
snp_id <- manhattan_data$SNP[min_index]

print(min_value)
print(snp_id)

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

# Oletetaan että mice.map on data.frame ja siinä on sarake snp_id
# Esimerkki: rs1234_C, rs5678_T jne.

# Poimi "_" jälkeinen alleeli
mice.map$allele <- sub(".*_", "", mice.map$snp_id)

# Tarkista tulos
head(mice.map$allele)

# Laske frekvenssit
allele.freq <- table(mice.map$allele)

# Muunna prosenttiosuuksiksi
allele.freq.prop <- prop.table(allele.freq)

# Tulosta
print(allele.freq)
print(round(allele.freq.prop, 3))

# Asenna tarvittaessa
install.packages("HardyWeinberg")

library(HardyWeinberg)

# Esimerkki: genotyyppien lukumäärät
genotypes <- c(AA = 30, AB = 50, BB = 20)

# Testaa HWE
HW.test <- HWChisq(genotypes)

print(HW.test)

# Alleelien lukumäärät
alleles <- c(A = 3103, C = 1351, G = 5357, T = 535)

# Laske kokonaismäärä
total_alleles <- sum(alleles)

# Alleelifrekvenssit
allele_freq <- alleles / total_alleles
round(allele_freq, 4)


# Genotyyppien frekvenssit HWE:n mukaan
genotype_freq <- outer(allele.freq, allele.freq)

# Tee symmetrisestä matriisista genotyyppien lukumäärät
# Oletetaan että yksilöitä on total_alleles / 2 (koska jokaisella on 2 alleelia)
n_individuals <- total_alleles / 2
# Luo kaikki mahdolliset genotyyppiparit
genotype_list <- expand.grid(allele1 = names(allele_freq), allele2 = names(allele_freq))

# Laske genotyyppien frekvenssit
genotype_list$freq <- mapply(function(a1, a2) {
  if (a1 == a2) {
    allele_freq[a1]^2
  } else {
    2 * allele_freq[a1] * allele_freq[a2]
  }
}, genotype_list$allele1, genotype_list$allele2)

# Luo genotyyppinimi (esim. AG tai GA → AG)
genotype_list$genotype <- apply(genotype_list[, c("allele1", "allele2")], 1, function(x) paste(sort(x), collapse = ""))

# Yhdistä duplikaatit (esim. AG ja GA)
genotype_summary <- aggregate(freq ~ genotype, data = genotype_list, sum)

# Laske lukumäärät
genotype_summary$count <- round(genotype_summary$freq * n_individuals)

# Tulosta
print(genotype_summary)

####################################################################################
# Varmistetaan että sukupuoli on faktori
pheno$GENDER <- factor(pheno$GENDER, levels = c("F", "M"))

# Vektori BMI:lle
Y <- BMI  # jos sarake on "obesity.BMI"
GENDER <- pheno$GENDER

head(GENDER)

# Funktio p-arvon laskentaan per SNP
get_p <- function() {
  fit <- lm(Y ~ mice.X + GENDER)
  summary(fit)$coefficients["geno", "Pr(>|t|)"]
}

# Luo dataframe tuloksia varten
manhattan_data <- data.frame(
  SNP = colnames(mice.X),
  p_value = p_values,
  neg_log10_p = -log10(p_values)
)

# Tulosta muutama ensimmäinen rivi
print(head(manhattan_data))

# Manhattan-kuvio
plot(manhattan_data$neg_log10_p,
     main = "Manhattan-kuvio: BMI ~ GENO + GENDER",
     xlab = "SNP Markkerin tunniste",
     ylab = "-log10(p-arvo)",
     col = ifelse(GENDER[1] == "M", "blue", "red"))

plot(manhattan_data$neg_log10_p,
     main = "Manhattan-kuvio: BMI ~ GENO + GENDER",
     xlab = "SNP Markkerin tunniste",
     ylab = "-log10(p-arvo)",
     pch = 20,
     col = "steelblue")

# Viivat merkitsevyystasoille
abline(h = -log10(0.05), col = "red")
abline(h = -log10(0.0001), col = "pink")
abline(h = -log10(0.05 / 10346), col = "blue")
abline(h = -log10(0.05 / 2816), col = "green")

# Tallennus PDF:ksi
ggsave("manhattan_gender_additive.pdf", width = 6, height = 4)


# Oletetaan, että mice.X on data.frame, jossa sarakkeet ovat SNP:t ja rivit henkilöt
# Y = BMI, GENDER = sukupuoli (factor)

# Erota miehet ja naiset
male_idx <- which(pheno$GENDER == "M")
female_idx <- which(pheno$GENDER == "F")

# P-arvot miehille
p_male <- sapply(1:ncol(mice.X), function(i) {
  fit <- lm(Y[male_idx] ~ mice.X[male_idx, i])
  summary(fit)$coefficients[2, "Pr(>|t|)"]
})

# P-arvot naisille
p_female <- sapply(colnames(mice.X), function(snp) {
  fit <- lm(Y[female_idx] ~ mice.X[female_idx, snp])
  summary(fit)$coefficients[2, "Pr(>|t|)"]
})

# Manhattan-plot
plot(-log10(p_male), col = "blue", pch = 20,
     main = "Manhattan: BMI ~ GENO sukupuolen mukaan",
     xlab = "SNP", ylab = "-log10(p-arvo)",
     ylim = range(c(-log10(p_male), -log10(p_female))))
points(-log10(p_female), col = "red", pch = 20)
legend("topright", legend = c("Miehet", "Naiset"),
       col = c("blue", "red"), pch = 20)

class(mice.X)   

dim(mice.X)          # rivit = henkilöt, sarakkeet = SNP:t
length(Y)            # täsmää rivimäärään
length(pheno$GENDER) 

summary_table

unique(pheno$GENDER)
pheno$GENDER <- toupper(trimws(pheno$GENDER))
pheno$GENDER <- factor(pheno$GENDER, levels = c("F", "M"))
levels(pheno$GENDER)
table(pheno$GENDER)
t_test_result <- t.test(BMI ~ GENDER, data = pheno)
t_test_result$conf.int


library(dplyr)

summary_table <- pheno %>%
  group_by(GENDER) %>%
  summarise(
    mean_BMI = mean(BMI, na.rm = TRUE),
    sd_BMI = sd(BMI, na.rm = TRUE),
    n = sum(!is.na(BMI))
  ) %>%
  mutate(
    se = sd_BMI / sqrt(n),
    ci_lower = mean_BMI - qt(0.975, df = n-1) * se,
    ci_upper = mean_BMI + qt(0.975, df = n-1) * se
  )

summary_table

t_test_result <- t.test(BMI ~ GENDER, data = pheno)
t_test_result$conf.int

pheno$GENDER <- toupper(trimws(pheno$GENDER))
pheno$GENDER <- factor(pheno$GENDER, levels = c("F", "M"))

summary_table <- pheno %>%
  group_by(GENDER) %>%
  summarise(
    mean_BMI = mean(BMI, na.rm = TRUE),
    sd_BMI = sd(BMI, na.rm = TRUE),
    n = sum(!is.na(BMI))
  ) %>%
  mutate(
    se = sd_BMI / sqrt(n),
    ci_lower = mean_BMI - qt(0.975, df = n-1) * se,
    ci_upper = mean_BMI + qt(0.975, df = n-1) * se
  )

summary_table


mean_F <- summary_table$mean_BMI[summary_table$GENDER == "F"]
mean_M <- summary_table$mean_BMI[summary_table$GENDER == "M"]

# Absoluuttinen ero
diff_BMI <- mean_M - mean_F

# Prosentuaalinen ero
percent_diff <- (diff_BMI / mean_F) * 100

diff_BMI
percent_diff

unique(pheno$GENDER)
table(pheno$GENDER, useNA = "ifany")


significant_snps <- subset(p_values, p_values < 5e-8)
head(significant_snps)

length(p_values)
nrow(mice.map)

head(GENDER)
mice.map$p_values <- p_values

head(sig_snps)
sig_snps <- subset(mice.map, p_values < 5e-8)
intersect_snps <- intersect(colnames(mice.X), sig_snps$snp_id)
geno_sig <- mice.X[, intersect_snps, drop = FALSE]

library(pheatmap)

pheatmap(as.matrix(geno_sig),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100)
         
significant_snps <- subset(mice.map, chr %in% c(1, 4))         

library(ggplot2)

# Laske -log10(p)
significant_snps$neglogp <- -log10(significant_snps$p_values)

tail(significant_snps)
library(ggplot2)

ggplot(significant_snps, aes(x = mbp, y = factor(chr), color = neglogp)) +
  geom_point(size = 3) +
  scale_color_gradient(low="white", high="red") +
  labs(
    x = "Genomic position (Mb)",
    y = "Chromosome",
    color = "-log10(p)",
    title = "GWAS – Chromosomes 1,4,7,15"
  ) +
  theme_minimal()