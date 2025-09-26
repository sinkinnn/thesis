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

#havaintomäärät
obs <- c(
  AA = 465,
  AC = 810,
  AG = 3213,
  AT = 321,
  CC = 88,
  CG = 1399,
  CT = 140,
  GG = 1387,
  GT = 554,
  TT = 14
)

# χ²-testi
HWChisq(obs, cc = 0)

# (vaihtoehto, tarkempi) eksakti testi
HWExact(obs)

geno <- c(465,3213,1387)
N <- sum(geno)
f <- sum(geno*c(0,1,2))/(2*N)
print(f)

hwe.prop <- c((1-f)^2, 2*f*(1-f), f^2)
rbind(obs=geno/N, hwe = hwe.prop)

# For testing HWE we use chi-square test even though counts are quite small in last cell:
hwe.test <- sum((geno-N*hwe.prop)^2/(N*hwe.prop)) # HWE test statistic
hwe.p = pchisq(hwe.test, df = 1, lower = FALSE) # P-value from the test
barplot(geno, main = paste("rs3669806; HWE P=", signif(hwe.p, 3)),
        names = c(0, 1, 2), xlab = "genotype", col = "skyblue")
