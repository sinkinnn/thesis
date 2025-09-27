# Asenna ja lataa simstudy-paketti
install.packages("simstudy")
library(simstudy)

# Määritellään BMI ja markkereiden simulaatio
# BMI:n arvot perustuvat sinun aiemmin lähettämiisi esimerkkeihin
mean_bmi <- mean(c(0.3762641, 0.4699584, 0.4298125, 0.4691828, 0.3435240, 0.5079082))
sd_bmi <- sd(c(0.3762641, 0.4699584, 0.4298125, 0.4691828, 0.3435240, 0.5079082))

# Simuloidaan data
def_bmi <- defData(varname = "BMI", dist = "normal", formula = mean_bmi, variance = sd_bmi^2)
# Simuloidaan myös muita markkereita, esimerkiksi glucose ja insulin
def_glucose <- defData(varname = "glucose", dist = "normal", formula = 100, variance = 20)
def_insulin <- defData(varname = "insulin", dist = "normal", formula = 15, variance = 5)

# Generoidaan 1000 yksilöä
sim_data <- genData(1000)
sim_data <- addColumns(def_bmi, sim_data)
sim_data <- addColumns(def_glucose, sim_data)
sim_data <- addColumns(def_insulin, sim_data)

# Tarkastellaan ensimmäisiä rivejä simuloidusta datasta
head(sim_data)

# Piirretään histogrammi BMI-arvoista
library(ggplot2)
ggplot(sim_data, aes(x = BMI)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "Simuloidut BMI-arvot", x = "BMI", y = "Frekvenssi") +
  theme_minimal()