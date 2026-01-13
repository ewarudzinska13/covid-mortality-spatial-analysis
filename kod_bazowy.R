# =============================================================================
# Wczytanie wymaganych bibliotek
# =============================================================================
library(sf)
library(sp)
library(spdep)
library(spatialreg)
library(RColorBrewer)
library(classInt)
library(ggplot2)
library(dplyr)
library(readxl)
library(lmtest)
library(texreg)

# =============================================================================
# Ustawienie katalogu roboczego
# =============================================================================
setwd("C:\\Users\\ewaru\\Desktop\\e_p")
getwd()

# =============================================================================
# Wczytanie danych
# =============================================================================
dane <- read_excel("zakazenia.xlsx")

# Uzupełnienie brakujących wartości NA średnią danej kolumny (tylko dla kolumn numerycznych)
for (col in names(dane)) {
  if (is.numeric(dane[[col]])) {
    dane[[col]][is.na(dane[[col]])] <- mean(dane[[col]], na.rm = TRUE)
  }
}

# =============================================================================
# Wczytanie mapy powiatów
# =============================================================================
POW <- st_read("powiaty.shp")
POW <- st_transform(POW, 4326)  # Konwersja do układu WGS84

# Konwersja mapy do obiektu sp
pow <- as_Spatial(POW, cast = TRUE, IDs = "jpt_kod_je")
pow <- spTransform(pow, CRS("+proj=longlat +datum=NAD83"))

# =============================================================================
# Przygotowanie macierzy wag przestrzennych
# =============================================================================
# Macierz wag wg wspólnej granicy (kontyguacja)
cont.sf <- poly2nb(POW)
cont.listw <- nb2listw(cont.sf, style = "W")

# Macierz wag – odwrotna odległość (dla uzupełnienia, choć nie jest używana przy testach Morana)
crds <- st_centroid(POW)
k_val <- round(sqrt(nrow(dane)))  # liczba sąsiadów = round(sqrt(nrow(dane)))
pov.knn <- knearneigh(as.matrix(st_geometry(crds)), k = k_val)
pov.nb <- knn2nb(pov.knn)
dist <- nbdists(pov.nb, crds)
dist1 <- lapply(dist, function(x) 1/x)
dist.listw <- nb2listw(pov.nb, glist = dist1)

# =============================================================================
# Przygotowanie danych do modelu
# =============================================================================
# Ustawienie zmiennej zależnej: liczba zgonów
dane$y <- dane$zgony  

# Na podstawie intuicji usuwamy zmienne, które nie powinny być razem:
# - Zamiast liczby przypadków (x1) używamy liczby na 10 tys. mieszkańców (x2)
# - Zamiast liczby ozdrowieńców (x4) używamy liczby osób objętych kwarantanną (x5)
# - Zamiast liczby testów (x6 i x8) używamy liczby testów z wynikiem pozytywnym (x7)
# - Zamiast populacji (x14) używamy gęstości zaludnienia (x10)
# - Zamiast liczby łóżek (x15) używamy dochodu (x12)
# - Pozostawiamy mediana wieku (x13)
dane$x2  <- dane$liczba_na_10_tys_mieszkancow  
dane$x5  <- dane$liczba_osob_objetych_kwarantanna  
dane$x7  <- dane$liczba_testow_z_wynikiem_pozytywnym  
dane$x10 <- dane$ludnosc_na_1_km2  
dane$x12 <- dane$dochod  
dane$x13 <- dane$mediana_wieku.1  

# Logarytmowanie wybranych zmiennych (dla stabilności modelu)
dane$log_y <- log(dane$y + 1)
for (var in c("x2", "x5", "x7", "x10", "x12", "x13")) {
  dane[[paste0("log_", var)]] <- log(dane[[var]] + 1)
}


# =============================================================================
# Budowa modelu – regresja krokowa (stepwise selection)
# =============================================================================
full_model <- lm(log_y ~ log_x2 + log_x5 + log_x10 + log_x12 + log_x13, data = dane)
best_model <- step(full_model, direction = "backward", 
                   k = qchisq(0.01, 1, lower.tail = FALSE))
summary(best_model)

bptest(best_model)
resettest(best_model)


# =============================================================================
# Wizualizacja zmiennych użytych w modelu
# =============================================================================

library(RColorBrewer)
library(sp)

# zmienna y w logarytmie
pow$log_y <- dane$log_y
range(pow$log_y) 
rng <-  seq(min(pow$log_y), max(pow$log_y), length.out = 8) 
cls_red <- colorRampPalette(c("lightpink", "red"))(7)  # Gradient czerwony
spplot(pow, "log_y", col.regions = cls_red, at = rng, main = "Logarytm zgonów")


# zmienna x2 w logarytmie
pow$log_x2 <- dane$log_x2
range(pow$log_x2) 
rng <- seq(min(pow$log_x2), max(pow$log_x2), length.out = 8)
cls_green <- colorRampPalette(c("lightgreen", "darkgreen"))(7)
spplot(pow, "log_x2", col.regions = cls_green, at = rng, main = "Logarytm liczby zachorowań")


# zmienna x5 w logarytmie
pow$log_x5 <- dane$log_x5
range(pow$log_x5) 
rng <- seq(min(pow$log_x5), max(pow$log_x5), length.out = 8) # Podział na 7 przedziałów
cls <- brewer.pal(7, "Blues")
spplot(pow, "log_x5", col.regions = cls, at = rng, main = "Logarytm liczby osób objętych kwarantanną")


# zmienna x10 w logarytmie
pow$log_x10 <- dane$log_x10
range(pow$log_x10) 
rng <- seq(min(pow$log_x10), max(pow$log_x10), length.out = 8)  
cls_orange <- colorRampPalette(c("lightyellow", "darkorange"))(7)  # Gradient pomarańczowy
spplot(pow, "log_x10", col.regions = cls_orange, at = rng, main = "Logarytm ludności na 1 km2")


# zmienna 13 w logarytmie
pow$log_x13 <- dane$log_x13
range(pow$log_x13) 
rng <- seq(min(pow$log_x13), max(pow$log_x13), length.out = 8)  
cls_orange <- colorRampPalette(c("lavender", "purple"))(7)  # Gradient pomarańczowy
spplot(pow, "log_x13", col.regions = cls_orange, at = rng, main = "Logarytm mediany wieku")


# =============================================================================
# Dopasowanie macierzy wag do obserwacji użytych w modelu
# =============================================================================
# Pobranie indeksów obserwacji użytych w modelu
used_obs <- as.numeric(rownames(best_model$model))
# Utworzenie wektora logicznego o długości równej liczbie jednostek w cont.sf
logicalUsed <- rep(FALSE, length(cont.sf))
logicalUsed[used_obs] <- TRUE
# Wyodrębnienie podzbioru listy sąsiedztwa
cont.sf_subset <- subset.nb(cont.sf, logicalUsed)
# Utworzenie nowej macierzy wag
cont.listw_subset <- nb2listw(cont.sf_subset, style = "W")
cat("Liczba jednostek w nowej macierzy wag:", length(cont.listw_subset$neighbours), "\n")

# =============================================================================
# Testy Morana (analiza autokorelacji przestrzennej)
# =============================================================================
# Globalny test Morana dla reszt modelu OLS
moran_test1 <- lm.morantest(best_model, cont.listw_subset)
print(moran_test1)

# Klasyczny test Morana dla reszt
res <- best_model$residuals
moran_test2 <- moran.test(res, cont.listw_subset)
print(moran_test2)

# =============================================================================
# Wizualizacja reszt na mapie
# =============================================================================
brks <- c(min(res), mean(res) - sd(res), mean(res), mean(res) + sd(res), max(res))
cols <- c("steelblue4", "lightskyblue", "thistle1", "plum3")
plot(pow, col = cols[findInterval(res, brks, all.inside = TRUE)])
title(main = "Reszty w modelu MNK")
legend("bottomleft", legend = c("<mean-sd", "(mean-sd, mean)", "(mean, mean+sd)", ">mean+sd"),
       fill = cols, bty = "n")


# =============================================================================
# Scatterplot dla Morana
# =============================================================================
png(file.path(output_dir, "scatterplot.png"), width = 800, height = 600)
# Reszty modelu OLS
res <- residuals(best_model)

# Wyliczenie reszt sąsiednich (lagged residuals)
lag_res <- lag.listw(cont.listw_subset, res)

# Scatterplot Morana
plot(
  res, lag_res,
  xlab = "res", 
  ylab = "spatially lagged res",
  main = "Wykres rozproszenia testu Morana dla modelu OLS",
  pch = 1,                # Typ pustych okręgów
  col = "black",          # Kolor okręgów
  cex = 1.2               # Wielkość punktów
)

# Dodanie linii regresji
abline(h = 0, v = 0, col = "grey", lty = 2) # Linie poziome i pionowe
lm_fit <- lm(lag_res ~ res) # Dopasowanie regresji
abline(lm_fit, col = "red", lwd = 2) # Linia regresji

# Dodanie etykiet dla punktów odstających
outliers <- which(abs(res) > 0.4 | abs(lag_res) > 0.4) # Progi dla odstających punktów
text(
  res[outliers], lag_res[outliers], 
  labels = names(outliers), 
  pos = 4, cex = 0.8, col = "black"
)
dev.off()


##########################
# Przyjęte wcześniej założenia:
# - Dane: "dane" zawiera zmienne, m.in. log_y, log_x2, log_x5, log_x10 i log_x13.
# - Macierz wag przestrzennych: "cont.listw" (dla pełnego zbioru jednostek) została poprawnie utworzona.
##########################

# Definicja równania modelu – wykorzystujemy najlepszy model z regresji krokowej:
eq <- log_y ~ log_x2 + log_x5 + log_x10 + log_x13

# =============================================================================
# Estymacja modeli przestrzennych
# =============================================================================

# 1. Manski model (full specification) – zawiera przestrzenny lag zmiennej zależnej (rho), przestrzenny lag zmiennych objaśniających (theta) oraz błąd przestrzenny (lambda)
GNS_1 <- sacsarlm(eq, data = dane, listw = cont.listw, type = "sacmixed", method = "LU")
summary(GNS_1)

# 2. SAC / SARAR model – model zawierający przestrzenny lag Y oraz błąd przestrzenny (bez lagów X)
SAC_1 <- sacsarlm(eq, data = dane, listw = cont.listw)
summary(SAC_1)

# 3. SDEM – Spatial Durbin Error Model: model błędu przestrzennego, w którym dodatkowo aktywowana jest opcja etype="emixed" (przestrzenne lagowanie zmiennych X)
SDEM_1 <- errorsarlm(eq, data = dane, listw = cont.listw, etype = "emixed")
summary(SDEM_1)

# 4. SEM – Spatial Error Model: model z przestrzennym błędem (bez przestrzennych lagów X)
SEM_1 <- errorsarlm(eq, data = dane, listw = cont.listw)
summary(SEM_1)

# 5. SDM – Spatial Durbin Model: model z przestrzennym lagiem Y i z przestrzennie opóźnionymi zmiennymi X.
SDM_1 <- lagsarlm(eq, data = dane, listw = cont.listw, type = "mixed")
summary(SDM_1)

# 6. SAR – Spatial Lag Model: model z przestrzennym lagiem Y (bez przestrzennych lagów X).
SAR_1 <- lagsarlm(eq, data = dane, listw = cont.listw)
summary(SAR_1)

# 7. SLX – Model SLX: klasyczny model OLS rozszerzony o przestrzennie opóźnione zmienne objaśniające
SLX_1 <- lmSLX(eq, data = dane, listw = cont.listw)
summary(SLX_1)

# 8. OLS – Klasyczny model regresji liniowej bez składników przestrzennych
OLS_1 <- lm(eq, data = dane)
summary(OLS_1)

# =============================================================================
# Porównanie modeli przy użyciu pakietu texreg
# =============================================================================

screenreg(list(GNS_1, SAC_1, SDEM_1, SEM_1, SDM_1, SAR_1, SLX_1, OLS_1), 
          custom.model.names = c("GNS_1", "SAC_1", "SDEM_1", "SEM_1", "SDM_1", "SAR_1", "SLX_1", "OLS_1"))

# =============================================================================
# Dodatkowa analiza – LR testy (przykłady) oraz estymacja efektów (impacts)
# =============================================================================

LR.Sarlm(GNS_1, SDM_1)
LR.Sarlm(GNS_1, SDEM_1)
LR.Sarlm(GNS_1, SLX_1)
LR.Sarlm(SDM_1, SAR_1)
LR.Sarlm(SDM_1, SEM_1)
LR.Sarlm(SDM_1, SLX_1)

# =============================================================================
#Dodatkowe testy - one chyba powinny być przed modelami
# =============================================================================

eq <- lm(log_y ~ log_x2 + log_x5 + log_x10 + log_x13, data = dane)

# Przeprowadzenie wszystkich testów naraz
test_results <- lm.RStests(eq, cont.listw_subset, test = "all")
print(test_results)

# Wynik dla konkretnego testu

# Dla "LMerr":
test_LMerr <- test_results$LMerr
print(test_LMerr)

# Dla "LMlag":
test_LMlag <- test_results$LMlag
print(test_LMlag)

# Dla SARMA:
test_SARMA <- test_results$SARMA
print(test_SARMA)

#WYNIKI: 

#RSerr (LM test for spatial error model)
#Istotność statystyczna wskazuje na obecność przestrzennego autokorelacji błędu, sugerując, że model może być ulepszony przez dodanie 
#składnika błędu przestrzennego (Spatial Error Model - SEM).

#RSlag (LM test for spatial lag model)
#Wysoka wartość p-wartość wskazuje na brak konieczności rozważenia przestrzennego efektu opóźnienia (Spatial Lag Model - SLM) - SAR

#adjRSerr (adjusted LM test for spatial error model)
#Statystyka ponownie istotna statystycznie co sugeruje, że nawet po dostosowaniu, model błędu przestrzennego pozostaje właściwy do rozważenia.

#adjRSlag (adjusted LM test for spatial lag model)
#Statystyka nie jest istotna, co sugeruje, że po dostosowaniu model z przestrzennym opóźnieniem może nie być niezbędny.

#SARMA (test for spatial autoregressive moving average model)
#wskazuje na możliwość zastosowania modelu SARMA, który łączy zarówno efekty opóźnienia, jak i błędu przestrzennego, 
#dostarczając potencjalnie najpełniejszego modelowania przestrzennych zależności w danych.

# =============================================================================
# najlepszy model SEM 
# =============================================================================

# =============================================================================
output_dir <- "C:\\Users\\ewaru\\Desktop\\e_p"
#ols
png(file.path(output_dir, "przestrzenny_rozklad_reszt_MNK.png"), width = 800, height = 600)
# Wyznaczenie reszt
res <- OLS_1$residuals

# Wyznaczenie zakresu wartości reszt
min_res <- min(res, na.rm = TRUE)
max_res <- max(res, na.rm = TRUE)

# Definiowanie liczby przedziałów
num_classes <- 5  # Możesz zmienić liczbę klas, jeśli chcesz

# Tworzenie przedziałów
brks <- seq(min_res, max_res, length.out = num_classes + 1)

# Definiowanie kolorów dla gradientu
cols <- colorRampPalette(c("darkgrey", "white", "darkred"))(num_classes)

# Tworzenie mapy
plot(pow, col = cols[findInterval(res, brks, all.inside = TRUE)], main = "Przestrzenny rozkład reszt z modelu OLS")

# Dodanie legendy
legend("bottomleft",
       legend = paste0("[", round(brks[-length(brks)], 2), " - ", round(brks[-1], 2), "]"),
       fill = cols,
       title = "Wartość reszt",
       bty = "n")
dev.off()


#SEM
png(file.path(output_dir, "przestrzenny_rozklad_reszt_SEM.png"), width = 800, height = 600)
# Wyznaczenie reszt
res <- SEM_1$residuals

# Wyznaczenie zakresu wartości reszt
min_res <- min(res, na.rm = TRUE)
max_res <- max(res, na.rm = TRUE)

# Definiowanie liczby przedziałów
num_classes <- 5  # Możesz zmienić liczbę klas, jeśli chcesz

# Tworzenie przedziałów
brks <- seq(min_res, max_res, length.out = num_classes + 1)

# Definiowanie kolorów dla gradientu
cols <- colorRampPalette(c("darkgrey", "white", "darkred"))(num_classes)

# Tworzenie mapy
plot(pow, col = cols[findInterval(res, brks, all.inside = TRUE)], main = "Przestrzenny rozkład reszt z modelu SEM")

# Dodanie legendy
legend("bottomleft",
       legend = paste0("[", round(brks[-length(brks)], 2), " - ", round(brks[-1], 2), "]"),
       fill = cols,
       title = "Wartość reszt",
       bty = "n")
dev.off()



