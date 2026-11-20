# =====================================================================
# 3o TRABALHO PRATICO - HIDROLOGIA ESTOCASTICA (EHR806)
# QUESTAO 10 - Esbocar a ACF de  Xt = -0,5 Xt-1 + Zt - 0,8 Zt-1   (ARMA(1,1))
# Baseado em: Aulas-11/03-ARMA-ARIMA.pdf (tabela: ARMA(1,1) com alfa1 < 0 ->
#             ACF com decrescimento OSCILATORIO apos o lag 1)
#             Aulas-11/02-ARMA.R (serie gerada com for, acf, pacf, arima.sim)
#             Aulas/10/01-AR.R   (ACF teorica com for + plot type="h")
# ATENCAO: a formula numerica do rho(1) do ARMA(1,1) NAO esta nos slides
# (la so tem o comportamento). Ela e a formula classica de Box & Jenkins,
# e eu a confiro com a simulacao e com a funcao ARMAacf do R.
# =====================================================================
library(tseries)
graphics.off()
set.seed(13679821)       # semente do professor (Aulas/09/ProcEstocasticos.r)

alfa1 <- -0.5            # coeficiente do AR  (Xt-1)
beta1 <- -0.8            # coeficiente do MA  (Zt-1)

# --- 1) o modelo pode existir? (estacionario e inversivel) --------------
cat("Estacionario: |alfa1| < 1 ?", abs(alfa1) < 1, "\n")
cat("Inversivel:   |beta1| < 1 ?", abs(beta1) < 1, "\n")

# --- 2) ACF TEORICA (a resposta da questao) -----------------------------
# rho(1) = (1 + alfa1*beta1)(alfa1 + beta1) / (1 + 2*alfa1*beta1 + beta1^2)
# rho(k) = alfa1 * rho(k-1)    para k >= 2   (cada passo multiplica por alfa1)
rho1 <- (1 + alfa1*beta1)*(alfa1 + beta1)/(1 + 2*alfa1*beta1 + beta1^2)
acft <- numeric(21)      # lags 0 a 20
acft[1] <- 1             # rho(0) = 1
acft[2] <- rho1          # rho(1)
for (k in 2:20) {
  acft[k+1] <- alfa1*acft[k]
}
cat("rho(1) =", round(rho1, 4), "\n")
print(round(acft[1:6], 4))   # lags 0 a 5
lag <- 0:20
plot(lag, acft, type = "h", col = "blue", ylim = c(-1, 1),
     main = "ACF teorica ARMA(1,1): alfa1=-0,5 e beta1=-0,8")
abline(h = 0)

# --- 3) CONFERENCIA 1: simulando a serie (como no 02-ARMA.R) ------------
xt <- 0
zt <- 0
zt <- rnorm(1,0,1)
for (i in 2:5001) {
  zt[i] <- rnorm(1,0,1)
  xt[i] <- alfa1*xt[i-1] + zt[i] + beta1*zt[i-1]     # a formula do enunciado
}
xt <- xt[2:5001]
kpss.test(xt, null = "L")        # H0: estacionaria (p alto = ok)
par(mfrow = c(1, 3))
plot(xt[1:200], type = "l", col = "red", main = "ARMA(1,1) (200 primeiros)")
acf(xt, type = "correlation", lag.max = 20, main = "ACF simulada")
pacf(xt, lag.max = 20, main = "PACF simulada")
par(mfrow = c(1, 1))
a <- acf(xt, lag.max = 5, plot = FALSE)
cat("ACF simulada lags 1 a 5:", round(a$acf[2:6], 3), "\n")

# --- 4) CONFERENCIA 2: funcao do R (extra, nao e do professor) ----------
print(round(ARMAacf(ar = alfa1, ma = beta1, lag.max = 5), 4))

# FRASE TECNICA: o processo e um ARMA(1,1) estacionario (|alfa1|<1) e inversivel
# (|beta1|<1); sua ACF teorica tem rho(1) = -0,746 e, para k >= 2, rho(k) =
# alfa1*rho(k-1), isto e, decai exponencialmente alternando de sinal
# (decrescimento oscilatorio apos o lag 1, pois alfa1 < 0).
